###########################################################################
#
# Generalizable Avian Movement and Energetics Model
#--------------------------------------------------------------------------
#
# Created by Sarah Jacobi, Eric Lonsdorf, Kevin Aagaard, Wayne Thogmartin, 
# Vicky Hunt, and Tim Jones in collaboration with the Integrated Waterbird 
# Management and Monitoring Program
#
# Modified:  
Sys.time()
# By: Sarah Jacobi 
# 
# Version edits - converting code to data.table format
#
###########################################################################

# This script gathers the necessary data for the GAME model.

# Initialize settings (workspace, directory, etc.) ------------------------

#Clear the workspace
rm(list=ls())

#Clear the console
cat("\014")

# For use when library path directory is not working 
#(permissioning issue)
.libPaths(.libPaths()[2])

#Include necessary libraries
x = c("compiler","fields","maptools","data.table","plyr","foreign","dplyr",
      "ggplot2","scales","stringr","rstudioapi","animation","profvis")

# If you don't know if you've installed the packages for some of these 
# libraries, run this:
# install.packages(x)

lapply(x, library, character.only=TRUE)
rm(x)

#Set working directory (should be universal)
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path)) 
getwd()

#Enable the just in time compiler to speed up analysis
enableJIT(3)

#Grab current time
timeit=proc.time()

#Set the significant digits to 15
options(digits=16)


#
# Define global landscape/movement parameters -----------------------------
# How far birds can move within nodes, in meters
node_move = 3000 
# How fast birds can fly, km/hr
flight_vel = 75
# Cost of flight per hour, kJ; 39,700 kJ per kg
flight_cost_hr = ((226.4 * 1000) / 39700)
# Cost of flight per km, kJ; 39,700 kJ per kg
flight_cost_km = ((3.1 * 1000) / 39700)
# Number of days in migration
fall_start = as.Date("2014-07-01", tz="UTC")
spring_end = strptime("2015-05-31", format="%Y-%m-%d", tz="UTC")
num_days = as.numeric(difftime(as.POSIXct(spring_end), 
                               as.POSIXct(fall_start, tz="UTC"), 
                               units="days")) + 1

# Set day
days=c(1:num_days)

#Threshold for weather severity
wsi_cutoff = 7.5

# Number of land cover classes currently analyzed
num_land_cover_class = 3

# Chose year to pull data
year_from_record = 1960
NonbreedingStart = paste0(year_from_record,"-07-01")

# Function to convert NA-elements to 0
NAFunc = function(x) {x[is.na(x)]<-0; x}

# Function to normalize landscape values
norm_func = function(x) {
  x / sum(x)
}

# Zero-to-one normalization
ZO_norm = function(x) (x-min(x))/(max(x)-min(x))

# Converstion from kilocalories to kiloJoules
kckJConv = 4.184 # 1 kilocalorie = 4.184 kiloJoules

# Energetic parameters
LowCritTemp = 20 #Celsius

flight_dist= 2253.083 #km
range_limit=0.8


#
# Input data --------------------------------------------------------------
# Load node specific data
# #Without WSI data
# NODE_DATA = fread("NodeSpecifData.txt")
# 
# #With WSI data, without Daily Gain
# NODE_DATA = fread("NodeSpecifData_WSI.txt")
# 
# setkey(NODE_DATA,Code,Y_INDEX,X_INDEX,Longitude,Latitude)
# 
# N_0 = NODE_DATA[,sum(N_Abund_0)]
# 
# DailyGain_by_day = fread("NodeDailyGain.txt")
# DailyGain_by_day[,FID:=NULL]
# setkey(DailyGain_by_day,Code,Y_INDEX,X_INDEX,Longitude,Latitude)
# 
# NODE_DATA = NODE_DATA[DailyGain_by_day]
# rm(DailyGain_by_day)
# fwrite(NODE_DATA,file="NodeSpecifData_meanWSI_DG.txt")
# 
#With WSI data and Daily Gain
NODE_DATA = fread("NodeSpecifData_meanWSI_DG.txt")

setkey(NODE_DATA,Code,Y_INDEX,X_INDEX,Longitude,Latitude)

N_0 = NODE_DATA[,sum(N_Abund_0,na.rm=T)]

NODE_DATA[,c(1:18,353),with=F]


#
# # ONLY RUN COMMENTED LINES IF "NodeSpecifData.txt" IS LOST ----------------
{
  # # Load node data and associated files ------------------------------------
  # NODE_DATA = fread("NA_NODE_LC2.txt")
  # 
  # # Organize column headers for easier interpretation
  # setnames(NODE_DATA,names(NODE_DATA),
  #          c("FID","Y_INDEX","X_INDEX","Code","Longitude","Latitude",
  #            "Shoreline","OpenWater","IceSnow","OpenSpace","LowIntDevel",
  #            "MedIntDevel","HighIntDevelop","Barren","DecidForest",
  #            "EvergreenForest","MixedForest","DwarfScrub","ShrubScrub",
  #            "HerbGrass","SedgeHerb","Moss","PastureHay","Crops",
  #            "WoodyWetlands","HerbWetlands","MEAN_DistToCrops",
  #            "SD_DistToCrops","MEAN_WoodyWetlands","SD_WoodyWetlands",
  #            "MEAN_HerbWetlands","SD_HerbWetlands","Prop_82","Prop_90",
  #            "Prop_95","cals","sum","winter","breeding"))
  # 
  # NODE_DATA[,`:=`(Prop_82=NULL,Prop_90=NULL,Prop_95=NULL,cals=NULL,sum=NULL,
  #                 winter=NULL,breeding=NULL,Code=NULL)]
  # 
  # # Load mallard distribution data
  # PopDistribution = data.table(
  #   read.dbf("NorthAmerica_20mi_grid_wAK_BPOP_NSmallard_join.dbf")
  # )
  # 
  # PopDistribution[,`:=`(NatSer_des=NULL,Rowid_1=NULL,ROWID_=NULL)]
  # 
  # # Put Shape data into temporary variable. Note that Natureserve = 3
  # # is non-breeding only; 2 is breeding only
  # PopDistribution_temp =
  #   PopDistribution[,.(ID,Y_INDEX,X_INDEX,NORMPOP,COUNT,Natureserv)]
  # 
  # # Set keys of data table for merge
  # setkey(PopDistribution_temp,Y_INDEX,X_INDEX)
  # setkey(NODE_DATA,Y_INDEX,X_INDEX)
  # 
  # # Combine node data with population distribution
  # NODE_DATA = PopDistribution_temp[NODE_DATA]
  # setnames(NODE_DATA,c("NORMPOP","Natureserv"),
  #          c("NormPop","NatureServCategory"))
  # 
  # rm(PopDistribution,PopDistribution_temp)
  # 
  # # Invert Y_INDEX order to go from small to large
  # NODE_DATA[,`:=`(Y_INDEX=max(Y_INDEX)+1-Y_INDEX,
  #                 X_INDEX=X_INDEX-(min(X_INDEX)-1))]
  # 
  # # Coordinatess from Node data
  # NODE_DATA[,Code:= paste0(NODE_DATA[,X_INDEX],"_",NODE_DATA[,Y_INDEX])]
  # 
  # # Determine number of nodes
  num_nodes = nrow(NODE_DATA)
  # #
  # # Set node size in km (20 miles = 32.1869 km)
  node_size = 32.1869
  # 
  # # Identify which nodes are in the breeding area
  # breeding_nodes = ifelse(NODE_DATA[,NatureServCategory==2],
  #                         NODE_DATA[,NormPop],0)
  # breeding_nodes[is.na(breeding_nodes)] = 0
  # 
  # # Initial population size
  # N_0 = 1000000*sum(breeding_nodes)
  # 
  # 
  # #
  # # Determine distribution of birds and forage/roosting data ----------------
  # # Set up matrix to hold habitat data
  # hab_dist=matrix(0,num_nodes,num_land_cover_class)
  # 
  # # Incorporate distance effects. Create normal distribution using mean and
  # # stdev dist to cover
  # hab_dist[,1] = pnorm(node_move,NODE_DATA[,MEAN_DistToCrops],
  #                      sapply(NODE_DATA[,SD_DistToCrops],
  #                             function(x) max(x,30)))
  # hab_dist[,2] = pnorm(node_move,NODE_DATA[,MEAN_WoodyWetlands],
  #                      sapply(NODE_DATA[,SD_WoodyWetlands],
  #                             function(x) max(x,30)))
  # hab_dist[,3] = pnorm(node_move,NODE_DATA[,MEAN_HerbWetlands],
  #                      sapply(NODE_DATA[,SD_HerbWetlands],
  #                             function(x) max(x,30)))
  # hab_dist[which(NODE_DATA[,Crops]==0),1] = 0
  # hab_dist[which(NODE_DATA[,WoodyWetlands]==0),2] = 0
  # hab_dist[which(NODE_DATA[,HerbWetlands]==0),3] = 0
  # 
  # # Convert bird habitat data from m^2 cm^2
  # bird_hab =
  #   cbind(NODE_DATA[,.(FID,Y_INDEX,X_INDEX,Code,Longitude,
  #                      Latitude,
  #                      Shoreline=as.numeric(Shoreline)/900)],
  #         (NODE_DATA[,.(Crops,WoodyWetlands,HerbWetlands)] *
  #            hab_dist)/900)
  # 
  # # Remove unneeded data
  # rm(hab_dist)
  # 
  # # Distribute the birds among breeding nodes (shoreline + herb/grass)
  # # Use this to keep land cover types separate going forward
  # cals_by_cover = bird_hab[,.(Shoreline,Crops,WoodyWetlands,HerbWetlands)]
  # 
  # # 33945 kcal = 143k kJ, based on energy in breeding habitat
  # NODE_DATA[,cals_sum:=rowSums(cals_by_cover)*33945 ]
  # 
  # # Identify breeding and wintering nodes
  # NODE_DATA[NatureServCategory==2,
  #           `:=`(BPOP_wt=NormPop/sum(NormPop,na.rm=TRUE),
  #                cal_wt=cals_sum/sum(cals_sum,na.rm=TRUE))]
  # 
  # # Calculate normalized calories and weight over breeding nodes
  # NODE_DATA[,breed_wt:=(BPOP_wt*cal_wt)/sum(BPOP_wt*cal_wt,na.rm=TRUE)]
  # 
  # # Distribute the birds based on the breeding weight
  # NODE_DATA[,N_Abund_0:=round(N_0*breed_wt)]
  # nacols = "N_Abund_0"
  # NODE_DATA[,c(nacols):=lapply(.SD,NAFunc),.SDcols=nacols]
  # 
  # # Clean up temp variables
  # rm(cals_by_cover)
  # gc()
  # 
  # # The initial population is distributed among the breeding_nodes
  # N_0 = NODE_DATA[NatureServCategory==2,sum(N_Abund_0)]
  # 
  # # DISTANCE TO BREEDING AREA
  # # Distance from each node to each breeding node
  # distance_to_breed_nodes =
  #   rdist(NODE_DATA[,.(X_INDEX,Y_INDEX)],
  #         NODE_DATA[NatureServCategory==2,.(X_INDEX,Y_INDEX)]) *
  #   node_size
  # 
  # # Determine distance to nearest breeding node
  # NODE_DATA[,DistBreed :=
  #             as.matrix(apply(distance_to_breed_nodes,1,min))]
  # 
  # # Garbage collector
  # rm(distance_to_breed_nodes)
  # gc()
  # 
  # # ROOSTING HABITAT
  # # Percent of total habitat that is shoreline (dabbler roosting habitat)
  # target_sums = c("Shoreline","OpenWater","IceSnow","OpenSpace",
  #                 "LowIntDevel","MedIntDevel","HighIntDevelop","Barren",
  #                 "DecidForest","EvergreenForest","MixedForest",
  #                 "DwarfScrub","ShrubScrub","HerbGrass","SedgeHerb",
  #                 "Moss","PastureHay","Crops","WoodyWetlands",
  #                 "HerbWetlands")
  # NODE_DATA[,RoostShoreline := (Shoreline / rowSums(.SD,na.rm=T)),
  #           .SDcols=target_sums]
  # 
  # # HUMAN DISTURBANCE EXTENT
  # # Disturbance will affect the movement/expenditures to forage
  # # May affect leave/stay decisions
  # # Disturbance will decrease daily gain, maybe enough force departure
  # 
  # # Create different distrubance weights for each cover type (open space,
  # # low intensity, medium intensity, high intensity development)
  # disturb_weights = c(0.4, 0.6, 0.8, 1)
  # 
  # # Disturbance is calculated as sum of open space, low/med/high intensity
  # # development divided by total habitat weight by disturbance weights
  # disturb_mtx = as.matrix(NODE_DATA[,.(OpenSpace,LowIntDevel,MedIntDevel,
  #                                      HighIntDevelop)])
  # 
  # NODE_DATA[,HumanDisturb := disturb_mtx %*% (disturb_weights) /
  #             rowSums(.SD,na.rm=T),
  #           .SDcols=target_sums]
  # 
  # # FORAGE AVAILABILITY
  # # Parameter scenarios
  # # c1 param ID; c2 shoreline cals; c3 crop cals; c4 herb wet
  # # cals; c5 woody wet cals; c6 tank cals; c7 flight cost; c8 flight speed;
  # # c9 flight fuel efficiency (c8/c7); c10 flight range (c9*c6)
  # 
  # param_scenarios = fread("param_scenarios_5_8_14.txt")
  # gc()
  # 
  # # Parameter row
  # rr=1
  # 
  # # Keep cals separate by cover type
  # cals.Shoreline = as.matrix(as.numeric(bird_hab[,Shoreline]) *
  #                              as.matrix(param_scenarios[rr,shoreline]))
  # cals.Crops = as.matrix(as.numeric(bird_hab[,Crops]) *
  #                          as.matrix(param_scenarios[rr,crops]))
  # cals.WoodyWetlands =
  #   as.matrix(as.numeric(bird_hab[,WoodyWetlands]) *
  #               as.matrix(param_scenarios[rr,woody_wet]))
  # cals.HerbWetlands =
  #   as.matrix(as.numeric(bird_hab[,HerbWetlands]) *
  #               as.matrix(param_scenarios[rr,herbaceous_wet]))
  # 
  # # Add node id data to input data
  # NODE_DATA[,`:=`(ForageShoreline=cals.Shoreline,
  #                 ForageCrops=cals.Crops,
  #                 ForageWoodyWetlands=cals.WoodyWetlands,
  #                 ForageHerbWetlands=cals.HerbWetlands)]
  # 
  # 
  # #
  # # Clean up table containing node-specific data ----------------------------
  # NODE_DATA[,`:=`(
  #   FID=NULL,Shoreline=NULL,OpenWater=NULL,IceSnow=NULL,OpenSpace=NULL,
  #   LowIntDevel=NULL,MedIntDevel=NULL,HighIntDevelop=NULL,Barren=NULL,
  #   DecidForest=NULL,EvergreenForest=NULL,MixedForest=NULL,DwarfScrub=NULL,
  #   ShrubScrub=NULL,HerbGrass=NULL,SedgeHerb=NULL,Moss=NULL,PastureHay=NULL,
  #   Crops=NULL,WoodyWetlands=NULL,HerbWetlands=NULL,MEAN_DistToCrops=NULL,
  #   SD_DistToCrops=NULL,MEAN_WoodyWetlands=NULL,SD_WoodyWetlands=NULL,
  #   MEAN_HerbWetlands=NULL,SD_HerbWetlands=NULL,cals_sum=NULL,BPOP_wt=NULL,
  #   cal_wt=NULL,breed_wt=NULL,NatureServCategory=NULL
  #   )]
  # 
  # 
  # #
  # # Clean up environment ----------------------------------------------------
  # rm(bird_hab,cals.Crops,cals.HerbWetlands,cals.Shoreline,cals.WoodyWetlands,
  #    disturb_mtx,disturb_weights,wsi_cutoff,target_sums,spring_end,
  #    fall_start,breeding_nodes,nacols,num_land_cover_class)
  # gc()
  # 
  # 
  # #
  # # Save node specific data table -------------------------------------------
  # fwrite(NODE_DATA,file="NodeSpecifData.txt")
}

# # RUN THESE LINES SEPARATELY FROM THE ABOVE SECTION, AS THEY REQUIRE
# # ALL AVAILABLE MEMORY
{
  # # Calculate cummulative distance-dependent movement probabilities -------
  # # Determine probability of going from one node to each other node based
  # # on gamma distribution
  # NODE_DATA_DISTANCE = as.matrix(NODE_DATA[,.(X_INDEX,Y_INDEX)])
  # rownames(NODE_DATA_DISTANCE) = NODE_DATA[,Code]
  # NODE_DATA_DISTANCE = rdist(NODE_DATA_DISTANCE) * node_size
  # gamma_b = flight_dist-flight_dist*range_limit
  # gamma_x = flight_dist - NODE_DATA_DISTANCE
  # 
  # rm(NODE_DATA_DISTANCE)
  # gc()
  # 
  # prelim_gamma_dist_prob = ((exp(1) ^ (- (gamma_x / gamma_b)) /
  #                              gamma_b)) * (gamma_x < 0)
  # 
  # rm(gamma_x)
  # gc()
  # 
  # gamma_dist_prob = prelim_gamma_dist_prob / sum(prelim_gamma_dist_prob)
  # 
  # rm(prelim_gamma_dist_prob)
  # gc()
  # 
  # NODE_DATA[,GammaMoveProb:=rowSums(gamma_dist_prob)]
  # rm(gamma_dist_prob)
  # gc()
  # 
  # 
  # #
  # # Save node specific data table -----------------------------------------
  # # fwrite(NODE_DATA,file="NodeSpecifData.txt")
}

# # RUN THESE LINES TO RECREATE THE NODE-SPECIFIC WSI DATA
{
  # WSI data ----------------------------------------------------------------
  # # Choose a given year from the historical record
  # #WSI_DATA = fread("NorAm_historical_wsi_node_info.txt")
  # 
  # # Use the average of the record
  # WSI_DATA[,`:=`(Y_INDEX=max(Y_INDEX)+1-Y_INDEX,
  #                X_INDEX=X_INDEX-(min(X_INDEX)-1))]
  # WSI_DATA_frame = data.frame(WSI_DATA)
  # data_columns = colnames(WSI_DATA_frame)
  # 
  # # Subset out a particular migration period (don't want to base this on
  # # calendar year, but rather on annual migratory cycle; i.e., fall one year
  # # to spring the next)
  # #
  # # Use this with historical data:
  # # 1972 was the most recent "average" (relative to 1951-1980) year,
  # # differing from the mean during that time span by 0.01 deg C.
  # # target_year = year_from_record
  # # nonbreeding_year = c(as.character(target_year),
  # #                      as.character(target_year+1))
  # # target_dates=c(paste0(nonbreeding_year[1],c("_07_","_08_","_09_","_10_",
  # #                                             "_11_","_12_")),
  # #                paste0(nonbreeding_year[2],c("_01_","_02_","_03_","_04_",
  # #                                             "_05_")))
  # # 
  # # Use this with mean data:
  # WSI_DATA = fread("NorAm_mean_wsi_node_info.txt")
  # WSI_DATA[,"02_29":=NULL]
  # WSI_DATA[,`:=`(Y_INDEX=max(Y_INDEX)+1-Y_INDEX,
  #                X_INDEX=X_INDEX-(min(X_INDEX)-1))]
  # WSI_DATA_frame = data.frame(WSI_DATA)
  # data_columns = colnames(WSI_DATA_frame)
  # 
  # target_dates=c("07_","08_","09_","10_","11_","12_","01_","02_",
  #                "03_","04_","05_")
  # 
  # # Use this with either dataset:
  # nonbreeding_period = data_columns[grepl(target_dates[1],data_columns)]
  # for(i in 2:length(target_dates)){
  #   nonbreeding_period =
  #     c(nonbreeding_period,data_columns[grepl(target_dates[i],
  #                                             data_columns)])
  # }
  # 
  # WSI_DATA.nonbreeding_period =
  #   data.table(WSI_DATA_frame[,c('X_INDEX','Y_INDEX',nonbreeding_period)])
  # WSI_DATA.nonbreeding_period[,Code:=paste0(
  #   WSI_DATA.nonbreeding_period[,X_INDEX],"_",
  #   WSI_DATA.nonbreeding_period[,Y_INDEX])]
  # setcolorder(WSI_DATA.nonbreeding_period,c(num_days+3,1:(num_days+2)))
  # setnames(WSI_DATA.nonbreeding_period,names(WSI_DATA.nonbreeding_period),
  #          c("Code","X_INDEX","Y_INDEX",paste0("WSI_day_",1:num_days)))
  # 
  # setkey(WSI_DATA.nonbreeding_period,Code,X_INDEX,Y_INDEX)
  # 
  # # ggplot(data=WSI_DATA.nonbreeding_period,
  # #        aes(x=X_INDEX,y=Y_INDEX,col=log(WSI_day_1))) +
  # #   geom_point(size=1.75,shape=15) +
  # #   scale_color_gradientn(colors=tim.colors(10)) +
  # #   theme_minimal() +
  # #   theme(
  # #     axis.text=element_blank(),
  # #     axis.title=element_blank(),
  # #     panel.grid.major=element_blank(),
  # #     panel.grid.minor=element_blank()
  # #   )
  # 
  # setkey(NODE_DATA,Code,X_INDEX,Y_INDEX)
  # NODE_DATA = NODE_DATA[WSI_DATA.nonbreeding_period]
  # 
  # rm(WSI_DATA,WSI_DATA.nonbreeding_period,WSI_DATA_frame,data_columns)
  # 
  # fwrite(NODE_DATA,file="NodeSpecifData_WSI.txt")
}


# # RUN THESE LINES TO RECREATE THE NODE-SPECIFIC DAILY GAIN FUNCTION
{
  # Calculate node daily gain -----------------------------------------------
  # # THIS ONLY NEEDS TO BE RUN TO IDENTIFY THE RELATION; ONCE FOUND, NO NEED
  # # TO UNCOMMENT THIS
  # #
  # # Calculate factor by which to multiply base 1 gram Temperature Dependent
  # # BMR by. That is, to calculate the temperature-dependent daily gain for
  # # a 1.14 kilogram bird, multiply the TDBMR by the output of
  # # `TD_DG_BodyMass_Relation(1.14)`.
  #
  # # Calculate temperature dependent BMR for 1 gram
  # BodyMass = 0.001
  # # Energy loss functions
  # # BMR in kJ, Heitmeyer 2010
  # BasalMetRate = 422 * (BodyMass^0.74)
  #
  # # BMR below LCT, from Prince 1979
  # BMR_LCT = function(temperature) {
  #   ((2.06 - (0.03 * temperature)) * kckJConv) * BasalMetRate
  # }
  #
  # TempDepMetRate =
  #   function(temperature) {
  #     (BMR_LCT(temperature) ^ (temperature<LowCritTemp)) *
  #       (BasalMetRate^(temperature>=LowCritTemp))
  #   }
  #
  # BodyMass = c(0.001,seq(round(range(body_mass)[1]/1000,1),
  #                        round(range(body_mass)[2]/1000,1),(1/1000)))
  # BasalMetRate = 422 * (BodyMass^0.74)
  #
  # temperaturerange=c(-60:31)
  #
  # TDBMR=vector("list")
  # for(i in 1:length(BodyMass)){
  #   TDBMR[[i]] =
  #     (((2.06 - (0.03 * temperaturerange)) * kckJConv) *
  #        BasalMetRate[i]) ^ (temperaturerange<LowCritTemp) *
  #     (BasalMetRate[i] ^ (temperaturerange>=LowCritTemp))
  # }
  #
  # # Energy gain function
  # # Energy gain functions; 9 kcal/g, Ricklefs 1974 via Whyte and Bolen 1988
  # # MeanDailyMassIncrease = 4 #g, Hanson et al. 1990
  # # RWB_coeff = MeanDailyMassIncrease * (9 * kckJConv)
  # #
  # # LK_coeff = 4.6              #Lindstrom and Kvist 1995; passerines
  # # KL_coeff = runif(1,3,10)    #Kvist and Lindstrom 2003; waders, ~3-10
  #
  # Empir_coeff = 1.8           #range of MetCostLandscape = 1:~2.2XBMR
  #
  # # DailyMetabolizableEnergy = RWB_coeff + MetCostLandscape
  # # DailyMetabolizableEnergy = LK_coeff * BasalMetRate
  # # DailyMetabolizableEnergy = KL_coeff * BasalMetRate
  # DailyMetabolizableEnergy = Empir_coeff * BasalMetRate
  #
  # HarvestDisturb = 0
  # FDR = DailyMetabolizableEnergy
  # HumanDisturb = 0.5
  # TD_DG=vector("list")
  # for(i in 1:length(BodyMass)){
  #   TD_DG[[i]] = ((1-HumanDisturb) * (1-HarvestDisturb) * FDR[i]) -
  #     TDBMR[[i]]
  # }
  #
  # TD_DG_Change = vector()
  # for(i in 1:(length(TD_DG)-1)) {
  #   TD_DG_Change[i] = (TD_DG[[i+1]] / TD_DG[[1]])[1]
  # }
  #
  # TD_DG_BodyMass_Relation_eq=
  #   nls(TD_DG_Change~a*BodyMass[2:length(BodyMass)]^2 +
  #         b*BodyMass[2:length(BodyMass)] + c,
  #       start=list(a=0,b=0,c=1))
  #
  # plot(TD_DG_Change~BodyMass[2:length(BodyMass)])
  # summary(TD_DG_BodyMass_Relation_eq)$parameters
  #
  # # TD DG Relation = (-14.80800820381551 * BodyMass^2 +
  # #                     153.24219961825537 * BodyMass +
  # #                     27.45433867668212)
  #
  # TD_DG_BodyMass_Relation = function(x) {
  #   (summary(TD_DG_BodyMass_Relation_eq)$parameters[1] * (x ^ 2)) +
  #     (summary(TD_DG_BodyMass_Relation_eq)$parameters[2] * x) +
  #     summary(TD_DG_BodyMass_Relation_eq)$parameters[3]
  # }
  #
  # TD_DG_BodyMass_Relation = function(x) {
  #   (-14.80800820381551 * (x ^ 2)) +
  #     (153.24219961825537 * x) +
  #     27.45433867668212
  # }
  # 
  # # Daily gain is calculated with the inputs:
  # #   HumanDisturb (Degree of human disturbance--urbanization)
  # #   HarvestDisturb (Degree of human disturbance--hunter harvest)
  # #   FDR (Fuel deposition rate, kJ ingested per day)
  # #   TD_BMR (Temperature dependent basal metabolic rate, kJ per day)
  # 
  # # DG = (HumanDisturb * HarvestDisturb * FDR) - TD_BMR
  # 
  # # Read in temperature data per node
  # TempData = fread("NorAm_historical_air_temperature_node_info.txt")
  # TempData[,`:=`(Y_INDEX=max(Y_INDEX)+1-Y_INDEX,
  #                X_INDEX=X_INDEX-(min(X_INDEX)-1))]
  # TempDataFrame = data.frame(TempData)
  # data_columns = colnames(TempDataFrame)
  # 
  # #subset out a particular migration period (don't want to base this on
  # #calendar year, but rather on annual migratory cycle; i.e., fall one year
  # #to spring the next)
  # # 1972 was the most recent "average" (relative to 1951-1980) year,
  # # differing from the mean during that time span by 0.01 deg C.
  # target_year = 1972
  # nonbreeding_year = c(as.character(target_year),
  #                      as.character(target_year+1))
  # target_dates=c(paste0(nonbreeding_year[1],c("_07_","_08_","_09_","_10_",
  #                                             "_11_","_12_")),
  #                paste0(nonbreeding_year[2],c("_01_","_02_","_03_","_04_",
  #                                             "_05_")))
  # 
  # nonbreeding_period = data_columns[grepl(target_dates[1],data_columns)]
  # for(i in 2:length(target_dates)){
  #   nonbreeding_period =
  #     c(nonbreeding_period,data_columns[grepl(target_dates[i],
  #                                             data_columns)])
  # }
  # 
  # TempDT_NonbreedingPeriod =
  #   data.table(TempDataFrame[,c('X_INDEX','Y_INDEX',nonbreeding_period)])
  # TempDT_NonbreedingPeriod[,Code:=paste0(
  #   TempDT_NonbreedingPeriod[,X_INDEX],"_",
  #   TempDT_NonbreedingPeriod[,Y_INDEX])]
  # TempDT_NonbreedingPeriod[,`:=`(X_INDEX=NULL,Y_INDEX=NULL)]
  # setcolorder(TempDT_NonbreedingPeriod,c(num_days+1,1:num_days))
  # setnames(TempDT_NonbreedingPeriod,names(TempDT_NonbreedingPeriod),
  #          c("Code",paste0("TD_BMR_day_",1:num_days)))
  # 
  # # Calculate temperature dependent BMR by node by day
  # BodyMass = 0.001
  # 
  # # Energy gain function
  # # Energy gain functions; 9 kcal/g, Ricklefs 1974 via Whyte and Bolen 1988
  # # MeanDailyMassIncrease = 4 #g, Hanson et al. 1990
  # # RWB_coeff = MeanDailyMassIncrease * (9 * kckJConv)
  # #
  # # LK_coeff = 4.6              #Lindstrom and Kvist 1995; passerines
  # # KL_coeff = runif(1,3,10)    #Kvist and Lindstrom 2003; waders, ~3-10
  # 
  # Empir_coeff = 1.8           #range of MetCostLandscape = 1:~2.2XBMR
  # 
  # # DailyMetabolizableEnergy = RWB_coeff + MetCostLandscape
  # # DailyMetabolizableEnergy = LK_coeff * BasalMetRate
  # # DailyMetabolizableEnergy = KL_coeff * BasalMetRate
  # DailyMetabolizableEnergy = Empir_coeff * BasalMetRate
  # 
  # # Metabolic landscape
  # TargetCols = names(TempDT_NonbreedingPeriod)[
  #   which(names(TempDT_NonbreedingPeriod) %like% "TD_BMR_day_")]
  # MetCostDT_NonbreedingPeriod =
  #   cbind(Code=TempDT_NonbreedingPeriod[,Code],
  #         TempDT_NonbreedingPeriod[,lapply(.SD,FUN=TempDepMetRate),
  #                                  .SDcols=TargetCols])
  # setkey(MetCostDT_NonbreedingPeriod,Code)
  # Mean_MetCost =
  #   round(mean(
  #     as.matrix(MetCostDT_NonbreedingPeriod[,(TargetCols),with=F])
  #   ),2)
  # 
  # # Calculate daily gain in terms of grams
  # setkey(NODE_DATA,Code)
  # NODE_DATA = NODE_DATA[MetCostDT_NonbreedingPeriod]
  # 
  # NODE_DATA[,`:=`(HarvestDisturb = 0,
  #                 FDR = DailyMetabolizableEnergy)]
  # 
  # InCols = names(NODE_DATA)[which(names(NODE_DATA) %like% "TD_BMR_day_")]
  # OutCols = paste0("DailyGain_day",1:num_days)
  # NODE_DATA[,c(OutCols) :=
  #             ((1-HumanDisturb) * (1-HarvestDisturb) * FDR) - .SD,
  #           .SDcols=InCols]
  # 
  # DailyGain_by_day = cbind(NODE_DATA[,.(FID,Y_INDEX,X_INDEX,Code,Longitude,
  #                                       Latitude)],
  #                          NODE_DATA[,which(names(NODE_DATA) %like%
  #                                             "DailyGain"),with=F])
  # 
  # Rm_Cols = names(NODE_DATA)[which(names(NODE_DATA) %like% "TD_BMR_day_")]
  # NODE_DATA[,c(Rm_Cols):=NULL]
  # 
  # # fwrite(DailyGain_by_day,file="NodeDailyGain.txt")
}


#
# Generate body mass/body condition classes -------------------------------
# Average body mass, in grams
mean_body_mass = 1140
min_body_mass = 800
max_body_mass = 1600
beta_mean_body_mass = mean(c(mean_body_mass,max_body_mass))

a=10
b=1.5

body_mas_pop = N_0
#body_mas_pop = 10000
body_mass = rbeta(body_mas_pop,a,b) * beta_mean_body_mass
range(body_mass)
body_mass = body_mass[body_mass > min_body_mass]

# plot(density(body_mass))

bc_bins = c(1:20)

body_mass_optim = density(body_mass)$x[which.max(density(body_mass)$y)]
# body_mass_discrete = c(seq(range(body_mass)[1],body_mass_optim,
#                            length.out=length(bc_bins)),
#                        range(body_mass)[2])
body_mass_discrete = seq(range(body_mass)[1],range(body_mass)[2],
                         length.out=length(bc_bins)+1)
body_mass_norm = ZO_norm(body_mass_discrete)
body_mass_lower = body_mass_norm[1:length(body_mass_norm)-1]
body_mass_upper = body_mass_norm[2:length(body_mass_norm)]

integrate_beta = function(l,h) {
  integrate(Vectorize(function(x) dbeta(x,a,b)),
            lower=l,upper=h)$value
}

body_mass_dist = sapply(body_mass_upper,integrate_beta,l=body_mass_lower)

BodyCondition_Table = 
  data.table(BC_bins=bc_bins,
             BodyMass_g=body_mass_discrete[-1],
             BodyMass_dist_discrete_day0=body_mass_dist - 
               data.table::shift(body_mass_dist,n=1L,fill=0,type="lag"),
             BodyFat_g=
               body_mass_discrete[-1]*seq(0,0.2,
                                          length.out=length(bc_bins)))

#Parameter for gamma distribution to determine movement probability
BodyCondition_Table[,BodyFat_gamma_b:=(flight_dist * ZO_norm(BodyFat_g))]

#Bottom bin has 99.1% daily survivorship
min_survival = .9975 #0.991
#Top bin has 99.997% daily survivorship
max_survival = .99997

daily_survival = vector()
for (i in 1:length(bc_bins)) {
  daily_survival[i] = ((((max_survival - min_survival) / 
                           length(bc_bins)) * i) 
                       + min_survival)
}

BodyCondition_Table[,Survivorship:=daily_survival]

#Possible flight distance by body condition (in km)
BodyCondition_Table[,MaxLike_flightdist:=BodyFat_g/flight_cost_km]
BodyCondition_Table[20,MaxLike_flightdist:=flight_dist]

#Calculate Pr(flying distance, x) based on body condition
g=0
j=flight_dist

integrate_unif = function(l,h) {
  integrate(Vectorize(function(x) dunif(x,g,j)),
            lower=l,upper=h)$value
}

flightdist_distrib = sapply(BodyCondition_Table[,MaxLike_flightdist],
                            integrate_unif,l=0)

# Body condition-dependent flight distance probabilities
BC_distprob = 1/BodyCondition_Table[,MaxLike_flightdist]

BC_dist_seq_ls = list()
for(i in 1:length(bc_bins)){
  BC_dist_seq_ls[[i]] = rep(BC_distprob[i],
                            (((flightdist_distrib[i]*flight_dist)+1))/node_size)
}

BC_dist_seq_ls_maxln = max(sapply(BC_dist_seq_ls,length))
BC_dist_seq_mtx = 
  do.call(cbind,lapply(BC_dist_seq_ls,
                       function(x){
                         length(x) = BC_dist_seq_ls_maxln
                         x
                       }
  ))

BC_dist_seq_mtx = NAFunc(BC_dist_seq_mtx)

# Calculate the cost of traveling each unit (node) of a flight
cummul_flightcost = 
  data.table(dist_km = seq(0,flight_dist,node_size),
    cost_km=c(0,cumsum(rep((flight_cost_km*node_size),(flight_dist/node_size)))))

# Calculate the number of bins by which birds will be reduced for flying
# different distances
BC_transition_g_ls = list()
for(i in length(bc_bins):2){
  BC_transition_g_ls[[i]] = BodyCondition_Table[i,BodyFat_g] - 
    BodyCondition_Table[1:(i-1),BodyFat_g]
}

BC_transition_g_ls_maxln = max(sapply(BC_transition_g_ls,length))
BC_transition_g_mtx = 
  do.call(cbind,lapply(BC_transition_g_ls,
                       function(x){
                         length(x) = BC_transition_g_ls_maxln
                         x
                       }
  ))

BC_transition_g_mtx = NAFunc(BC_transition_g_mtx)
BC_transition_g_mtx = rbind(BC_transition_g_mtx,rep(0,length(bc_bins)-1))
BC_transition_g_mtx = cbind(rep(0,length(bc_bins)),BC_transition_g_mtx)

BC_transitions = 
  sort(BC_transition_g_mtx[,20][which(BC_transition_g_mtx[,20]>0)])

BC_decrement = vector()
BC_decrement_total = rep(0,nrow(cummul_flightcost))
for(i in 1:length(BC_transitions)){
  BC_decrement=
    ifelse(cummul_flightcost[,cost_km]<BC_transitions[i],0,1)
  BC_decrement_total = BC_decrement + BC_decrement_total
}

cummul_flightcost[,BC_decrementvalue:=BC_decrement_total]                  


#
# Impose movement rules ---------------------------------------------------
# Weights for WSI, body condition, and distance to breeding node
weight_WSI = 1/3
weight_BC = 1/3
weight_DB = 1/3

# Probability of departure component related to WSI
# The relation below comes from Schummer et al. 2010
WSIProb = function(x) -0.6282+(0.0547*x)+(0.0046*((x)^2))

# Probability of departure component related to body condition
ProbDepart_BC = c(0,(weight_BC * ((bc_bins^8) / 
                                    ((bc_bins^8) + 
                                       ((length(bc_bins)-3)^8))))[-1])
ProbDepart_BC_mtx = matrix(rep(ProbDepart_BC,nrow(NODE_DATA)),
                           ncol=length(bc_bins),byrow=T)

# Probability of departure component related to distance to nearest 
# breeding node
breedingdistances = NODE_DATA[,DistBreed]
n_breedingdistances = max(breedingdistances)
ProbDepart_DB = weight_DB * (breedingdistances ^ 5) / 
  ((breedingdistances ^ 5) + 
     ((n_breedingdistances-(n_breedingdistances)*0.15)) ^ 5)
ProbDepart_DB_mtx = matrix(rep(ProbDepart_DB,length(bc_bins)),
                           ncol=length(bc_bins),byrow=T)

ProbMort = 1 - BodyCondition_Table[,Survivorship]
ProbMort_mtx = matrix(rep(ProbMort,nrow(NODE_DATA)),ncol=length(bc_bins),
                      byrow=T)

# BodyMass_DailyGain_Relation function (see commented code above for
# parameter value generation)
TD_DG_BodyMass_Relation = function(x) {
  (-14.80800820381551 * ((x/1000) ^ 2)) + 
    (153.24219961825537 * (x/1000)) + 
    27.45433867668212
}
TD_DG_BM_conversion = 
  TD_DG_BodyMass_Relation(BodyCondition_Table[,BodyMass_g])

# Natural forage decay rates by landcover type
Shoreline_decay = 0.9998
Crops_decay = 0.997
WoodyWetlands_decay = 0.9965
HerbWetlands_decay = 0.991

# Normalize roosting, distance to breeding
NODE_DATA[,`:=`(
  NormRoostShoreline = RoostShoreline/ sum(RoostShoreline),
  NormDistBreed = (max(DistBreed) - DistBreed) / 
    sum(max(DistBreed) - DistBreed)
)]

# Set values for CD exponents
W_Forage = 0.55
W_Roost = 0.25
W_Breed = 0.8
W_Gamma = 0.1

# Establish inflection point--date on which importance of landscape 
# variables switch directions
inflection = floor(runif(1,184,258))

# Generate sequence of weights
inflection_seq = c(seq(1,0.1,length.out=inflection),
                   seq(0.1,1,length.out=num_days-inflection))

WeightScheme = data.table(
  ForageAvailComponent = ((((-inflection_seq^3) * W_Forage) + 
                             (max(inflection_seq)^3)*W_Forage) + 0.03),
  RoostingHabComponent = ((((-inflection_seq^3) * W_Roost) + 
                             (max(inflection_seq)^3) * W_Roost) + 0.03),
  BreedingDistComponent = (((inflection_seq^3) * W_Breed) + 0.03),
  GammaMoveComponent = W_Gamma + 0.01
)

deadbirds = vector()

Sys.time()
start_timer = Sys.time()
for(i in 1:length(days)) {
  # Reduce forage material based on mass- and temperature-dependent daily 
  # gain for each body mass/condition class and landcover-specific decay 
  # rates, with place holders for decay rates
  DayAbund = 
    as.matrix(NODE_DATA[,paste0("N_Abund_",i-1),with=F])
  DailyGain = 
    as.matrix(NODE_DATA[,paste0("DailyGain_day",i),with=F])
  ReducedForage = paste0("Forage_day_",i)
  
  BM_PopDist = 
    as.vector(as.matrix(
      BodyCondition_Table[,paste0("BodyMass_dist_discrete_day",i-1),
                          with=F]))
  
  NODE_DATA[,(ReducedForage):=(((ForageShoreline * (Shoreline_decay^i)) +
                                  (ForageCrops * (Crops_decay^i)) +
                                  (ForageWoodyWetlands * 
                                     (WoodyWetlands_decay^i)) +
                                  (ForageHerbWetlands * 
                                     (HerbWetlands_decay^i)))) - 
              (rowSums(matrix((DailyGain %o% TD_DG_BM_conversion),
                              ncol=length(bc_bins)) * 
                         matrix((DayAbund %o% BM_PopDist),
                                ncol=length(bc_bins)),na.rm=T))]
  
  NormReducedForage = paste0("NormForage_day_",i)
  incols = ReducedForage
  NODE_DATA[,(NormReducedForage) := lapply(.SD,norm_func),
            .SDcols=incols]
  
  # Calculate kilojoules gained per body condition to redistribute 
  # individuals among body condition classes after foraging
  DG_BC_Node = matrix(((DailyGain/39.7) %o% (TD_DG_BM_conversion/39.7)), 
                      ncol=length(bc_bins))
  
  DG_BC_Transition = BodyCondition_Table[,outer(BodyFat_g,BodyFat_g,'-')]
  
  DG_BC_Node_Transition = vector("list")
  for(j in 1:length(bc_bins)){
    DG_BC_Node_Transition[[j]] = findInterval(DG_BC_Node[,j],
                                              DG_BC_Transition[,j])
  }
  
  DG_BC_Node_Transition_mtx = 
    matrix(unlist(DG_BC_Node_Transition), ncol=length(bc_bins), byrow=F)
  
  BC_Abund = matrix((DayAbund %o% BM_PopDist),ncol=length(bc_bins))
  
  BC_Abund_Gain_func = function(x) {
    sum(BC_Abund[which(DG_BC_Node_Transition_mtx==x)])
  }
  
  BC_Abund_Next = sapply(1:20,BC_Abund_Gain_func)
  BC_Abund_New = paste0("BodyMass_dist_discrete_day",i-1)
  BodyCondition_Table[,(BC_Abund_New):= 
                        BC_Abund_Next/sum(BC_Abund_Next)]
  
  BM_PopDist = 
    as.vector(as.matrix(
      BodyCondition_Table[,paste0("BodyMass_dist_discrete_day",i-1),
                          with=F]))
  
  rm(DailyGain,ReducedForage)
  
  # Determine WSI- and body condition-dependent probability of departure
  WSIProbDepart = 
    as.matrix(NODE_DATA[,paste0("WSI_day_",i),with=F])

  ProbDepart_WSI = 
    NAFunc((weight_WSI * (((WSIProb(WSIProbDepart) + wsi_cutoff) ^ 3) / 
                            (((WSIProb(WSIProbDepart) + wsi_cutoff) ^ 3) + 
                            (wsi_cutoff ^ 3))) * 
             (WSIProbDepart > wsi_cutoff)))

  ProbDepart_WSI_mtx = matrix(rep(ProbDepart_WSI,length(bc_bins)),
                              ncol=length(bc_bins))
  
  ProbDepart = ProbDepart_WSI_mtx + ProbDepart_BC_mtx + 
    ProbDepart_DB_mtx
  
  AbundDepart_vec = rowSums(
    matrix((DayAbund %o% BM_PopDist),
           ncol=20) * ProbDepart)

  rm(ProbDepart,ProbDepart_WSI,ProbDepart_WSI_mtx)
  
  # Use Cobb-Douglas function to calculate node-specific
  # attractiveness.
  NodeAttract = paste0("NodeAttract_day_",i)
  
  NODE_DATA[,(NodeAttract):= as.numeric(
    (WSIProbDepart<=wsi_cutoff) * 
      (NODE_DATA[,(NormReducedForage),with=F]^
         WeightScheme[i,ForageAvailComponent]) * 
      (NormRoostShoreline^WeightScheme[i,RoostingHabComponent]) *
      (NormDistBreed^WeightScheme[i,BreedingDistComponent]) * 
      (GammaMoveProb^WeightScheme[i,GammaMoveComponent]))]
  
  NODE_DATA[,(NodeAttract) := lapply(.SD,NAFunc),.SDcols=NodeAttract]
  NODE_DATA[,(NodeAttract) := lapply(.SD,norm_func),
            .SDcols=NodeAttract]
  
  # Determine next day's abundance based on probability of staying
  # and departing.
  NextAbund = paste0("N_Abund_",i)
  NodeAttract = 
    as.matrix(NODE_DATA[,paste0("NodeAttract_day_",days[i]),with=F])
  
  CurrAbund = 
    (DayAbund - AbundDepart_vec) + 
    (rowSums(NodeAttract %o% AbundDepart_vec,na.rm=TRUE))
  
  # Redistribute departing population among body condition classes based 
  # on movement among nodes
  BC_DepartAbund = as.vector(
    colSums(AbundDepart_vec %o% BM_PopDist))
  BC_DepartAbund_mtx = matrix(rep(BC_DepartAbund,nrow(BC_dist_seq_mtx)),
                              ncol=length(bc_bins),byrow=TRUE)
  
  BC_DepartAbundPropor_dist = 
    data.table(cbind((BC_DepartAbund_mtx * BC_dist_seq_mtx),
                     cummul_flightcost[-1,BC_decrementvalue]))
  
  target_cols=paste0("V",1:20)
  BC_DepartAbundDist = BC_DepartAbundPropor_dist[,lapply(.SD,sum),by=V21,
                                                 .SDcols=target_cols]
  
  # Calculate number of inidivuals leaving each body condition
  BC_Abund_Loss = BC_DepartAbund - 
    colSums(BC_DepartAbundDist[2:nrow(BC_DepartAbundDist),c(2:21),
                               with=F])
  
  # Calculate number of individuals entering each body condition
  BC_DepartAbundDist_mtx = 
    as.matrix(BC_DepartAbundDist[,2:ncol(BC_DepartAbundDist)])
  element_indicator = row(BC_DepartAbundDist_mtx) - 
    col(BC_DepartAbundDist_mtx)
  
  BC_DepartAbundDist_diags = 
    split(BC_DepartAbundDist_mtx,element_indicator)
  
  BC_Abund_Gain = 
    unlist(lapply(BC_DepartAbundDist_diags,sum))[1:length(bc_bins)]
  
  BC_Abund_New = paste0("BodyMass_dist_discrete_day",i)
  BodyCondition_Table[,(BC_Abund_New):=
                        (BC_Abund_Loss + BC_Abund_Gain) / 
                        sum(BC_Abund_Loss + BC_Abund_Gain)]
  
  BM_PopDist = 
    as.vector(as.matrix(
      BodyCondition_Table[,paste0("BodyMass_dist_discrete_day",i),
                          with=F]))

  # Incur mortality
  ProporDie = rowSums(floor(
    matrix((CurrAbund %o% BM_PopDist),
           ncol=20)) * ProbMort_mtx)
  
  # Calulcate remaining abundance after mortality
  NODE_DATA[,(NextAbund) := CurrAbund - ProporDie]
  NODE_DATA[,sum(CurrAbund-ProporDie)]
  
  deadbirds[i] = sum(CurrAbund) - sum(NODE_DATA[,NextAbund,with=F])
  mortality = deadbirds[i]/N_0
  print(deadbirds[i])

  rm(DayAbund,NextAbund,NodeAttract,AbundDepart_vec,CurrAbund,
     ProporDie,NormReducedForage)
}
Sys.time()
end_timer=Sys.time()
total_time = end_timer - start_timer
total_time

plot(1:335,ZO_norm(deadbirds),type="l",ylab="Normalized value (0-1)",
     xlab="Day of migration")

# test=names(NODE_DATA)[which(names(NODE_DATA) %like% "WSI_day_")]
# example=NODE_DATA[,lapply(.SD,FUN=mean),.SDcols=test]
# lines(1:335,ZO_norm(example),col="red")
# 
# test=names(NODE_DATA)[which(names(NODE_DATA) %like% "NormForage_day_")]
# example=NODE_DATA[,lapply(.SD,FUN=mean),.SDcols=test]
# lines(1:335,ZO_norm(example),col="blue")

test=names(NODE_DATA)[which(names(NODE_DATA) %like% "DailyGain_day")]
example=NODE_DATA[,lapply(.SD,FUN=mean),.SDcols=test]
lines(1:335,ZO_norm(example),col="green")

# test=names(NODE_DATA)[which(names(NODE_DATA) %like% "NodeAttract_day")]
# example=NODE_DATA[,lapply(.SD,FUN=mean),.SDcols=test]
# lines(1:335,ZO_norm(example),col="gold")

text(c(250,250),c(1,0.95),c("Mortality","Daily Gain"),
     col=c("black","green"))

deadbirds = NODE_DATA[,sum(N_Abund_0)-sum(N_Abund_335)]

mortality = deadbirds/N_0
deadbirds
mortality

# deadbirds
# 10111486.42496241
# mortality
# 0.5092276733449994
#
# Diagnostics on metrics so far -------------------------------------------
# Node attractiveness
ggplot(data=NODE_DATA,aes(x=X_INDEX,y=Y_INDEX,col=ZO_norm(N_Abund_335))) +
  geom_point(size=2,shape=15) +
  scale_color_gradientn(colors=tim.colors(10)) +
  theme_minimal() +
  theme(
    axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  )
#
# 
# fwrite(NODE_DATA,file="NodeSpecifData_DAILYABUND.txt")
# 
# # Daily abundance
# NODE_DATA = fread("NodeSpecifData_DAILYABUND.txt")

# Animation for daily node abundance (variable scale)
DailyAbundAnimate_func = function() {
  for(i in 1:length(days)){
    print(
      ggplot(data=NODE_DATA,aes_string(x="X_INDEX",y="Y_INDEX",
                                       col=paste0("N_Abund_",i-1))) +
        geom_point(size=2,shape=15) +
        scale_color_gradientn(
          name=paste0("Abundance\non ",
                      format(as.Date(i,origin=NonbreedingStart),
                             format="%b %d")),
          colors=tim.colors(10)) + #,limits=c(0,70000)) +
        theme_minimal() +
        theme(
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
        )
    )
  }  
}

saveVideo(DailyAbundAnimate_func(),interval=0.1,
          video.name="DailyAbundance_VarScale.mp4")

# Animation for daily node abundance (0-1 scale)
DailyAbundAnimate_func = function() {
  for(i in 1:length(days)){
    DayAbund = 
      ZO_norm(as.matrix(NODE_DATA[,paste0("N_Abund_",i-1),with=F]))
    
    print(
      ggplot(data=NODE_DATA,aes_string(x="X_INDEX",y="Y_INDEX",
                                       col=DayAbund)) +
        geom_point(size=2,shape=15) +
        scale_color_gradientn(
          name=paste0("Abundance\non ",
                      format(as.Date(i,origin=NonbreedingStart),
                             format="%b %d")),
          colors=tim.colors(10)) + #,limits=c(0,70000)) +
        theme_minimal() +
        theme(
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
        )
    )
  }  
}

saveVideo(DailyAbundAnimate_func(),interval=0.1,
          video.name="DailyAbundance_NormScale.mp4")

# Animation for daily body condition class proportional abundance
DailyBCPAAnimate_func = function() {
  for(i in 1:length(days)){
    DayProportion = 
      ZO_norm(as.matrix(
        BodyCondition_Table[,paste0("BodyMass_dist_discrete_day",i),
                            with=F]))
    print(
      ggplot(data=BodyCondition_Table,
             aes(x=BC_bins,y=DayProportion)) +
        geom_line(lwd=1.2) + 
        xlab("Body Condition Class") + 
        ylab("Normalized Proportional Abundance") +
        ggtitle(paste0(format(as.Date(i,origin=NonbreedingStart),
                              format="%b %d"))) +
        theme_bw() +
        theme(
          axis.text.x=element_text(size=16,face="bold",colour="black"),
          axis.text.y=element_text(size=16,face="bold",colour="black"),
          axis.title=element_text(size=24,face="bold",colour="black"),
          plot.title=element_text(size=24,face="bold",colour="black")
          
        )
    )
  }  
}

saveVideo(DailyBCPAAnimate_func(),interval=0.1,
          video.name="DailyBCPA.mp4")


#







