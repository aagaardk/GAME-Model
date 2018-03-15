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
# By: Kevin Aagaard
# 
# Version edits - converting for-loop to data.table format
#
###########################################################################

# Initialize settings (workspace, directory, etc.) ------------------------

#Clear the workspace
rm(list=ls())

#Clear the console
cat("\014")

# #If Rtools isn't working right:
# Sys.setenv(PATH="%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin")
# #And test:
# Sys.getenv('PATH')

# For use when library path directory is not working 
#(permissioning issue)
# .libPaths(.libPaths()[2])

#Include necessary libraries
x = c("compiler","fields","maptools","data.table","plyr","foreign",
      "dplyr","ggplot2","scales","stringr","rstudioapi","animation",
      "profvis","directlabels","SDMTools","raster","sp","rgdal","rgeos",
      "viridis","geosphere","mapproj","mapdata","parallel","snow",
      'ggtern','doSNOW','doParallel','bit64')

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

# #Time all of the code
# Sys.time()
# start_timer = Sys.time()


#
# Define global landscape/movement parameters -----------------------------
# Define global geographical projection from LCC to WGS84
geo_prj =
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# How far birds can move within nodes, in meters
node_move = 3000 
# How fast birds can fly, km/hr
flight_vel = 75
# Cost of flight, kg per hour; 39,700 kJ per kg
flight_cost_hr = (226.4 / 39700) * 1000
# Cost of flight, kg per km; 39,700 kJ per kg
flight_cost_km = (3.1 / 39700) * 1000
# Number of days in non-breeding period
fall_start = as.Date("2014-09-01", tz="UTC")
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
year_from_record = 1957:2014
NonbreedingStart = paste0(year_from_record,"-09-01")

# Function to define temperature-dependent metabolic rate
TempDepMetRate =
  function(temperature,basalmetabolicrate) {
    (BMR_LCT(temperature,basalmetabolicrate) ^ 
       (temperature<LowCritTemp)) *
      (basalmetabolicrate^(temperature>=LowCritTemp))
  }

# BMR below LCT, from Prince 1979
BMR_LCT = function(temperature,basalmetabolicrate) {
  ((2.06 - (0.03 * temperature)) * kckJConv) * basalmetabolicrate
}

# Energy loss functions
# BMR in kJ, Heitmeyer 2010
BasalMetRate_func = function(bodymass) 422 * (bodymass^0.74)

# Energy gain functions (I THINK THIS VALUE REALLY MATTERS)
# 9 kcal/g, Ricklefs 1974 via Whyte and Bolen 1988
# MeanDailyMassIncrease = 4 #g, Hanson et al. 1990
# RWB_coeff = MeanDailyMassIncrease * (9 * kckJConv)
#
# LK_coeff = 4.6              #Lindstrom and Kvist 1995; passerines
# KL_coeff = 10               #Kvist and Lindstrom 2003; waders, ~3-10
# 
Empir_coeff = 13#1.8           #range of MetCostLandscape = 1:~2.2XBMR
# 
# Fuel deposition rate
# FDR = RWB_coeff + MetCostLandscape
# FDR = LK_coeff * BasalMetRate
FDR_func = function(basalmetabolicrate) Empir_coeff * basalmetabolicrate
# FDR = Empir_coeff * BasalMetRate

# Function to convert NA-elements to 0
NAFunc = function(x) {x[is.na(x)]<-0; x}

# Function to normalize landscape values
norm_func = function(x) {x / sum(x)}

# Zero-to-one normalization
ZO_norm = function(x) (x-min(x))/(max(x)-min(x))

# Z-score standardization
ZScore_func = function(x) (x - mean(x)) / sd(x)

# Converstion from kilocalories to kiloJoules
kckJConv = 4.184 # 1 kilocalorie = 4.184 kiloJoules

# Energetic parameters
LowCritTemp = 20 #Celsius

flight_dist= 2253.083 #km
range_limit=0.8

# Proportion of landscape that is inhospitable (WSI > wsi_cutoff)
inhospitable = function(x) length(which(x>wsi_cutoff))/num_nodes

# Proportion of landscape that is hospitable (WSI < wsi_cutoff)
hospitable = function(x) length(which(x<wsi_cutoff))/num_nodes

# Mean WSI within proportion of landscape that is hospitable
# (WSI < wsi_cutoff)
availmeanWSI = function(x) mean(x[x<wsi_cutoff],na.rm=T)

# Median WSI within proportion of landscape that is hospitable
# (WSI < wsi_cutoff)
availmedianWSI = function(x) median(x[x<wsi_cutoff],na.rm=T)

# Mean WSI across landscape 
ALLmeanWSI = function(x) mean(x,na.rm=T)

# Median WSI across landscape
ALLmedianWSI = function(x) median(x,na.rm=T)

# Top 2% most populace nodes per day
node_populace = function(x) which(x>quantile(x,probs=0.98))

# Center of mass of population per day
Pop_CoM = function(wt,x,y,z) COGravity(x,y,z,wt)

# Mask operation for arrays
fun.mask = function(x) {
  x[!arr_dummy_mask] = NA
  return(x)
}

# Rounding function
round_func = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

# Temperature mean in available habitat function
availmeanTemp = function(x) mean(x[y<wsi_cutoff])

# Gamma movement functions for body condition transitions
gamma_x_func = function(x) ((flight_dist-x)*(x>0))
t1_func = function(x) x^(gamma_a-1)

# Gamma alpha-paramter for body condition transitions
gamma_a = 1


#
# Input data --------------------------------------------------------------
# Load node specific data
# #Without Temperature or WSI data
NODE_DATA = fread("NodeSpecifData.txt")
setkey(NODE_DATA,Code,ID,Y_INDEX,X_INDEX,Longitude,Latitude)

NODE_DATA[,`:=`(
  NormRoostShoreline = RoostShoreline/ sum(RoostShoreline),
  NormDistBreed = (max(DistBreed) - DistBreed) /
    sum(max(DistBreed) - DistBreed),
  HarvestDisturb = 0)]

N_0 = NODE_DATA[,sum(N_Abund_0)]

# Determine number of nodes
num_nodes = nrow(NODE_DATA)

# Set node size in km (20 miles = 32.1869 km)
node_size = 32.1869

# Average body mass, in grams
mean_body_mass = 1140
min_body_mass = 800
max_body_mass = 1600
beta_mean_body_mass = mean(c(mean_body_mass,max_body_mass))

a=10
b=1.5

body_mass_pop = N_0
#body_mass_pop = 10000
body_mass = rbeta(body_mass_pop,a,b) * beta_mean_body_mass
range(body_mass)
body_mass = body_mass[body_mass > min_body_mass]

# plot(density(body_mass))


#
# ONLY RUN COMMENTED LINES IF "NodeSpecifData.txt" IS LOST ----------------
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

# RUN THESE LINES SEPARATELY FROM THE ABOVE SECTION, AS THEY REQUIRE
# ALL AVAILABLE MEMORY
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

# Create node-specific temperature measures -------------------------------
# Read in temperature data per node
# TempData = fread("NorAm_historical_air_temperature_node_info.txt")
# TempData[,`:=`(Y_INDEX=max(Y_INDEX)+1-Y_INDEX,
#                X_INDEX=X_INDEX-(min(X_INDEX)-1))]
# TempDataFrame = data.frame(TempData)
# data_columns = colnames(TempDataFrame)
# #
# for(i in 1:(length(year_from_record)-1)){
#   #subset out a particular migration period (don't want to base this on
#   #calendar year, but rather on annual migratory cycle; i.e., fall one year
#   #to spring the next)
#   # 1972 was the most recent "average" (relative to 1951-1980) year,
#   # differing from the mean during that time span by 0.01 deg C.
#   target_year = year_from_record[i]
#   nonbreeding_year = c(as.character(target_year),
#                        as.character(target_year+1))
#   target_dates=c(paste0(nonbreeding_year[1],c("_09_","_10_",
#                                               "_11_","_12_")),
#                  paste0(nonbreeding_year[2],c("_01_","_02_","_03_","_04_",
#                                               "_05_")))
# 
#   nonbreeding_period = data_columns[grepl(target_dates[1],data_columns)]
#   for(j in 2:length(target_dates)){
#     nonbreeding_period =
#       c(nonbreeding_period,data_columns[grepl(target_dates[j],
#                                               data_columns)])
#   }
#   nonbreeding_period =
#     nonbreeding_period[!(nonbreeding_period %like% "_02_29")]
# 
#   TempDT_NonbreedingPeriod =
#     data.table(TempDataFrame[,c('X_INDEX','Y_INDEX',nonbreeding_period)])
#   TempDT_NonbreedingPeriod[,Code:=paste0(
#     TempDT_NonbreedingPeriod[,X_INDEX],"_",
#     TempDT_NonbreedingPeriod[,Y_INDEX])]
#   TempDT_NonbreedingPeriod[,`:=`(X_INDEX=NULL,Y_INDEX=NULL)]
#   setcolorder(TempDT_NonbreedingPeriod,c(num_days+1,1:num_days))
#   setnames(TempDT_NonbreedingPeriod,names(TempDT_NonbreedingPeriod),
#            c("Code",paste0("Temperature_day_",1:num_days)))
# 
#   # Calculate daily gain in terms of grams
#   NODE_DATA_TEMPER = NODE_DATA
#   setkey(NODE_DATA_TEMPER,Code)
#   setkey(TempDT_NonbreedingPeriod,Code)
#   NODE_DATA_TEMPER = NODE_DATA_TEMPER[TempDT_NonbreedingPeriod]
# 
#   fwrite(NODE_DATA_TEMPER,
#          file=paste0("NodeSpecifData_Temperature",year_from_record[i],".txt"))
#   
#   rm(NODE_DATA_TEMPER)
# 
#   print(i)
# }
# 
# rm(TempData,TempDataFrame,TempDT_NonbreedingPeriod)
# 
# 
#
# Create node-specific WSI layers -----------------------------------------
# WSI_DATA = fread("NorAm_historical_wsi_node_info.txt")
# WSI_DATA[,`:=`(Y_INDEX=max(Y_INDEX)+1-Y_INDEX,
#                X_INDEX=X_INDEX-(min(X_INDEX)-1))]
# remove_leapday = "_02_29"
# WSI_DATA[,grep(remove_leapday,names(WSI_DATA)):=NULL]
# WSI_DATA_frame = data.frame(WSI_DATA)
# data_columns = colnames(WSI_DATA_frame)
# 
# # Use this with mean data:
# {
#   # WSI_DATA = fread("NorAm_mean_wsi_node_info.txt")
#   # WSI_DATA[,"02_29":=NULL]
#   # WSI_DATA[,`:=`(Y_INDEX=max(Y_INDEX)+1-Y_INDEX,
#   #                X_INDEX=X_INDEX-(min(X_INDEX)-1))]
#   # WSI_DATA_frame = data.frame(WSI_DATA)
#   # data_columns = colnames(WSI_DATA_frame)
#   # target_dates=c("09_","10_","11_","12_","01_","02_",
#   #                "03_","04_","05_")
#   # # Use this with either dataset:
#   # # Subset out a particular migration period (don't want to base this on
#   # # calendar year, but rather on annual migratory cycle; i.e., fall one year
#   # # to spring the next)
#   # nonbreeding_period = data_columns[grepl(target_dates[1],data_columns)]
#   # for(i in 2:length(target_dates)){
#   #   nonbreeding_period =
#   #     c(nonbreeding_period,data_columns[grepl(target_dates[i],
#   #                                             data_columns)])
#   # }
#   #
#   # WSI_DATA.nonbreeding_period =
#   #   data.table(WSI_DATA_frame[,c('X_INDEX','Y_INDEX',nonbreeding_period)])
#   # WSI_DATA.nonbreeding_period[,Code:=paste0(
#   #   WSI_DATA.nonbreeding_period[,X_INDEX],"_",
#   #   WSI_DATA.nonbreeding_period[,Y_INDEX])]
#   # setcolorder(WSI_DATA.nonbreeding_period,c(num_days+3,1:(num_days+2)))
#   # setnames(WSI_DATA.nonbreeding_period,names(WSI_DATA.nonbreeding_period),
#   #          c("Code","X_INDEX","Y_INDEX",paste0("WSI_day_",1:num_days)))
#   #
#   # setkey(WSI_DATA.nonbreeding_period,Code,X_INDEX,Y_INDEX)
#   #
#   # # ggplot(data=WSI_DATA.nonbreeding_period,
#   # #        aes(x=X_INDEX,y=Y_INDEX,col=log(WSI_day_1))) +
#   # #   geom_point(size=1.75,shape=15) +
#   # #   scale_color_gradientn(colors=tim.colors(10)) +
#   # #   theme_minimal() +
#   # #   theme(
#   # #     axis.text=element_blank(),
#   # #     axis.title=element_blank(),
#   # #     panel.grid.major=element_blank(),
#   # #     panel.grid.minor=element_blank()
#   # #   )
#   #
#   # setkey(NODE_DATA,Code,X_INDEX,Y_INDEX)
#   # NODE_DATA = NODE_DATA[WSI_DATA.nonbreeding_period]
#   #
#   # rm(WSI_DATA,WSI_DATA.nonbreeding_period,WSI_DATA_frame,data_columns)
#   #
#   # fwrite(NODE_DATA,file="NodeSpecifData_MeanWSI.txt")
# }
# 
# NODE_DATA_LIST = NULL
# for(i in 1:(length(year_from_record)-1)){
#   # Use this with historical data:
#   target_year = year_from_record[i]
#   nonbreeding_year = c(as.character(target_year),
#                        as.character(target_year+1))
#   target_dates=c(paste0(nonbreeding_year[1],c("_09_","_10_",
#                                               "_11_","_12_")),
#                  paste0(nonbreeding_year[2],c("_01_","_02_","_03_","_04_",
#                                               "_05_")))
# 
#   nonbreeding_period = data_columns[grepl(target_dates[1],data_columns)]
#   for(j in 2:length(target_dates)){
#     nonbreeding_period =
#       c(nonbreeding_period,data_columns[grepl(target_dates[j],
#                                               data_columns)])
#   }
# 
#   WSI_DATA.nonbreeding_period =
#     data.table(WSI_DATA_frame[,c('X_INDEX','Y_INDEX',nonbreeding_period)])
#   WSI_DATA.nonbreeding_period[,Code:=paste0(
#     WSI_DATA.nonbreeding_period[,X_INDEX],"_",
#     WSI_DATA.nonbreeding_period[,Y_INDEX])]
#   setcolorder(WSI_DATA.nonbreeding_period,c(num_days+3,1:(num_days+2)))
#   setnames(WSI_DATA.nonbreeding_period,names(WSI_DATA.nonbreeding_period),
#            c("Code","X_INDEX","Y_INDEX",paste0("WSI_day_",1:num_days)))
# 
#   setkey(WSI_DATA.nonbreeding_period,Code,X_INDEX,Y_INDEX)
# 
#   NODE_DATA =
#     fread(paste0("NodeSpecifData_Temperature",year_from_record[i],".txt"))
#   setkey(NODE_DATA,Code,X_INDEX,Y_INDEX)
#   NODE_DATA_LIST = NODE_DATA[WSI_DATA.nonbreeding_period]
# 
#   fwrite(NODE_DATA_LIST,
#          file=paste0("NodeSpecifData_WSI_TEMPERATURE",year_from_record[i],".txt"))
# 
#   rm(NODE_DATA_LIST)
# 
#   print(i)
# }
# 
# 
#
# Clean up environment ----------------------------------------------------
# gc()
# rm(WSI_DATA.nonbreeding_period,WSI_DATA_frame,WSI_DATA,data_columns,OutCols)
# gc()
# 
# 
#
# Generate body mass/body condition classes -------------------------------
n_bins = 21
bc_bins = c(1:20)

body_mass_optim = density(body_mass)$x[which.max(density(body_mass)$y)]
# body_mass_discrete = c(seq(range(body_mass)[1],body_mass_optim,
#                            length.out=length(bc_bins)),
#                        range(body_mass)[2])
body_mass_discrete = seq(range(body_mass)[1],range(body_mass)[2],
                           length.out=n_bins)
body_mass_norm = ZO_norm(body_mass_discrete)
body_mass_lower = body_mass_norm[1:length(body_mass_norm)-1]
body_mass_upper = body_mass_norm[2:length(body_mass_norm)]

integrate_beta = 
  Vectorize(
    function(l,h) {
      integrate(function(x) dbeta(x,a,b),lower=l,upper=h)$value
    }
  )

body_mass_dist = sapply(body_mass_upper,integrate_beta,l=body_mass_lower)

BodyCondition_Table = 
  data.table(BC_bins=c(1,bc_bins+1),
             BodyMass_g=body_mass_discrete,
             BodyMass_dist_discrete_day0=
               c(0,body_mass_dist[1,] - 
                   data.table::shift(body_mass_dist[1,],
                                     n=1L,fill=0,type="lag")),
             BodyFat_g=
               body_mass_discrete*
               seq(0,0.13,length.out=n_bins),
             BodyMass_BMR=BasalMetRate_func(body_mass_discrete))

rm(body_mass)

#Parameter for gamma distribution to determine movement probability
BodyCondition_Table[,BodyFat_gamma_b:=(flight_dist * ZO_norm(BodyFat_g))]

#Bottom bin has 99.1% daily survivorship
min_survival = 0.9974 #0.991
#Top bin has 99.997% daily survivorship
max_survival = .99997

daily_survival = round(c(0,seq(min_survival,max_survival,
                               length.out=length(bc_bins))),5)

BodyCondition_Table[,Survivorship:=daily_survival]

#Possible flight distance by body condition (in km)
BodyCondition_Table[,MaxLike_flightdist:=
                      range_limit * (BodyFat_g/flight_cost_km)]

#Calculate Pr(flying distance, x) based on body condition
g=0
j=flight_dist

integrate_unif = function(l,h) {
  integrate(Vectorize(function(x) dunif(x,g,j)),
            lower=l,upper=h)$value
}

flightdist_distrib = sapply(BodyCondition_Table[,MaxLike_flightdist],
                            integrate_unif,l=0)

# Calculate the cost of traveling each unit (node) of a flight
cummul_flightcost = 
  data.table(dist_km = seq(0,flight_dist,node_size),
             cost_km=c(0,cumsum(rep((flight_cost_km*node_size),
                                    (flight_dist/node_size)))))

# Body condition-dependent flight distance probabilities
BC_distprob = BodyCondition_Table[,ifelse(MaxLike_flightdist==0,0,
                                          MaxLike_flightdist)]

# function to identify the possible movement distances for each body condition;
# since only birds with a probability of departing > 0 will have this applied
# to them, we remove consideration of the probability of moving 0 km
prob_dist_func = Vectorize(function(x) {
  cummul_flightcost[dist_km<x,dist_km[dist_km>0]]})

BC_ProbMoveDist = prob_dist_func(BC_distprob)

BC_gamma_x = lapply(BC_ProbMoveDist,gamma_x_func)
BC_gamma_b = flight_dist-BC_distprob


t1 = lapply(BC_gamma_x,t1_func)
t2 = vector('list')
for(i in 1:n_bins){
  t2[[i]] = BC_gamma_x[[i]]/BC_gamma_b[i]
}
t3 = gamma(gamma_a)

BC_gamma_moveprob = vector('list')
for(i in 1:n_bins){
  BC_gamma_moveprob[[i]] = ((t1[[i]])*(exp(1)^(-t2[[i]])/(BC_gamma_b[i]^gamma_a*t3)))
}

# Convert ragged BC_gamma_moveprob list to data.frame
rowMax = max(sapply(BC_gamma_moveprob, length)) 
BC_gamma_moveprob_dt = data.table(
  do.call(cbind, lapply(BC_gamma_moveprob, 
                        function(x){ 
                          length(x) <- rowMax 
                          x 
                        })))
setnames(BC_gamma_moveprob_dt,names(BC_gamma_moveprob_dt),
         c(paste0('BC_',1:n_bins)))

# Normalize the gamma values for relativized movement probabilities
BC_gamma_moveprob_dt = NAFunc(BC_gamma_moveprob_dt)
BC_gamma_moveprob_mtx = apply(BC_gamma_moveprob_dt,2,norm_func)
BC_gamma_moveprob_mtx = NAFunc(BC_gamma_moveprob_mtx)

# Add the possible movement distances for reference
BC_gamma_moveprob_dt = data.table(BC_gamma_moveprob_mtx)
BC_gamma_moveprob_dt[,PossibleDistances:=BC_ProbMoveDist[[n_bins]]]
setcolorder(BC_gamma_moveprob_dt,c(n_bins+1,1:n_bins))

# Calculate the number of bins by which birds will be reduced for flying
# different distances
BC_transition_g_ls = list()
for(i in n_bins:2){
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
BC_transition_g_mtx = rbind(BC_transition_g_mtx,
                            rep(0,n_bins-1))
BC_transition_g_mtx = cbind(rep(0,n_bins),
                            BC_transition_g_mtx)

BC_transitions = 
  sort(BC_transition_g_mtx[,n_bins][
    which(BC_transition_g_mtx[,n_bins]>0)])

BC_decrement = vector()
BC_decrement_total = rep(0,nrow(cummul_flightcost))
for(i in 1:length(BC_transitions)){
  BC_decrement=
    ifelse(cummul_flightcost[,cost_km]<BC_transitions[i],0,1)
  BC_decrement_total = BC_decrement + BC_decrement_total
}

cummul_flightcost[,BC_decrementvalue:=BC_decrement_total]   

# Add body condition decrement per km-bin traveled to movement probability
# matrix
BC_gamma_moveprob_decrements = 
  cummul_flightcost[cummul_flightcost[,dist_km] %in% 
                      BC_gamma_moveprob_dt[,PossibleDistances],
                    BC_decrementvalue]

BC_gamma_moveprob_dt[,BC_decrementvalue:=BC_gamma_moveprob_decrements]
BC_gamma_moveprob_dt[,PossibleDistances:=NULL]

sum_cols = paste0('BC_',1:n_bins)
BC_transitionprobs = BC_gamma_moveprob_dt[,lapply(.SD,sum),.SDcols=sum_cols,
                                          by=BC_decrementvalue]

# Calculate number of individuals entering each body condition using
# a transition matrix. Need to convert from data table with 'from' in 
# columns and 'to' as: 
#           column number - body condition decrement value of row
# to transition matrix with from in columns and to in rows

# fill in matrix based on row and column id
BC_Transition_mtx=matrix(NA,ncol=21,nrow=21)
for(i in 1:nrow(BC_Transition_mtx)){
  for(j in 1:ncol(BC_Transition_mtx)){
    BC_Transition_mtx[i,j] = 
      as.numeric(BC_transitionprobs[BC_decrementvalue==j-i,
                                    (j+1),with=F])
  }
}

BC_Transition_mtx = NAFunc(BC_Transition_mtx)


#
# Define parameters and rules for movement outside of loop ----------------
# Weights for WSI, body condition, and distance to breeding node for
# probability of departure functions
weight_WSI = 0.33
weight_BC = 0.33
weight_DB = 0.34

# Probability of departure component related to WSI
# The relation below comes from Schummer et al. 2010
WSIProb = function(x) -0.6282+(0.0547*x)+(0.0046*((x)^2))

# Probability of departure component related to body condition
ProbDepart_BC = c(0,(weight_BC * ((c(1:n_bins)^8) / 
                                    ((c(1:n_bins)^8) + 
                                       ((n_bins-3)^8))))[-1])
ProbDepart_BC_mtx = matrix(rep(ProbDepart_BC,nrow(NODE_DATA)),
                           ncol=n_bins,byrow=T)

# Probability of departure component related to distance to nearest 
# breeding node
breedingdistances = NODE_DATA[,DistBreed]
n_breedingdistances = max(breedingdistances)
ProbDepart_DB = weight_DB * (breedingdistances ^ 5) / 
  ((breedingdistances ^ 5) + 
     ((n_breedingdistances-(n_breedingdistances)*0.15)) ^ 5)
ProbDepart_DB_mtx = matrix(rep(ProbDepart_DB,n_bins),
                           ncol=n_bins,byrow=T)

ProbMort = 1 - BodyCondition_Table[,Survivorship]
ProbMort_mtx = matrix(rep(ProbMort,nrow(NODE_DATA)),ncol=n_bins,
                      byrow=T)

# Natural forage decay rates by landcover type
# Convert this to a matrix with dimensions dependent upon the number
# of landcover types used
Shoreline_decay = 0.9998
Crops_decay = 0.997
WoodyWetlands_decay = 0.9965
HerbWetlands_decay = 0.991

# Set values for CD exponents
W_Forage = 0.55
W_Roost = 0.25
W_Breed = 0.8
W_Gamma = 0.1

# Establish inflection point--date on which importance of landscape
# variables switch directions
set.seed(7)
inflection = floor(runif(1,(num_days/2)-15,(num_days/2)+15))

# Generate sequence of weights
inflection_seq = c(seq(1,0.1,length.out=inflection),
                   seq(0.1,1,length.out=num_days-inflection))

WeightScheme = data.table(
  DayofNonBreedingPeriod = 1:num_days,
  ForageAvailComponent = ((((-inflection_seq^3) * W_Forage) +
                             (max(inflection_seq)^3)*W_Forage) + 0.03),
  RoostingHabComponent = ((((-inflection_seq^3) * W_Roost) +
                             (max(inflection_seq)^3) * W_Roost) + 0.03),
  BreedingDistComponent = (((inflection_seq^3) * W_Breed) + 0.03),
  GammaMoveComponent = W_Gamma + 0.01
)

WeightScheme_Time =
  melt(WeightScheme,id.vars="DayofNonBreedingPeriod")

# WeightScheme_Plot =
# ggplot(data=WeightScheme_Time, aes(x=DayofNonBreedingPeriod,y=value,
#                                    linetype=variable)) +
# geom_line(lwd=1.25) +
# scale_linetype_manual(values=c("solid","longdash",
#                                "dotdash","dotted"),
#                       guide='none') +
# geom_text(aes(x=inflection,y=0),label="Breeding",size=5) +
# geom_text(aes(x=inflection,y=0.15),label="Gamma",size=5) +
# geom_text(aes(x=inflection,y=0.625),label="Forage",size=5) +
# geom_text(aes(x=inflection,y=0.325),label="Roosting",size=5) +
# ylim(0,1) +
# labs(x="Day of Nonbreeding Period",y="Relative weight") +
# theme_minimal() +
# theme(
#   axis.text=element_text(face="bold",size=16),
#   axis.title=element_text(face="bold",size=20)
# )
# ggsave("CobbDouglasWeights.jpeg",WeightScheme_Plot,width=5,height=5)


#
# Impose movement rules ---------------------------------------------------
rm(NODE_DATA)

deadbirds = vector()
COMPLETE_DATA = NULL

Sys.time()
start_timer = Sys.time()

n_cores = detectCores() - 1
start_cluster = makeCluster(n_cores)
registerDoParallel(start_cluster)

foreach(j = 1:(length(year_from_record)-1),
        .packages = c('data.table'),
        .verbose = T) %dopar% {
          
          COMPLETE_DATA =
            fread(paste0("NodeSpecifData_WSI_TEMPERATURE",
                         year_from_record[j],".txt"))
          
          for(i in 1:length(days)) {
            # Reduce forage material based on mass- and temperature-dependent
            # daily gain for each body mass/condition class and landcover-specific
            # decay rates, with place holders for decay rates
            DayAbund =
              as.matrix(COMPLETE_DATA[,paste0("N_Abund_",i-1),with=F])
            DayTemper =
              as.matrix(COMPLETE_DATA[,paste0("Temperature_day_",i),with=F])
            DailyGain = paste0("DailyGain_day_",i)
            ReducedForage = paste0("Forage_day_",i)
            
            BM_PopDist =
              as.vector(as.matrix(
                BodyCondition_Table[,paste0("BodyMass_dist_discrete_day",i-1),
                                    with=F]))
            
            # Daily gain is calculated with the inputs:
            #   HumanDisturb (Degree of human disturbance--urbanization)
            #   HarvestDisturb (Degree of human disturbance--hunter harvest)
            #   FDR (Fuel deposition rate, kJ ingested per day)
            #   TD_BMR (Temperature dependent basal metabolic rate, kJ per day)
            # DG = ((1 - HumanDisturb) * (1 - HarvestDisturb) * FDR) - TD_BMR
            
            # Calculate the number of birds in each body condition class
            # in each node on the day
            BC_Node_Abund =
              matrix((DayAbund %o% BM_PopDist),
                     ncol=n_bins)
            
            # Calculate the basal metabolic rate for each body
            # condition class in each node on the day
            BMR = BodyCondition_Table[,matrix(rep(BodyMass_BMR,each=num_nodes),
                                              ncol=n_bins)]
            
            # Calculate the disturbance-dependent FDR
            Disturb_mtx =
              COMPLETE_DATA[,matrix(
                rep(((1 - HumanDisturb) * (1 - HarvestDisturb)),n_bins),
                ncol=n_bins)]
            
            Disturb_FDR = Disturb_mtx * FDR_func(BMR)
            
            rm(Disturb_mtx)
            
            # Calculate the temperature-dependent basal metabolic
            # rate for each body condition class in each node on the day
            DayTemper_BC = matrix(rep(DayTemper,n_bins),ncol=n_bins)
            TD_BMR = TempDepMetRate(DayTemper_BC,BMR)
            rm(DayTemper_BC)
            
            # Calculate the daily gain in each node on the day
            Node_Disturb_FDR = rowSums(Disturb_FDR)
            COMPLETE_DATA[,(DailyGain):= Node_Disturb_FDR]
            
            # Calculate the forage material reduction
            COMPLETE_DATA[,(ReducedForage):=
                            (((ForageShoreline * (Shoreline_decay^i)) +
                                (ForageCrops * (Crops_decay^i)) +
                                (ForageWoodyWetlands *
                                   (WoodyWetlands_decay^i)) +
                                (ForageHerbWetlands *
                                   (HerbWetlands_decay^i)))) -
                            Node_Disturb_FDR]
            
            NormReducedForage = paste0("NormForage_day_",i)
            incols = ReducedForage
            COMPLETE_DATA[,(NormReducedForage) := lapply(.SD,norm_func),
                          .SDcols=incols]
            
            # Calculate kilojoules gained per body condition to redistribute
            # individuals among body condition classes after foraging
            DG_BC_Node = (Disturb_FDR - TD_BMR) / 39700
            
            DG_BC_Transition = BodyCondition_Table[,outer(BodyFat_g,BodyFat_g,
                                                          '-')]
            
            DG_BC_Node_Transition = vector("list")
            for(k in 1:n_bins){
              DG_BC_Node_Transition[[k]] = findInterval(DG_BC_Node[,k],
                                                        DG_BC_Transition[,k])
            }
            
            DG_BC_Node_Transition_mtx =
              matrix(unlist(DG_BC_Node_Transition), ncol=n_bins,
                     byrow=F)
            
            BC_Abund = matrix((DayAbund %o% BM_PopDist),ncol=n_bins)
            
            BC_Abund_Gain_func = function(x) {
              sum(BC_Abund[which(DG_BC_Node_Transition_mtx==x)])
            }
            
            BC_Abund_Next = sapply(1:n_bins,BC_Abund_Gain_func)
            BC_Abund_New = paste0("BodyMass_dist_discrete_day",i)
            BodyCondition_Table[,(BC_Abund_New):=
                                  BC_Abund_Next/sum(BC_Abund_Next)]
            
            rm(DailyGain,ReducedForage,DG_BC_Node,BC_Abund_Next,
               BC_Abund_New)
            
            # Determine WSI- and body condition-dependent probability of departure
            WSIProbDepart =
              as.matrix(COMPLETE_DATA[,paste0("WSI_day_",i),with=F])
            
            ProbDepart_WSI =
              NAFunc((weight_WSI * (((WSIProb(WSIProbDepart) + wsi_cutoff) ^ 3) /
                                      (((WSIProb(WSIProbDepart) + wsi_cutoff) ^ 3) +
                                         (wsi_cutoff ^ 3))) *
                        (WSIProbDepart > wsi_cutoff)))
            
            ProbDepart_WSI_mtx = matrix(rep(ProbDepart_WSI,n_bins),
                                        ncol=n_bins)
            
            ProbDepart = ProbDepart_WSI_mtx + ProbDepart_BC_mtx +
              ProbDepart_DB_mtx
            
            AbundDepart_vec = rowSums(
              matrix((DayAbund %o% BM_PopDist),
                     ncol=n_bins) * ProbDepart)
            
            rm(ProbDepart,ProbDepart_WSI,ProbDepart_WSI_mtx)
            
            # Use Cobb-Douglas function to calculate node-specific
            # attractiveness.
            NodeAttract = paste0("NodeAttract_day_",i)
            
            COMPLETE_DATA[,(NodeAttract):= as.numeric(
              (WSIProbDepart<=wsi_cutoff) *
                (COMPLETE_DATA[,(NormReducedForage),with=F]^
                   WeightScheme[i,ForageAvailComponent]) *
                (NormRoostShoreline^WeightScheme[i,RoostingHabComponent]) *
                (NormDistBreed^WeightScheme[i,BreedingDistComponent]) *
                (GammaMoveProb^WeightScheme[i,GammaMoveComponent]))]
            
            COMPLETE_DATA[,(NodeAttract) := lapply(.SD,NAFunc),
                          .SDcols=NodeAttract]
            COMPLETE_DATA[,(NodeAttract) := lapply(.SD,norm_func),
                          .SDcols=NodeAttract]
            
            rm(WSIProbDepart)
            
            # Determine next day's abundance based on probability of staying
            # and departing.
            NextAbund = paste0("N_Abund_",i)
            NodeAttract =
              as.matrix(COMPLETE_DATA[,paste0("NodeAttract_day_",days[i]),
                                      with=F])
            
            gc()
            CurrAbund =
              (DayAbund - AbundDepart_vec) +
              (rowSums(NodeAttract %o% AbundDepart_vec,na.rm=TRUE))
            gc()
            
            # Redistribute departing population among body condition classes based
            # on movement among nodes
            BC_AbundDepart_vec = as.vector(colSums(AbundDepart_vec %o% BM_PopDist))
            BC_AbundStay_vec = colSums(BC_Abund) - BC_AbundDepart_vec
            
            BC_AbundDepart_mtx =
              matrix(rep(BC_AbundDepart_vec,each=n_bins),ncol=n_bins)
            
            BC_AbundArrive_mtx = BC_Transition_mtx * BC_AbundDepart_mtx
            BC_AbundArrive_vec = rowSums(BC_AbundArrive_mtx)
            
            BC_Abund_Move = norm_func(BC_AbundStay_vec + BC_AbundArrive_vec)
            
            BM_PopDist =
              as.vector(as.matrix(
                BodyCondition_Table[,paste0("BodyMass_dist_discrete_day",i),
                                    with=F]))
            
            BC_Abund_Next = (BM_PopDist + BC_Abund_Move) / 2
            
            BC_Abund_New = paste0("BodyMass_dist_discrete_day",i)
            BodyCondition_Table[,(BC_Abund_New):=BC_Abund_Next]
            
            BM_PopDist =
              as.vector(as.matrix(
                BodyCondition_Table[,paste0("BodyMass_dist_discrete_day",i),
                                    with=F]))
            
            # Incur mortality
            ProporDie = rowSums(floor(
              matrix((CurrAbund %o% BM_PopDist),
                     ncol=n_bins)) * ProbMort_mtx)
            
            # Calulcate remaining abundance after mortality
            COMPLETE_DATA[,(NextAbund) := CurrAbund - ProporDie]
            COMPLETE_DATA[,sum(CurrAbund-ProporDie)]
            
            deadbirds[i] = sum(CurrAbund) - sum(COMPLETE_DATA[,NextAbund,
                                                              with=F])
            mortality = deadbirds[i]/N_0
            print(mortality*100)
            
            # Remove objects from environment to optimize memory use
            rm(DayAbund,NextAbund,NodeAttract,AbundDepart_vec,CurrAbund,
               ProporDie,NormReducedForage,BC_Abund)
            
            print(i)
          }
  
          # Save COMPLETE_DATA and BodyCondition_Table elements and remove
          # to reduce memory usage
          fwrite(COMPLETE_DATA,
                 file=paste0(year_from_record[j],"_Prediction.txt"))
          fwrite(BodyCondition_Table,
                 file=paste0(year_from_record[j],"_BCTable.txt"))
          rm(COMPLETE_DATA)
          
          print(j)
        }

stopCluster(start_cluster)

rm(ProbDepart_BC_mtx,ProbDepart_DB_mtx,ProbMort_mtx,DG_BC_Node_Transition,
   DG_BC_Node_Transition_mtx,breedingdistances,ProbDepart_DB,
   integrate_beta,integrate_unif,BC_DepartAbund_mtx,WeightScheme,
   BC_DepartAbundPropor_dist,BC_dist_seq_ls,BC_dist_seq_ls_maxln,
   BC_dist_seq_mtx,BC_distprob,BC_Abund_Gain,BC_DepartAbundDist,
   BC_DepartAbundDist_diags,BC_DepartAbundDist_mtx)


Sys.time()
end_timer=Sys.time()
total_time = end_timer - start_timer
total_time


#
# Read in predictive migration dynamics for analysis ----------------------
# availhabitat = MeanWSI_availhabitat = MedianWSI_availhabitat = 
#   MeanWSI_allhabitat = MedianWSI_allhabitat = vector('list')
# mortality = autumn_mortality = winter_mortality = spring_mortality =
#   annmean_availhabitat = annsd_availhabitat =annmin_availhabitat =
#   median_forage = vector()
# annualavailhabitat = c(1:num_days)
# for(i in 1:(length(year_from_record)-1)){
#   COMPLETE_DATA = fread(paste0(year_from_record[i],"_Prediction.txt"))
# 
#   # Mortality over non-breeding season
#   last_day_abund = paste0("N_Abund_",length(days))
#   mortality[i] = COMPLETE_DATA[,sum(N_Abund_0) - sum(.SD),
#                                .SDcols=last_day_abund]
# 
#   # Mortality over autumn migration
#   autumn_mortality[i] =
#     COMPLETE_DATA[,sum(N_Abund_0) - sum(N_Abund_90)]
# 
#   # Mortality over winter
#   winter_mortality[i] =
#     COMPLETE_DATA[,sum(N_Abund_91) - sum(N_Abund_180)]
# 
#   # Mortality over spring migration
#   spring_mortality[i] = COMPLETE_DATA[,sum(N_Abund_181) - sum(.SD),
#                                .SDcols=last_day_abund]
# 
#   # Availability of habitat across non-breeding period (based on WSI)
#   target_cols =
#     names(COMPLETE_DATA)[which(names(COMPLETE_DATA) %like% "WSI_day_")]
#   availhabitat[[i]] =
#     COMPLETE_DATA[,lapply(.SD,hospitable),.SDcols=target_cols]
#   annmean_availhabitat[i] = mean(as.matrix(availhabitat[[i]]))
#   annsd_availhabitat[i] = sd(as.matrix(availhabitat[[i]]))
#   annmin_availhabitat[i] = min(as.matrix(availhabitat[[i]]))
# 
#   # Gather proportion of available habitat through time
#   annualavailhabitat =
#     cbind(annualavailhabitat,
#           as.vector(as.matrix(unlist(availhabitat[[i]]))))
# 
#   # Average WSI of available habitat
#   MeanWSI_availhabitat[[i]] =
#     as.matrix(COMPLETE_DATA[,lapply(.SD,availmeanWSI),
#                             .SDcols=target_cols])
#   
#   # Median WSI of available habitat
#   MedianWSI_availhabitat[[i]] =
#     as.matrix(COMPLETE_DATA[,lapply(.SD,availmedianWSI),
#                             .SDcols=target_cols])
#   # Average WSI of all habitat
#   MeanWSI_allhabitat[[i]] =
#     as.matrix(COMPLETE_DATA[,lapply(.SD,ALLmeanWSI),
#                             .SDcols=target_cols])
#   
#   # Median WSI of all habitat
#   MedianWSI_allhabitat[[i]] =
#     as.matrix(COMPLETE_DATA[,lapply(.SD,ALLmedianWSI),
#                             .SDcols=target_cols])
#   
#   # Clean up environment
#   rm(target_cols,COMPLETE_DATA)
# 
#   print(i)
# }
# 
# norm_mortality = ZO_norm(mortality)
# norm_autumn_mortality = ZO_norm(autumn_mortality)
# norm_winter_mortality = ZO_norm(winter_mortality)
# norm_spring_mortality = ZO_norm(spring_mortality)
# norm_annmean_availhabitat = ZO_norm(annmean_availhabitat)
# norm_annmin_availhabitat = ZO_norm(annmin_availhabitat)
# norm_AnnMeanWSI_availhabitat =
#   ZO_norm(unlist(lapply(MeanWSI_availhabitat,mean)))
# norm_AnnMedianWSI_availhabitat =
#   ZO_norm(unlist(lapply(MedianWSI_availhabitat,mean)))
# norm_AnnMeanWSI_allhabitat =
#   ZO_norm(unlist(lapply(MeanWSI_allhabitat,mean)))
# norm_AnnMedianWSI_allhabitat =
#   ZO_norm(unlist(lapply(MedianWSI_allhabitat,mean)))
# 
# SummaryStats =
#   data.table(Years=year_from_record[1:(length(year_from_record)-1)],
#              Mortality=mortality,
#              NormMortality=norm_mortality,
#              AutumnMortality=autumn_mortality,
#              NormAutumnMortality=norm_autumn_mortality,
#              WinterMortality=winter_mortality,
#              NormMortality=norm_winter_mortality,
#              SpringMortality=spring_mortality,
#              NormSpringMortality=norm_spring_mortality,
#              Mean_AvailHab=annmean_availhabitat,
#              StDev_AvailHab=annsd_availhabitat,
#              NormMean_AvailHab=norm_annmean_availhabitat,
#              Min_AvailHab=annmin_availhabitat,
#              NormAnnMinHab=norm_annmin_availhabitat,
#              MeanWSI_AvailHab=unlist(lapply(MeanWSI_availhabitat,mean)),
#              MedianWSI_AvailHab=unlist(lapply(MedianWSI_availhabitat,mean)),
#              MeanWSI_AllHab=unlist(lapply(MeanWSI_allhabitat,mean)),
#              MedianWSI_AllHab=unlist(lapply(MedianWSI_allhabitat,mean)),
#              NormAnnMeanWSI_AvailHab=norm_AnnMeanWSI_availhabitat,
#              NormAnnMedianWSI_AvailHab=norm_AnnMedianWSI_availhabitat,
#              NormAnnMeanWSI_AvailHab=norm_AnnMeanWSI_allhabitat,
#              NormAnnMedianWSI_AllHab=norm_AnnMedianWSI_allhabitat)
# 
# fwrite(SummaryStats,file="FullRun_SummaryStats.txt")

SummaryStats = fread("FullRun_SummaryStats.txt")

# Mean/Median mortality
SummaryStats[,mean((N_0-Mortality)/N_0)]
SummaryStats[,median((N_0-Mortality)/N_0)]
SummaryStats[,range((N_0-Mortality)/N_0)]
# SummaryStats[,mean((N_0-Mortality)/N_0) + 
#                ((1.96*sd((N_0-Mortality)/N_0))) / 
#                nrow(SummaryStats)]
# SummaryStats[,mean((N_0-Mortality)/N_0) - 
#                ((1.96*sd((N_0-Mortality)/N_0))) / 
#                nrow(SummaryStats)]

SummaryStats[,mean((N_0-AutumnMortality)/N_0)]
SummaryStats[,median((N_0-AutumnMortality)/N_0)]
SummaryStats[,range((N_0-AutumnMortality)/N_0)]

SummaryStats[,mean(((N_0-AutumnMortality)-WinterMortality) / 
                     (N_0-AutumnMortality))]
SummaryStats[,median(((N_0-AutumnMortality)-WinterMortality) / 
                       (N_0-AutumnMortality))]
SummaryStats[,range(((N_0-AutumnMortality)-WinterMortality) / 
                      (N_0-AutumnMortality))]

SummaryStats[,mean((((N_0-AutumnMortality)-WinterMortality) - 
                      SpringMortality) / 
                     ((N_0-AutumnMortality)-WinterMortality))]
SummaryStats[,median((((N_0-AutumnMortality)-WinterMortality) - 
                        SpringMortality) / 
                       ((N_0-AutumnMortality)-WinterMortality))]
SummaryStats[,range((((N_0-AutumnMortality)-WinterMortality) - 
                       SpringMortality) / 
                      ((N_0-AutumnMortality)-WinterMortality))]
# 
# Mort_Plot =
#   ggplot(data=SummaryStats,aes(x=Years,y=Mortality)) +
#   geom_point(size=2) + geom_line() +
#   xlab("Years") +
#   ylab("Mortality") +
#   theme_minimal() +
#   theme(
#     axis.text=element_text(face="bold",size=16),
#     axis.title=element_text(face="bold",size=20)
#   )
# ggsave(paste0("Mortality.jpeg"),
#        Mort_Plot,width=5,height=5)
# dev.off()
# 
# NormMort_Plot =
#   ggplot(data=SummaryStats,aes(x=Years,y=NormMortality)) +
#   geom_point() +
#   geom_smooth(method='glm', color='black') +
#   xlab("Years") +
#   ylab("Mortality") +
#   theme_minimal()
# ggsave(paste0("NormalizedMortality.jpeg"),
#        NormMort_Plot,width=5,height=5)
# dev.off()
# 
# NormMeanAvailHab_Plot =
#   ggplot(data=SummaryStats,aes(x=Years,y=NormMean_AvailHab)) +
#   geom_point(size=2) +
#   geom_smooth(method='glm', color='black') +
#   xlab("Years") +
#   ylab("Proporitonal Available Habitat") +
#   theme_minimal()
# ggsave(paste0("NormalizedAvailHab.jpeg"),
#        NormMeanAvailHab_Plot,width=5,height=5)
# dev.off()
# 
# MinAvailHab_Plot =
#   ggplot(data=SummaryStats,aes(x=Years,y=Min_AvailHab)) +
#   geom_point(size=2) + 
#   geom_smooth(method='glm', color='black') +
#   ylim(0,SummaryStats[,max(Min_AvailHab)]) + xlab("Years") +
#   ylab("Minimum Available Habitat") +
#   theme_minimal()
# ggsave(paste0("ProporMinAnnAvailHab.jpeg"),
#        MinAvailHab_Plot,width=5,height=5)
# dev.off()
# 
# NormAnnMinHab_Plot =
#   ggplot(data=SummaryStats,aes(x=Years,y=NormAnnMinHab)) +
#   geom_point(size=2) + geom_line() +
#   xlab("Years") +
#   ylab("Normalized Annual Minimum Available Habitat") +
#   theme_minimal()
# ggsave(paste0("NormalizedAnnMinHab.jpeg"),
#        NormAnnMinHab_Plot,width=5,height=5)
# dev.off()
# 
# Combined_Plot =
#   ggplot(data=SummaryStats,aes(x=Years,y=NormMortality)) +
#   geom_smooth(method='glm', color='darkred', fill='darkred') +
#   geom_smooth(aes(x=Years,y=NormAnnMinHab),
#               method='glm', color='darkblue', fill='darkblue') +
#   ylim(c(0,1)) +
#   xlab("Years") +
#   ylab("Normalized Metrics") +
#   theme_minimal()
# ggsave(paste0("AllVariables.jpeg"),
#        Combined_Plot,width=8,height=5)
# dev.off()
# 
# Mort_MinHab_Plot =
#   ggplot(data=SummaryStats,aes(x=NormAnnMinHab,y=NormMortality)) +
#   geom_point(size=2) +
#   xlab("Normalized Minimum Annual Available Habitat") +
#   ylab("Normalized Mortality") +
#   theme_minimal()
# ggsave(paste0("Mortality_v_ProporMinAnnAvailHab.jpeg"),
#        Mort_MinHab_Plot,width=5,height=5)
# dev.off()
# 
# Norm_Mort_MeanHab_Plot =
#   ggplot(data=SummaryStats,aes(x=NormMean_AvailHab,y=NormMortality)) +
#   geom_point(size=2) +
#   geom_smooth(method='glm', color='black') +
#   geom_vline(xintercept=0.62, col='gray50',lty='dashed') +
#   ylim(c(0,1)) +
#   xlab("Available Habitat") +
#   ylab("Mortality") +
#   theme_minimal()
# ggsave(paste0("Norm_Mortality_v_ProporMeanAnnAvailHab_GLM.jpeg"),
#        Norm_Mort_MeanHab_Plot,width=5,height=5)
# dev.off()
# 
# Norm_Mort_MeanHab_LoessPlot =
#   ggplot(data=SummaryStats,aes(x=NormMean_AvailHab,y=NormMortality)) +
#   geom_point(size=2) +
#   geom_smooth(method='loess', color='black') +
#   geom_vline(xintercept=0.62, col='gray50',lty='dashed') +
#   ylim(c(0,1)) +
#   xlab("Available Habitat") +
#   ylab("Mortality") +
#   theme_minimal()
# ggsave(paste0("Norm_Mortality_v_ProporMeanAnnAvailHab_LOESS.jpeg"),
#        Norm_Mort_MeanHab_LoessPlot,width=5,height=5)
# dev.off()
# 
# Mort_MeanHab_Plot =
#   ggplot(data=SummaryStats,aes(x=Mean_AvailHab,y=Mortality)) +
#   geom_point(size=2) +
#   xlab("Normalized Mean Annual Available Habitat") +
#   ylab("Normalized Mortality") +
#   theme_minimal()
# ggsave(paste0("Mortality_v_ProporMeanAnnAvailHab.jpeg"),
#        Mort_MeanHab_Plot,width=5,height=5)
# dev.off()
# 
MeanWSI_Plot =
  ggplot(data=SummaryStats,aes(x=Years,y=MeanWSI_AvailHab)) +
  geom_point(size=2) +
  geom_smooth(method='glm', color='black') +
  ylim(c(SummaryStats[,floor(min(MeanWSI_AvailHab))],
         SummaryStats[,ceiling(max(MeanWSI_AvailHab))])) +
  xlab("Years") +
  ylab("Mean WSI of Available Habitat") +
  theme_minimal()
ggsave(paste0("MeanWSIAvailHab.jpeg"),
       MeanWSI_Plot,width=5,height=5)
dev.off()
# 
MedianWSI_Plot =
  ggplot(data=SummaryStats,aes(x=Years,y=MedianWSI_AvailHab)) +
  geom_point(size=2) +
  geom_smooth(method='glm', color='black') +
  ylim(c(-8.5,-6.8)) +
  xlab("Years") +
  ylab("Median WSI of Available Habitat") +
  theme_minimal()
ggsave(paste0("MedianWSIAvailHab.jpeg"),
       MedianWSI_Plot,width=5,height=5)
dev.off()
# 
MeanWSIALL_Plot =
  ggplot(data=SummaryStats,aes(x=Years,y=MeanWSI_AllHab)) +
  geom_point(size=2) +
  geom_smooth(method='glm', color='black') +
  ylim(c(SummaryStats[,floor(min(MeanWSI_AllHab))],
         SummaryStats[,ceiling(max(MeanWSI_AllHab))])) +
  xlab("Years") +
  ylab("Mean WSI of Landscape") +
  theme_minimal()
ggsave(paste0("MeanWSIAllHab.jpeg"),
       MeanWSIALL_Plot,width=5,height=5)
dev.off()
# 
MedianWSIALL_Plot =
  ggplot(data=SummaryStats,aes(x=Years,y=MedianWSI_AllHab)) +
  geom_point(size=2) +
  geom_smooth(method='glm', color='black') +
  ylim(c(SummaryStats[,round(min(MedianWSI_AllHab),1)],
         SummaryStats[,round(max(MedianWSI_AllHab),1)])) +
  xlab("Years") +
  ylab("Median WSI of Landscape") +
  theme_minimal()
ggsave(paste0("MedianWSIAllHab.jpeg"),
       MedianWSIALL_Plot,width=5,height=5)
dev.off()

NormMeanWSI_Plot =
  ggplot(data=SummaryStats,aes(x=Years,y=NormMeanWSI_AvailHab)) +
  geom_point(size=2) +
  geom_smooth(method='glm', color='black') +
  xlab("Years") +
  ylab("Normalized Mean WSI of Available Habitat") +
  theme_minimal()
ggsave(paste0("NormMeanWSIAvailHab.jpeg"),
       NormMeanWSI_Plot,width=5,height=5)
dev.off()
# 
Mort_HabWSI_Plot =
  ggplot(data=SummaryStats,aes(x=NormMeanWSI_AvailHab,y=NormMortality)) +
  geom_point(size=2) +
  xlab("Normalized Mean WSI of Available Habitat") +
  ylab("Normalized Mortality") +
  theme_minimal()
ggsave(paste0("Mortality_v_MeanWSIAvailHab.jpeg"),
       Mort_HabWSI_Plot,width=5,height=5)
dev.off()
# 
# MeanHab_HabWSI_Plot =
#   ggplot(data=SummaryStats,
#          aes(x=NormMeanWSI_AvailHab,y=NormMean_AvailHab)) +
#   geom_point(size=2) +
#   xlab("Normalized Mean WSI of Available Habitat") +
#   ylab("Normalized Mean Annual Available Habitat") +
#   theme_minimal()
# ggsave(paste0("ProporMeanAnnAvailHab_v_MeanWSIAvailHab.jpeg"),
#        MeanHab_HabWSI_Plot,width=5,height=5)
# dev.off()
# 
# annualavailhabitat = data.table(annualavailhabitat)
# setnames(annualavailhabitat,names(annualavailhabitat),
#          c('MigDay',
#            paste0('AvailHab_',
#                   year_from_record[1:(length(year_from_record)-1)])))
# annualavailhabitat = melt(annualavailhabitat,id.vars='MigDay',
#                           variable.name='Year',value='ProporAvailHab')
# 
# Annual_AvailHab_Plot = 
#   ggplot(data=annualavailhabitat,
#          aes(x=MigDay,y=ProporAvailHab,group=Year)) +
#   geom_line(alpha=0.1) +
#   scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1)) +
#   xlab("Day of Non-breeding Period") +
#   ylab("Normalized Mean Annual Available Habitat") +
#   theme_minimal()
# ggsave(paste0("ProporAnnAvailHab_OverTime.jpeg"),
#        Annual_AvailHab_Plot,width=5,height=5)
# dev.off()
#
# # # With state lines
# # na_countries_sts = c("USA","CAN")
# # northamerica_sts =
# #   do.call("bind",lapply(na_countries_sts, function(x)
# #     getData("GADM",country=x,level=1)))
# # 
# # # Remove Alaska and Hawaii
# # Can_CONUS = c("Alaska","Hawaii")
# # CONUS_NA_sts =
# #   northamerica_sts[!(toupper(northamerica_sts$NAME_1) %like%
# #                        "ALASKA"|
# #                        toupper(northamerica_sts$NAME_1) %like%
# #                        "HAWAII"),]
# # 
# # writeOGR(obj=CONUS_NA_sts,dsn=getwd(),layer="states",
# #          driver="ESRI Shapefile")
# CONUS_NA_sts = readOGR("states",dsn=getwd())
# # plot(CONUS_NA_sts)
# #
# # Obtain spatial data for plotting
# CoM_DailyDist_Between = CoM_DailyDist_Among = CoM_DailyDir =
#   CenterofMass_coords = vector('list')
# SumCoM_AdjDist = MeanCoM_AdjDist_Between = MeanCoM_AdjDist_Among =
#   CoM_NS_Dist = vector()
# # setwd(paste0(getwd(),'/Output'))
# for(i in 1:(length(year_from_record)-1)){
#   COMPLETE_DATA = fread(paste0(year_from_record[i],"_Prediction.txt"))
# 
#   target_cols =
#     names(COMPLETE_DATA)[
#       which(names(COMPLETE_DATA) %like% "N_Abund_")][2:(num_days+1)]
# 
#   most_populace_node =
#     as.vector(as.matrix(
#       COMPLETE_DATA[,lapply(.SD,node_populace),.SDcols=target_cols]))
#   COMPLETE_DATA[most_populace_node,within_populace:=1]
#   COMPLETE_DATA[,within_populace:=NAFunc(within_populace)]
# 
#   Mig_route = COMPLETE_DATA[,.(Longitude,Latitude,within_populace)]
#   coordinates(Mig_route) = c("Longitude","Latitude")
#   gridded(Mig_route) = TRUE
#   Mig_route_grd = as(Mig_route, "SpatialGridDataFrame")
#   crs(Mig_route_grd) =
#     "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
#   +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#   Mig_route_rstr = raster(Mig_route_grd)
#   Mig_route_proj = projectRaster(Mig_route_rstr,crs=geo_prj,
#                                   res=0.1)
# 
#   CenterofMass = t(as.matrix(
#     COMPLETE_DATA[,lapply(.SD,Pop_CoM,x=Longitude,y=Latitude,z=NULL),
#                   .SDcols=target_cols][c(1,3)]))
#   CenterofMass_pts = SpatialPoints(CenterofMass)
# 
#   crs(CenterofMass_pts) =
#     "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
#       +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#   CenterofMass_proj = spTransform(CenterofMass_pts,CRS(geo_prj))
#   CenterofMass_coords[[i]] = data.table(coordinates(CenterofMass_proj))
#   CenterofMass_coords[[i]][,ColorIterate:=
#                              c(1:nrow(CenterofMass_coords[[i]]))]
# 
#   CoM_DailyDist_Between[[i]] =
#     sapply(2:nrow(CenterofMass_coords[[i]]),
#            function(x) {
#              distm(CenterofMass_coords[[i]][x-1,.(coords.x1,coords.x2)],
#                    CenterofMass_coords[[i]][x,.(coords.x1,coords.x2)]) /
#                1000
#            }
#     )
# 
#   CoM_DailyDist_Among[[i]] =
#     distm(CenterofMass_coords[[i]][,.(coords.x1,coords.x2)]) / 1000
#   CoM_DailyDist_Among[[i]] =
#     ifelse(CoM_DailyDist_Among[[i]]==0,NA,CoM_DailyDist_Among[[i]])
# 
#   CoM_DailyDir[[i]] =
#     CenterofMass_coords[[i]][,coords.x2 -
#                                data.table::shift(coords.x2,
#                                                  n=1L,fill=0,type="lag")][-1]
# 
# 
#   CoM_NS_Dist[i] =
#     distm(CenterofMass_coords[[i]][1,.(coords.x1,coords.x2)],
#           CenterofMass_coords[[i]][which.min(coords.x2),
#                                    .(coords.x1,coords.x2)]) / 1000
# 
#   SumCoM_AdjDist[i] = sum(CoM_DailyDist_Between[[i]])
# 
#   MeanCoM_AdjDist_Between[i] = mean(CoM_DailyDist_Between[[i]])
#   MeanCoM_AdjDist_Among[i] = mean(CoM_DailyDist_Among[[i]],na.rm=T)
# 
#   col_pal = colorRampPalette(c('darkred','deepskyblue4'))
#   CenterofMass_coords[[i]][,ColorGradient:=
#                              col_pal(10)[as.numeric(
#                                cut(ColorIterate,breaks = 10))]]
# 
#   # tiff(file=paste0(year_from_record[i],
#   #                  '_MigrationRoute_CoM_STPRVBRDR.tiff'),
#   #      bg="transparent",width=1500,height=900)
#   # par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#   # image(Mig_route_proj,col=c('white','gray75'),
#   #       bty='n',xlab=NA,ylab=NA,xaxt='n',yaxt='n')
#   # plot(CONUS_NA_sts,lty="dotted",border="gray",add=TRUE)
#   # # plot(CONUS_NA,lty="dotted",add=TRUE)
#   # points(CenterofMass_coords[[i]][,.(coords.x1,coords.x2)],
#   #        pch=16,cex=0.5,col=CenterofMass_coords[[i]][,ColorGradient])
#   # dev.off()
# 
#   rm(COMPLETE_DATA)
# }
# 
# DistanceBetweenMetric = data.table(do.call(cbind,CoM_DailyDist_Between))
# setnames(DistanceBetweenMetric,names(DistanceBetweenMetric),
#          c(paste0("DistanceBetween_",
#                   year_from_record[-which.max(year_from_record)])))
# DistanceAmongMetric = 
#   matrix(rep(MeanCoM_AdjDist_Among,each=num_days-1),
#          nrow=num_days-1)
# DistanceAmongMetric = data.table(DistanceAmongMetric)
# setnames(DistanceAmongMetric,names(DistanceAmongMetric),
#          c(paste0("DistanceAmong_",
#                   year_from_record[-which.max(year_from_record)])))
# DirectionMetric = data.table(do.call(cbind,CoM_DailyDir))
# setnames(DirectionMetric,names(DirectionMetric),
#          c(paste0("Direction_",
#                   year_from_record[-which.max(year_from_record)])))
# 
# MovementMetrics = cbind(DayofPeriod=days[-which.max(days)],
#                         DistanceBetweenMetric,DistanceAmongMetric,
#                         DirectionMetric)
# fwrite(MovementMetrics,file="MovementMetrics.txt")

MovementMetrics = fread("MovementMetrics.txt")

# for(i in 1:(length(year_from_record)-1)){
#   CenterofMass_coords[[i]][,`:=`(
#     Year=rep(year_from_record[i],nrow(CenterofMass_coords[[i]])),
#     MeanWSIAvail=as.vector(MeanWSI_availhabitat[[i]]),
#     MedianWSIAvail=as.vector(MedianWSI_availhabitat[[i]]),
#     MeanWSIAll=as.vector(MeanWSI_allhabitat[[i]]),
#     MedianWSIAll=as.vector(MedianWSI_allhabitat[[i]]),
#     Mean_AvailHab=as.vector(as.matrix(availhabitat[[i]])),
#     Norm_AvailHab=as.vector(ZO_norm(unlist(availhabitat[[i]]))))]
# }
# 
# CenterofMass_coords = data.table(do.call(rbind,CenterofMass_coords))
# fwrite(CenterofMass_coords,file="CenterofMass_coords.txt")

CenterofMass_coords = fread("CenterofMass_coords.txt")
CenterofMass_coords[,`:=`(
  meanDailyLat=mean(coords.x2),
  meanDailyWSIAvail=mean(MeanWSIAvail),
  medianDailyWSIAvail=mean(MedianWSIAvail),
  meanDailyWSIAll=mean(MeanWSIAll),
  medianDailyWSIAll=mean(MedianWSIAll),
  meanDailyAvailHab=mean(Mean_AvailHab)),
  by=ColorIterate]
CenterofMass_coords[,`:=`(ColorIterate=NULL,ColorGradient=NULL)]

CenterofMass_coords[,`:=`(
  AnnualMeanWSIAvail=mean(MeanWSIAvail),
  AnnualMedianWSIAvail=mean(MedianWSIAvail),
  AnnualMeanWSIAll=mean(MeanWSIAll),
  AnnualMedianWSIAll=mean(MedianWSIAll)),
  by=Year]

col_pal = colorRampPalette(c('red','blue'))
CenterofMass_coords[,`:=`(
  ColorGradient=col_pal(57)[as.numeric(cut(MeanWSIAvail,breaks = 10))],
  NonbreedingDay=rep(1:num_days,(length(year_from_record)-1)))]

scalemin = CenterofMass_coords[,round(min(AnnualMedianWSIAvail),1)]
scalemax = CenterofMass_coords[,round(max(AnnualMedianWSIAvail),1)]

LatAbundTime_plot =
  ggplot(data=CenterofMass_coords,
         aes(x=NonbreedingDay,y=coords.x2,group=Year,
             col=AnnualMedianWSIAvail)) +
  geom_line(stat='smooth',method='loess',se=FALSE,size=1.2,alpha=0.4) +
  scale_color_continuous(name='Annual Median \nWeather Severity',
                         breaks=c(seq(scalemin,scalemax,
                                      (scalemax-scalemin)/6)),
                         low='gray90',high='black') +
  xlab("Day of Nonbreeding Period") +
  ylab("Latitude") +
  theme_minimal()
ggsave("LatAbundTime.jpeg",LatAbundTime_plot,width=5,height=5)
dev.off()

LatAbundWSI = 
  CenterofMass_coords[,.(min(coords.x2),
                         MedianWSIAvail[which.min(coords.x2)]),
                      by=Year]
setnames(LatAbundWSI,names(LatAbundWSI),c('Year','Latitude','WSI'))

LatAbundWSI_Plot = 
  ggplot(data=LatAbundWSI,aes(x=WSI,y=Latitude)) +
  geom_point() +
  geom_smooth(method='glm',color='black') +
  xlab("Weather Severity") +
  ylab("Latitude") +
  theme_minimal()
ggsave("LatAbundWSI.jpeg",LatAbundWSI_Plot,width=5,height=5)
dev.off()

#Regression of Most populace minimum latitude with WSI
MPML_WSI = 
  CenterofMass_coords[,.(min(coords.x2),min(MeanWSIAvail),
                         mean(MeanWSIAvail),min(MedianWSIAvail),
                         mean(MedianWSIAvail),min(MeanWSIAll),
                         mean(MeanWSIAll),min(MedianWSIAll),
                         mean(MedianWSIAll),min(Mean_AvailHab),
                         mean(Mean_AvailHab)),by=Year]
setnames(MPML_WSI,names(MPML_WSI),
         c('Year','Latitude','Min_MeanWSIAvail','Mean_MeanWSIAvail',
           'Min_MedianWSIAvail','Mean_MedianWSIAvail','Min_MeanWSIAll',
           'Mean_MeanWSIAll','Min_MedianWSIAll','Mean_MedianWSIAll',
           'Min_Mean_AvailHab','AnnualMean_DailyMeanAvailableHabitat'))


# Evaluate and visualize results ------------------------------------------
MPML_Year_plot =
  ggplot(data=MPML_WSI,
         aes(x=Year,y=Latitude)) +
  geom_point() +
  geom_smooth(method='glm') +
  xlab("Year") +
  ylab("Latitude") +
  theme_minimal()
ggsave("MPML_Year.jpeg",MPML_Year_plot,width=5,height=5)
dev.off()

MPML_WSI_plot =
  ggplot(data=MPML_WSI,
         aes(x=Mean_MedianWSIAll,y=Latitude)) +
  geom_point() +
  geom_smooth(method='glm') +
  # scale_x_continuous(breaks=c(-9,-8.5,-8,-7.5,-7,-6.5,-6)) +
  xlab("WSI") +
  ylab("Latitude") +
  theme_minimal()
ggsave("MPML_WSI.jpeg",MPML_WSI_plot,width=5,height=5)
dev.off()

Year_WSI_plot =
  ggplot(data=MPML_WSI,
         aes(x=Year,y=Mean_WSI)) +
  geom_point() +
  geom_smooth(method='glm') +
  scale_y_continuous(breaks=c(-9,-8.5,-8,-7.5,-7,-6.5,-6)) +
  xlab("Year") +
  ylab("WSI") +
  theme_minimal()
ggsave("Year_WSI.jpeg",Year_WSI_plot,width=5,height=5)
dev.off()

# SummaryStats[,`:=`(TotaDist_Traveled=SumCoM_AdjDist,
#                    MeanDist_Between=MeanCoM_AdjDist_Between,
#                    MeanDist_Among=MeanCoM_AdjDist_Among,
#                    NorthSouthDist_Traveled=CoM_NS_Dist)]
# 
# fwrite(SummaryStats,file="FullRun_SummaryStats.txt")

SummaryStats = fread("FullRun_SummaryStats.txt")

Dist_Among =
  ggplot(data=SummaryStats,aes(x=Years,y=MeanDist_Among)) +
  geom_point(size=2) + geom_line(lty="dotted") +
  xlab("Years") +
  ylab("Mean Distance (km) Among all Center-of-mass Points") +
  theme_minimal()
ggsave(paste0("Among_Dist.jpeg"),
       Dist_Among,width=5,height=5)
dev.off()

Dist_Btwn =
  ggplot(data=SummaryStats,aes(x=Years,y=MeanDist_Between)) +
  geom_point(size=2) + geom_line(lty="dotted") +
  xlab("Years") +
  ylab("Mean Distance (km) Between Adjacent Center-of-mass Points") +
  theme_minimal()
ggsave(paste0("Between_Dist.jpeg"),
       Dist_Btwn,width=5,height=5)
dev.off()

Dist_NS =
  ggplot(data=SummaryStats,aes(x=Years,y=NorthSouthDist_Traveled)) +
  geom_point(size=2) + geom_line(lty="dotted") +
  xlab("Years") +
  ylab("Mean Distance (km) Between Northern and Southern Extremes") +
  theme_minimal()
ggsave(paste0("NS_Dist.jpeg"),
       Dist_NS,width=5,height=5)
dev.off()

# Abundance-weighted population center-of-mass across years
CenterofMass_data = data.table()
for(i in 1:(length(year_from_record)-1)){
  COMPLETE_DATA = fread(paste0(year_from_record[i],"_Prediction.txt"))

  CenterofMass = t(as.matrix(
    COMPLETE_DATA[,lapply(.SD,Pop_CoM,x=Longitude,y=Latitude,z=NULL),
                  .SDcols=target_cols][c(1,3)]))
  CenterofMass_pts = SpatialPoints(CenterofMass)

  crs(CenterofMass_pts) =
    "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
  +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  CenterofMass_proj = spTransform(CenterofMass_pts,CRS(geo_prj))
  CenterofMass_data_init = data.table(
    Year=rep(year_from_record[i],nrow(coordinates(CenterofMass_proj))),
    coordinates(CenterofMass_proj))
  setnames(CenterofMass_data_init,names(CenterofMass_data_init),
           c("Year","Longitude","Latitude"))
  CenterofMass_data = rbind(CenterofMass_data, CenterofMass_data_init)

  rm(COMPLETE_DATA)
  print(i)
}

fwrite(CenterofMass_data,file="CenterofMass_data.txt")

CenterofMass_data = fread("CenterofMass_data.txt")
CenterofMass_data[,NonBreedingDay:=rep(days,(length(year_from_record)-1))]
CenterofMass_data[,`:=`(mean_Longitude=mean(Longitude),
                        mean_Latitude=mean(Latitude)),
                  by=NonBreedingDay]

# Obtain most populace nodes
for(i in 1:(length(year_from_record)-1)){
  COMPLETE_DATA = fread(paste0(year_from_record[i],"_Prediction.txt"))

  # Migration route (based on most populous node per day)
  target_cols =
    names(COMPLETE_DATA)[
      which(names(COMPLETE_DATA) %like% "N_Abund_")][2:(num_days+1)]
  most_populace_node =
    as.vector(as.matrix(
      COMPLETE_DATA[,lapply(.SD,node_populace),.SDcols=target_cols]))
  COMPLETE_DATA[most_populace_node,within_populace:=1]

  Mig_route = COMPLETE_DATA[,.(Longitude,Latitude,within_populace)]
  coordinates(Mig_route) = c("Longitude","Latitude")
  gridded(Mig_route) = TRUE
  Mig_route_grd = as(Mig_route, "SpatialGridDataFrame")
  crs(Mig_route_grd) =
    "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
    +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  Mig_route_rstr = raster(Mig_route_grd)
  Mig_route_proj = projectRaster(Mig_route_rstr,crs=geo_prj,
                                 res=0.1)

  image(Mig_route_proj,col=plot_colors,bty='n',xlab=NA,ylab=NA,
        xaxt='n',yaxt='n',add=T)

  rm(COMPLETE_DATA)
  print(i)
}

plot(CONUS_NA_sts,lty="dotted",border="gray",add=TRUE)

for(i in 1:(length(year_from_record)-1)){
  lines(CenterofMass_data[Year %like% year_from_record[i],Longitude],
        CenterofMass_data[Year %like% year_from_record[i],Latitude],
        lwd=1.5,col='gray40')
}

lines(CenterofMass_data[,unique(mean_Longitude)],
      CenterofMass_data[,unique(mean_Latitude)],
      lwd=3,col='black')

dev.off()
dev.off()
  

#
# Animations to visualize output over time --------------------------------
# Create masking work-flow to improve efficiency
COMPLETE_DATA = fread(paste0(year_from_record[1],"_Prediction.txt"))

DayAbund = 
  ZO_norm(as.matrix(COMPLETE_DATA[,N_Abund_0]))

DayNodeAbund = cbind(COMPLETE_DATA[,.(Longitude,Latitude)],DayAbund)
coordinates(DayNodeAbund) = c("Longitude","Latitude")
gridded(DayNodeAbund) = TRUE
DayNodeAbund_grd = as(DayNodeAbund, "SpatialGridDataFrame")
crs(DayNodeAbund_grd) =
  "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
        +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
DayNodeAbund_rstr = raster(DayNodeAbund_grd)
DayNodeAbund_proj = projectRaster(DayNodeAbund_rstr,crs=geo_prj,
                                  res=0.1)

dummy_rstr = raster(nrows=nrow(DayNodeAbund_proj),
                    ncols=ncol(DayNodeAbund_proj),
                    xmn=xmin(DayNodeAbund_proj),
                    xmx=xmax(DayNodeAbund_proj),
                    ymn=ymin(DayNodeAbund_proj),
                    ymx=ymax(DayNodeAbund_proj))
dummy_rstr[] = 1

dummy_mask = mask(dummy_rstr,CONUS_NA_sts)
arr_dummy_mask = as.array(dummy_mask)[,,1]

arr_dummy_mask[which(is.na(arr_dummy_mask))] = 0
arr_dummy_mask[which(arr_dummy_mask != 0)] = 1

storage.mode(arr_dummy_mask) = "logical"

# Animation for daily node abundance (0-1 scale)
for(j in 1:(length(year_from_record)-1)) {
  COMPLETE_DATA = fread(paste0(year_from_record[j],"_Prediction.txt"))
  
  DailyAbundAnimate_func = function() {
    for(i in 1:length(days)){
      Sys.time()
      start_timer = Sys.time()
      DayAbund = 
        ZO_norm(as.matrix(COMPLETE_DATA[,paste0("N_Abund_",i-1),with=F]))
      
      DayNodeAbund = cbind(COMPLETE_DATA[,.(Longitude,Latitude)],DayAbund)
      coordinates(DayNodeAbund) = c("Longitude","Latitude")
      gridded(DayNodeAbund) = TRUE
      DayNodeAbund_grd = as(DayNodeAbund, "SpatialGridDataFrame")
      crs(DayNodeAbund_grd) =
        "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
        +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
      DayNodeAbund_rstr = raster(DayNodeAbund_grd)
      DayNodeAbund_proj = projectRaster(DayNodeAbund_rstr,crs=geo_prj,
                                        res=0.1)
      crop_bounds = extent(CONUS_NA_sts)
      DayNodeAbund_crp = crop(DayNodeAbund_proj,crop_bounds)
      DayNodeAbund_array = as.array(DayNodeAbund_proj)[,,1]
      DayNodeAbund_array_masked = fun.mask(DayNodeAbund_array)
      DayNodeAbund_msk = raster(DayNodeAbund_array_masked,
                                xmn=xmin(DayNodeAbund_proj),
                                xmx=xmax(DayNodeAbund_proj),
                                ymn=ymin(DayNodeAbund_proj),
                                ymx=ymax(DayNodeAbund_proj))
      
      target_cols =
        names(COMPLETE_DATA)[
          which(names(COMPLETE_DATA) %like% "N_Abund_")][2:(num_days+1)]
      
      CenterofMass = t(as.matrix(
        COMPLETE_DATA[,lapply(.SD,Pop_CoM,x=Longitude,y=Latitude,z=NULL),
                      .SDcols=target_cols][c(1,3)]))
      CenterofMass_pts = SpatialPoints(CenterofMass)
      
      crs(CenterofMass_pts) =
        "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
    +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
      CenterofMass_proj = spTransform(CenterofMass_pts,CRS(geo_prj))
      CenterofMass_coords = data.table(coordinates(CenterofMass_proj))
      CenterofMass_coords[,ColorIterate:=c(1:nrow(CenterofMass_coords))]
      
      col_pal = colorRampPalette(c('white','gray50'))
      CenterofMass_coords[,ColorGradient:=
                            col_pal(10)[as.numeric(
                              cut(ColorIterate,breaks = 10))]]

      par(mar=c(0.1,0.1,0.1,0.1),oma=c(0.1,0.1,0.1,0.1),bty='n')
      plot(DayNodeAbund_msk,col=viridis(20),zlim=c(0,1),
           xlim=c(xmin(CONUS_NA_sts),xmax(CONUS_NA_sts)),
           ylim=c(ymin(CONUS_NA_sts),ymax(DayNodeAbund_msk)-5),
           axes=FALSE,box=FALSE,useRaster=TRUE,
           legend.width=1,legend.shrink=0.2,
           legend.args=
             list(text=paste0('Normalized Abundance\non ',
                              format(as.Date(i-1,
                                             origin=NonbreedingStart[j]),
                                     format="%d %b %Y")),
                  side=3,font=2,line=1,cex=0.8))
      plot(CONUS_NA_sts,lty="dotted",border="gray",add=TRUE)
      points(CenterofMass_coords[c(1:i),.(coords.x1,coords.x2)],
             pch=16,cex=0.5,col=CenterofMass_coords[,ColorGradient])
      print(i)
      Sys.time()
      end_timer = Sys.time()
      timetaken=end_timer-start_timer
      print(timetaken)
    }  
  }
  
  saveVideo(DailyAbundAnimate_func(),interval=0.1,
            video.name=paste0("DailyAbundance_NormScale",
                              year_from_record[j],".mp4"))
  print(j)
}

# Animation for daily body condition class proportional abundance
for(j in 1:(length(year_from_record)-1)) {
  BC_Table = fread(paste0(year_from_record[j],"_BCTable.txt"))
  
  DailyBCPAAnimate_func = function() {
    for(i in 1:length(days)){
      DayProportion =
        ZO_norm(as.matrix(
          BC_Table[,paste0("BodyMass_dist_discrete_day",i),
                   with=F]))
      print(
        ggplot(data=BC_Table,
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
      print(i)
    }
  }
  saveVideo(DailyBCPAAnimate_func(),interval=0.1,
            video.name=paste0("DailyBCPA",year_from_record[j],
                              ".mp4"))
  print(j)
}


#







