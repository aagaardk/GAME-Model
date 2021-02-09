###########################################################################
#
# Generalizable Avian Movement and Energetics Model
#--------------------------------------------------------------------------
#
# Created by Kevin Aagaard, Wayne Thogmartin, and Eric Lonsdorf
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
# (permissioning issue)
# .libPaths(.libPaths()[2])

# # Install archived version of SDMTools as it is not maintained 
# # (as of 06/25/2020)
# require(devtools)
# install_version(
#   'SDMTools', 
#   version = '1.1-221.2', 
#   repos = 'http://cran.us.r-project.org'
# )

#Include necessary libraries
x = c('compiler','fields','maptools','data.table','plyr','foreign',
      'dplyr','scales','stringr','rstudioapi','animation','ggplot2',
      'profvis','directlabels','SDMTools','raster','sp','rgdal','rgeos',
      'viridis','geosphere','mapproj','mapdata','parallel','snow',
      'doSNOW','doParallel','bit64','piecewiseSEM','fGarch','arules')

# If you don't know if you've installed the packages for some of these 
# libraries, run this:
# install.packages(x)

lapply(x, library, character.only=TRUE)
rm(x)

#Set working directory (should be universal)
setwd(
  dirname(
    rstudioapi::callFun(
      'getActiveDocumentContext'
    )$path
  )
)

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
# Input data --------------------------------------------------------------
# Load node specific data
# Without Temperature or WSI data
NODE_DATA = 
  fread(
    "NodeSpecifData.txt"
  )

setkey(
  NODE_DATA,
  Code,
  ID,
  Y_INDEX,
  X_INDEX,
  Longitude,
  Latitude
)

# Convert to shapefile
sp_node_landscape = NODE_DATA

sp_node_landscape[
  ,
  ForageShoreline :=
    as.numeric(
      ForageShoreline
    )
  ]

coordinates(
  sp_node_landscape
) = 
  sp_node_landscape[
    ,
    .(
      Longitude,
      Latitude
    )
    ]

crs(
  sp_node_landscape
) = 
  CRS(
    '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 
    +y_0=0 +datum=NAD83 +units=m +no_defs'
  )

writeOGR(
  obj = sp_node_landscape,
  dsn = '.',
  layer = 'sp_node_landscape',
  driver = 'ESRI Shapefile',
  overwrite_layer = TRUE
)

# Continue with data input
NODE_DATA[
  ,
  `:=`(
    NormRoostShoreline = 
      RoostShoreline / 
      sum(
        RoostShoreline
      ),
    NormDistBreed = 
      (
        max(
          DistBreed
        ) - 
          DistBreed
      ) /
      sum(
        max(
          DistBreed
        ) - 
          DistBreed
      ),
    HarvestDisturb = 0
  )
  ]

N_0 = 
  NODE_DATA[
    ,
    sum(
      N_Abund_0
    )
    ]

# Determine number of nodes
num_nodes = 
  nrow(
    NODE_DATA
  )

# Set node size in km (20 miles = 32.1869 km)
node_size = 32.1869

# Converstion from kilocalories to kiloJoules
kckJConv = 4.184 # 1 kilocalorie = 4.184 kiloJoules

flight_dist = 2253.083 #km
range_limit = 0.8

# Gamma alpha-paramter for body condition transitions
gamma_a = 1


#
# Generic-use functions ---------------------------------------------------
# Function to convert NA-elements to 0
NAFunc = 
  function(
    x
  ) {
    x[
      is.na(
        x
      )
      ] = 0
    x
  }

# Function to normalize landscape values
norm_func = 
  function(
    x, 
    na.rm = TRUE
  ) {
    x / 
      sum(
        x, 
        na.rm = na.rm
      )
  }

# Zero-to-one standardization
ZO_std = 
  function(
    x, 
    na.rm = TRUE
  ) {
    (
      x - 
        min(
          x, 
          na.rm = na.rm
        )
    ) / 
      (
        max(
          x, 
          na.rm = na.rm
        ) - 
          min(
            x, 
            na.rm = na.rm
          )
      )
  }

# Z-score standardization
ZScore_func = 
  function(
    x, 
    na.rm = TRUE
  ) {
    (
      x - 
        mean(
          x, 
          na.rm = na.rm
        )
    ) / 
      sd(
        x, 
        na.rm = na.rm
      )
  }

# Proportion of landscape that is inhospitable (WSI > wsi_cutoff)
inhospitable = 
  function(
    x
  ) {
    length(
      which(
        x > wsi_cutoff
      )
    ) /
      num_nodes
  }

# Proportion of landscape that is hospitable (WSI < wsi_cutoff)
hospitable = 
  function(
    x
  ) {
    length(
      which(
        x < wsi_cutoff
      )
    ) /
      num_nodes
  }

# Mean WSI within proportion of landscape that is hospitable
# (WSI < wsi_cutoff)
availmeanWSI = 
  function(
    x, 
    na.rm = TRUE
  ) {
    mean(
      x[
        x < wsi_cutoff
        ], 
      na.rm = na.rm
    )
  }

# Median WSI within proportion of landscape that is hospitable
# (WSI < wsi_cutoff)
availmedianWSI = 
  function(
    x, 
    na.rm = TRUE
  ) {
    median(
      x[
        x < wsi_cutoff
        ],
      na.rm = na.rm
    )
  }

# Mean WSI across landscape 
ALLmeanWSI = 
  function(
    x, 
    na.rm = TRUE
  ) {
    mean(
      x,
      na.rm = na.rm
    )
  }

# Median WSI across landscape
ALLmedianWSI = 
  function(
    x, 
    na.rm = TRUE
  ) {
    median(
      x,
      na.rm = na.rm
    )
  }

# Top 2% most populace nodes per day
node_populace = 
  function(
    x
  ) {
    which(
      x > 
        quantile(
          x,
          probs = 0.98
        )
    )
  }

# Center of mass of population per day
Pop_CoM = 
  function(
    wt,
    x,
    y,
    z
  ) {
    COGravity(
      x,
      y,
      z,
      wt
    )
  }

# Mask operation for arrays
fun.mask = 
  function(
    x
  ) {
    x[
      !arr_dummy_mask
      ] = NA
    return(
      x
    )
  }

# Rounding function
round_func = 
  function(
    x, 
    n
  ) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z * posneg
  }

# Temperature mean in available habitat function
availmeanTemp = 
  function(
    x, 
    na.rm = TRUE
  ) {
    mean(
      x[
        y < 
          wsi_cutoff
        ]
    )
  }

# Gamma movement functions for body condition transitions
gamma_x_func = 
  function(
    x
  ) {
    (
      (
        flight_dist - x
      ) * (
        x > 0
      )
    )
  }

t1_func = 
  function(
    x
  ) {
    x ^ (
      gamma_a - 1
    )
  }

# Integration function for skewed normal distribution
integrate_snorm = 
  Vectorize(
    function(
      l,
      h
    ) {
      integrate(
        function(x) {
          dsnorm(
            x,
            sd_body_mass
          )
        },
        lower = l,
        upper = h
      )$value
    }
  )

# Integration function for uniform distribution
lower_bound = 0
upper_bound = flight_dist

integrate_unif = 
  Vectorize(
    function(
      lower_limit,
      upper_limit
    ) {
      integrate(
        function(
          x
        ) 
        {
          dunif(
            x,
            lower_bound,
            upper_bound
          )
        },
        lower = lower_limit,
        upper = upper_limit
      )$value
    }
  )

# Axis label rounding function
scaleFUN = 
  function(
    x
  ) {
    sprintf(
      '%.1f', 
      x
    )
  }


#
# Define global landscape/movement parameters -----------------------------
# Define global geographical projection from LCC to WGS84
geo_prj =
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# How far birds can move within nodes, in meters
node_move = 1000 

# # Apparent velocity of air molecules "moving" past bird, proxy for 
# # flight velocity (Vmr), km/hr
# true_air_vel = 75

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Body mass in kg, from Owen and Cook (1977). They list a max mass of 
# 1.58 kg. To yield the proper mean, we need to allow the distribution
# to extend past their max (a low-probability tail up to 2 kg)
min_body_mass = 0.5
field_max_body_mass = 1.58
max_body_mass = 2
sd_body_mass = 
  (
    sqrt(
      15602
    ) * 0.97
  ) / 100

mean_body_mass = 1.2

body_mass_probs = 
  rsnorm(
    N_0,
    mean = mean_body_mass,
    sd = sd_body_mass
  )

body_mass_samp = 
  scales::rescale(
    body_mass_probs,
    to = 
      c(
        min_body_mass,
        max_body_mass
      )
  )

# Number of days in non-breeding period
fall_start = 
  as.Date(
    "2014-09-01", 
    tz = "UTC"
  )

spring_end = 
  strptime(
    "2015-05-31", 
    format = "%Y-%m-%d", 
    tz = "UTC"
  )

num_days = 
  as.numeric(
    difftime(
      as.POSIXct(
        spring_end
      ),
      as.POSIXct(
        fall_start, 
        tz = "UTC"
      ), 
      units = "days"
    )
  ) + 1

# Set day
days = 1:num_days

#Threshold for weather severity
wsi_cutoff = 7.5

# Number of land cover classes currently analyzed
num_land_cover_class = 3

# Chose year to pull data
year_from_record = 
  1957:2019

NonbreedingStart = 
  paste0(
    year_from_record,
    "-09-01"
  )

# Function to define temperature-dependent metabolic rate
TempDepMetRate =
  function(
    temperature,
    basalmetabolicrate,
    LowCritTemp,
    body_mass
  ) {
    (
      BMR_LCT(
        temperature,
        basalmetabolicrate,
        body_mass,
        LowCritTemp
      ) ^ 
        (
          temperature < LowCritTemp
        )
    ) *
      (
        basalmetabolicrate ^ (
          temperature >= LowCritTemp
        )
      )
  }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BMR below LCT, from Bicudo et al. (2010)
BMR_LCT = 
  function(
    temperature,
    basalmetabolicrate,
    body_mass,
    LowCritTemp
  ) {
    (
      2.82 * body_mass * (
        LowCritTemp - temperature
      )
    ) * kckJConv
  }

# Energy loss functions
# BMR in kJ per day, Heitmeyer 2010
BasalMetRate_func = 
  function(
    body_mass
  ) {
    422 * (
      (
        body_mass
      ) ^ 0.74
    )
  }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Energy gain functions (I THINK THIS VALUE REALLY MATTERS)
# Lindstrom 2003 shows maximum FDR for a ~1 kg non-passerine caps out 
# at 2% of lean body mass (0.3 - 2). He cites Hanson et al. (1990)
# who show a FDR in terms of % LBM for mallards of 0.49 (LBM = 0.9 kg). 
# He also states that FDR changes with body mass even among individuals, 
# so this static representation is a known simplifying assumption.
FDR_func = 
  function(
    body_mass
  ) {
    (
      (
        1.1 /
          100
      ) * body_mass 
    )
  }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Parameter to penalize birds staying in areas with severe weather. High 
# WSI inhibits access to forage material, so FDR should be 0.
fdr_penalty_param = 0


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

# Create node-specific temperature, WSI, and air density measures ---------
# TempData =
#   fread(
#     "NorAm_historical_air_temperature_node_info.txt"
#   )
# 
# TempDataFrame =
#   data.frame(
#     TempData
#   )
# 
# WSIData =
#   fread(
#     "NorAm_historical_wsi_node_info.txt"
#   )
# 
# WSIDataFrame =
#   data.frame(
#     WSIData
#   )
# 
# AirDensityData =
#   fread(
#     "NorAm_historical_airdensity_node_info.txt"
#   )
# 
# AirDensityDataFrame =
#   data.frame(
#     AirDensityData
#   )
# 
# data_columns =
#   colnames(
#     TempDataFrame
#   )
# 
# for(i in 1:(length(year_from_record)-1)){
#   #subset out a particular migration period (don't want to base this on
#   #calendar year, but rather on annual migratory cycle; i.e., fall one year
#   #to spring the next)
#   # 1972 was the most recent "average" (relative to 1951-1980) year,
#   # differing from the mean during that time span by 0.01 deg C.
#   target_year = year_from_record[i]
# 
#   nonbreeding_year =
#     c(
#       as.character(
#         target_year
#       ),
#       as.character(
#         target_year + 1
#       )
#     )
# 
#   target_dates =
#     c(
#       paste0(
#         nonbreeding_year[1],
#         c(
#           "_09_",
#           "_10_",
#           "_11_",
#           "_12_"
#         )
#       ),
#       paste0(
#         nonbreeding_year[2],
#         c(
#           "_01_",
#           "_02_",
#           "_03_",
#           "_04_",
#           "_05_"
#         )
#       )
#     )
# 
#   nonbreeding_period =
#     data_columns[
#       grepl(
#         target_dates[1],
#         data_columns
#       )
#       ]
# 
#   for(j in 2:length(target_dates)){
#     nonbreeding_period =
#       c(
#         nonbreeding_period,
#         data_columns[
#           grepl(
#             target_dates[j],
#             data_columns
#           )
#           ]
#       )
#   }
# 
#   nonbreeding_period =
#     nonbreeding_period[
#       !(
#         nonbreeding_period %like% "_02_29"
#       )
#       ]
# 
#   # Temperature
#   TempDT_NonbreedingPeriod =
#     data.table(
#       TempDataFrame[
#         ,
#         c(
#           'X_INDEX',
#           'Y_INDEX',
#           nonbreeding_period
#         )
#         ]
#     )
# 
#   setkey(
#     TempDT_NonbreedingPeriod,
#     X_INDEX,
#     Y_INDEX
#   )
# 
#   TempDT_NonbreedingPeriod[
#     ,
#     `:=`(
#       X_INDEX = NULL,
#       Y_INDEX = NULL
#     )
#   ]
# 
#   setnames(
#     TempDT_NonbreedingPeriod,
#     names(
#       TempDT_NonbreedingPeriod
#     ),
#     c(
#       # "Code",
#       paste0(
#         "Temperature_day_",
#         1:num_days
#       )
#     )
#   )
# 
#   # WSI
#   WSIDT_NonbreedingPeriod =
#     data.table(
#       WSIDataFrame[
#         ,
#         c(
#           'X_INDEX',
#           'Y_INDEX',
#           nonbreeding_period
#         )
#         ]
#     )
# 
#   setkey(
#     WSIDT_NonbreedingPeriod,
#     X_INDEX,
#     Y_INDEX
#   )
# 
#   WSIDT_NonbreedingPeriod[
#     ,
#     `:=`(
#       X_INDEX = NULL,
#       Y_INDEX = NULL
#     )
#   ]
# 
#   setnames(
#     WSIDT_NonbreedingPeriod,
#     names(
#       WSIDT_NonbreedingPeriod
#     ),
#     c(
#       # "Code",
#       paste0(
#         "WSI_day_",
#         1:num_days
#       )
#     )
#   )
# 
#   # Air Density
#   AirDensityDT_NonbreedingPeriod =
#     data.table(
#       AirDensityDataFrame[
#         ,
#         c(
#           'X_INDEX',
#           'Y_INDEX',
#           nonbreeding_period
#         )
#         ]
#     )
# 
#   setkey(
#     AirDensityDT_NonbreedingPeriod,
#     X_INDEX,
#     Y_INDEX
#   )
# 
#   AirDensityDT_NonbreedingPeriod[
#     ,
#     `:=`(
#       X_INDEX = NULL,
#       Y_INDEX = NULL
#     )
#   ]
# 
#   setnames(
#     AirDensityDT_NonbreedingPeriod,
#     names(
#       AirDensityDT_NonbreedingPeriod
#     ),
#     c(
#       # "Code",
#       paste0(
#         "AirDensity_day_",
#         1:num_days
#       )
#     )
#   )
# 
#   # Merge together, and add NODE_DATA
#   setkey(
#     NODE_DATA,
#     X_INDEX,
#     Y_INDEX
#   )
# 
#   NODE_DATA_Merge =
#     cbind(
#       NODE_DATA,
#       TempDT_NonbreedingPeriod,
#       WSIDT_NonbreedingPeriod,
#       AirDensityDT_NonbreedingPeriod
#     )
# 
#   fwrite(
#     NODE_DATA_Merge,
#     file =
#       paste0(
#         "NodeSpecifData_Temperature_WSI_AirDensity",
#         year_from_record[i],
#         ".txt"
#       )
#   )
# 
#   rm(
#     NODE_DATA_TEMPWSIAD,
#     NODE_DATA_Merge
#   )
# 
#   print(i)
# }
# 
# rm(
#   TempData,
#   TempDataFrame,
#   TempDT_NonbreedingPeriod
# )


#
# Generate body mass/body condition classes -------------------------------
n_bins = 21
bc_bins = c(1:20)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Body Mass is from Owen and Cook (1977)
# Body mass, g
body_mass_discrete = 
  seq(
    min_body_mass,
    max_body_mass,
    length.out = n_bins
  )

body_mass_discrete_dist = 
  as.vector(
    table(
      findInterval(
        body_mass_samp,
        seq(
          min_body_mass,
          field_max_body_mass,
          length.out = n_bins
        )
      )
    ) / N_0
  )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LCT is now variable, not set at 20 C

# Body condition table
# ~10% to 16% of body mass is composed of lipids at maximum body condition (Dabbert et al. 
# 1997, Boos et al. 2007). Lipids account for ~81% to 84% of metabolizable energy (Boos et 
# al. 2007). Not all lipids are available for metabolic processes (some retained for 
# other purposes, not detailed; Boos et al. 2002, 2007). We assume the ~16% to 19% of 
# metabolizable energy provided by sources other than lipids is used for processes other
# than flight (basal metabolic rate, reproductive organs, cellular replacement, etc.).
# Therefore we assume that all energy directed toward powered flight relies on lipids
# as its source exclusively, and not all of the ~10% to 16% of body mass comprised of 
# lipids is available for powered flight processes. We set this to 2% below the extremes 
# (8% and 14%).
mean_body_fat_prop = 0.11#0.14#0.08

BodyCondition_Table = 
  data.table(
    BC_bins = 
      c(
        1, 
        bc_bins + 1
      ),
    BodyMass_kg = body_mass_discrete,
    BodyMass_dist_discrete_day0 = body_mass_discrete_dist,
    BodyFat_kg =
      body_mass_discrete *
      seq(
        0,
        mean_body_fat_prop,
        length.out = n_bins
      ),
    BodyMass_BMR = 
      BasalMetRate_func(
        body_mass_discrete
      ),
    LowCritTemp = 
      47.2 * (
        (
          body_mass_discrete * 1000
        ) ^ -0.18
      )
  )

# Read in simulated flight velocity and cost data from varying mass,
# wing area, wing span, air density, and true-air velocity, from
# "flight_velocity_relationship v2.R" code.
flight_metrics_table =
  fread(
    'flight_vel_cost_tbl.csv'
  )

flight_vel_body_mass_rel =
  flight_metrics_table[
    ,
    lm(
      flight_velocity_kmhr ~ body_mass_kg
    )
    ]$coefficients

flight_cost_body_mass_rel =
  flight_metrics_table[
    ,
    lm(
      flight_cost_kgkm ~ body_mass_kg
    )
    ]$coefficients

# Use simulated data to generate relationship between body mass
# and flight velocity, and body mass and flight cost. Divide
# flight velocity ~ body mass slope by 2 to account for 
# lower equivalent air speed than true air speed (see Flight 1.25 from
# Pennycuick)
BodyCondition_Table[
  ,
  `:=`(
    FlightVel = 
      (
        BodyMass_kg *
          flight_vel_body_mass_rel[2] +
          flight_vel_body_mass_rel[1]
      ),
    FlightCost_kgkm =
      (
        BodyMass_kg *
          flight_cost_body_mass_rel[2] +
          flight_cost_body_mass_rel[1]
      )
  )
  ]

#Parameter for gamma distribution to determine movement probability
BodyCondition_Table[
  ,
  BodyFat_gamma_b := 
    flight_dist * 
    ZO_std(
      BodyFat_kg
    )
  ]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Survivorship peaks at an optimal body mass (assumed to be 
# the mean for the population). This is inline with Optimal Body Mass 
# theory (Lima 1986, Witter and Cuthill 1993, Guillemain et al. 2004)
# Bottom bin has 99.74% daily survivorship
# min_survival = 0.9974 
#Top bin has 99.997% daily survivorship
max_survival = 0.99984#0.99997

BodyCondition_Table[
  ,
  Survivorship := 
    -0.003 * (
      (
        BodyMass_kg - 1.6
      ) ^ 2
    ) + max_survival
  ]

# survivorship_plot = 
#   ggplot(
#     data = BodyCondition_Table, 
#     aes(
#       x = BC_bins, 
#       y = Survivorship
#     )
#   ) + 
#   geom_line() + 
#   xlab(
#     'Body Condition Class'
#   ) + 
#   ylab(
#     'Daily Survivorship'
#   ) + 
#   theme_bw() + 
#   theme(
#     axis.text = 
#       element_text(
#         size = 14
#       ), 
#     text = 
#       element_text(
#         size = 20
#       )
#   )
# 
# ggsave(
#   'bodycondition_survivorship.jpeg',
#   survivorship_plot,
#   width = 5,
#   height = 5
# )

# Possible flight distance by body condition (in km)
BodyCondition_Table[
  ,
  MaxPoss_flightdist :=
    range_limit *
    (
      BodyFat_kg / 
        FlightCost_kgkm
    )
  ]

# Calculate Pr(flying 'x' kms | body condition)
flightdist_distrib =
  sapply(
    BodyCondition_Table[
      ,
      MaxPoss_flightdist
      ],
    integrate_unif,
    l = 0
  )

# Calculate the cost (kg) of traveling each unit (node) of a flight
cummul_flightcost =
  data.table(
    dist_km =
      seq(
        0,
        flight_dist,
        node_size
      )
  )

bc_dep_cost_km = 
  BodyCondition_Table[
    ,
    paste0(
      'cost_km_BC_',
      BC_bins
    )
    ]

for(i in 1:nrow(BodyCondition_Table)) {
  cummul_flightcost[
    ,
    (bc_dep_cost_km[i]) :=
      c(
        0,
        cumsum(
          rep(
            BodyCondition_Table[
              i,
              FlightCost_kgkm * 
                node_size
              ],
            (
              flight_dist / 
                node_size
            )
          )
        )
      )
    ]
}

# Body condition-dependent flight distance probabilities
BC_distprob =
  BodyCondition_Table[
    ,
    ifelse(
      MaxPoss_flightdist == 0,
      0,
      MaxPoss_flightdist
    )
    ]

# Function to identify the possible movement distances for each body condition;
# since only birds with a probability of departing > 0 will have this applied
# to them, we remove consideration of the probability of moving 0 km
prob_dist_func =
  Vectorize(
    function(x) {
      cummul_flightcost[
        dist_km < x,
        dist_km[
          dist_km > 0
          ]
        ]
    }
  )

BC_ProbMoveDist =
  prob_dist_func(
    BC_distprob
  )

BC_gamma_x =
  lapply(
    BC_ProbMoveDist,
    gamma_x_func
  )

BC_gamma_b =
  flight_dist -
  BC_distprob

t1 =
  lapply(
    BC_gamma_x,
    t1_func
  )

t2 =
  vector(
    'list'
  )

for(i in 1:n_bins){
  t2[[i]] =
    BC_gamma_x[[i]] / BC_gamma_b[i]
}

t3 =
  gamma(
    gamma_a
  )

BC_gamma_moveprob = vector('list')
for(i in 1:n_bins){
  BC_gamma_moveprob[[i]] =
    (
      t1[[i]] * (
        exp(1) ^ (-t2[[i]]) / (
          BC_gamma_b[i] ^ gamma_a * t3
        )
      )
    )
}

# Convert ragged BC_gamma_moveprob list to data.frame
rowMax =
  max(
    sapply(
      BC_gamma_moveprob,
      length
    )
  )

BC_gamma_moveprob_dt =
  data.table(
    do.call(
      cbind,
      lapply(
        BC_gamma_moveprob,
        function(x){
          length(
            x
          ) = 
            rowMax
          x
        }
      )
    )
  )

setnames(
  BC_gamma_moveprob_dt,
  names(BC_gamma_moveprob_dt),
  c(
    paste0(
      'BC_',
      1:n_bins
    )
  )
)

# Normalize the gamma values for relativized movement probabilities
BC_gamma_moveprob_dt = NAFunc(BC_gamma_moveprob_dt)

BC_gamma_moveprob_mtx = apply(BC_gamma_moveprob_dt,2,norm_func)

BC_gamma_moveprob_mtx = NAFunc(BC_gamma_moveprob_mtx)

# Add the possible movement distances for reference
BC_gamma_moveprob_dt = data.table(BC_gamma_moveprob_mtx)

BC_gamma_moveprob_dt[
  ,
  PossibleDistances :=
    BC_ProbMoveDist[[n_bins]]
  ]

setcolorder(
  BC_gamma_moveprob_dt,
  c(
    n_bins + 1,
    1:n_bins
  )
)

# Calculate the number of bins by which birds will be reduced for flying
# different distances
BC_transition_kg_ls = list()
for(k in n_bins:2){
  BC_transition_kg_ls[[k]] = 
    BodyCondition_Table[
      k,
      BodyFat_kg
      ] - 
    BodyCondition_Table[
      1:(k-1),
      BodyFat_kg
      ]
}

BC_transition_kg_ls_maxln = 
  max(
    sapply(
      BC_transition_kg_ls,
      length
    )
  )

BC_transition_kg_mtx = 
  do.call(
    cbind,
    lapply(
      BC_transition_kg_ls,
      function(
        x
      ) {
        length(
          x
        ) = 
          BC_transition_kg_ls_maxln
        x
      }
    )
  )

BC_transition_kg_mtx = 
  NAFunc(
    BC_transition_kg_mtx
  )

BC_transition_kg_mtx = 
  rbind(
    BC_transition_kg_mtx,
    rep(
      0,
      n_bins - 1
    )
  )

BC_transition_kg_mtx = 
  cbind(
    rep(
      0,
      n_bins
    ),
    BC_transition_kg_mtx
  )

BC_transitions = 
  sort(
    BC_transition_kg_mtx[
      ,
      n_bins
      ][
        which(
          BC_transition_kg_mtx[
            ,
            n_bins
            ] > 0
        )
        ]
  )

BC_decrement_total = vector('list')
for(i in 1:nrow(BodyCondition_Table)){
  BC_decrement_total[[i]] = 
    rep(
      0,
      nrow(
        cummul_flightcost
      )
    )
  
  BC_decrement = vector()
  for(j in 1:length(BC_transitions)){
    BC_decrement =
      ifelse(
        cummul_flightcost[
          ,
          (bc_dep_cost_km[i]),
          with = FALSE
          ] < 
          BC_transitions[j],
        0,
        1
      )
    
    BC_decrement_total[[i]] = BC_decrement + BC_decrement_total[[i]]
  }
}

BC_decrementvalue = 
  BodyCondition_Table[
    ,
    paste0(
      'BC_decrementvalue_',
      BC_bins
    )
    ]

cummul_flightcost[
  ,
  (BC_decrementvalue) := 
    BC_decrement_total
  ]   

# Add body condition decrement per km-bin traveled to movement probability matrix
for(i in 1:nrow(BodyCondition_Table)){
  BC_gamma_moveprob_decrements = 
    cummul_flightcost[
      cummul_flightcost[
        ,
        dist_km
        ] %in% 
        BC_gamma_moveprob_dt[
          ,
          PossibleDistances
          ],
      (BC_decrementvalue[i]),
      with = FALSE
      ]
  
  BC_gamma_moveprob_dt[
    ,
    (BC_decrementvalue[i]) := 
      BC_gamma_moveprob_decrements
    ]
}

BC_gamma_moveprob_dt[
  ,
  PossibleDistances := NULL
  ]

sum_cols = 
  paste0(
    'BC_',
    1:n_bins
  )

BC_transitionprobs = list()
for(i in 1:nrow(BodyCondition_Table)) {
  BC_transitionprobs[[i]] = 
    BC_gamma_moveprob_dt[
      ,
      lapply(
        .SD,
        sum
      ),
      .SDcols = sum_cols,
      by = 
        eval(
          BC_decrementvalue[i]
        )
      ]
}

# Calculate number of individuals entering each body condition using
# a transition matrix. Need to convert from data table with 'from' in 
# columns and 'to' as: 
#           column number - body condition decrement value of row
# to transition matrix with 'from' in columns and 'to' in rows
# fill in matrix based on row and column id
BC_Transition_mtx=
  matrix(
    NA,
    nrow = 21,
    ncol = 21
  )

for(i in 1:ncol(BC_Transition_mtx)) {
  BC_transitionprobs_id = 
    rbind(
      BC_transitionprobs[[i]],
      matrix(
        0,
        nrow = 
          n_bins - 
          nrow(
            BC_transitionprobs[[i]]), 
        ncol = 
          ncol(
            BC_transitionprobs[[i]])), 
      use.names = FALSE
    )
  
  BC_transitionprobs_id[
    ,
    (
      BC_decrementvalue[i]
    ) :=
      0:(n_bins - 1)
    ]
  
  for(j in 1:nrow(BC_Transition_mtx)) {
    BC_Transition_mtx[j,i] =
      as.numeric(
        BC_transitionprobs_id[
          get(
            BC_decrementvalue[i]
          ) == i - j,
          (i + 1),
          with = FALSE
          ]
      )
  }
}

BC_Transition_mtx = NAFunc(BC_Transition_mtx)


#
# Define parameters and rules for movement outside of loop ----------------
# Probability of departure component related to WSI
# The relation below comes from Schummer et al. 2010
WSIProb = 
  function(
    x
  ) {
    -0.6282 + (
      0.0547 * x
    ) + (
      0.0046 * (
        x ^ 2
      )
    )
  }

# Probability of departure component related to distance to nearest 
# breeding node
breedingdistances = 
  NODE_DATA[
    ,
    DistBreed
    ]

n_breedingdistances = 
  max(
    breedingdistances
  )

# Mortality probability matrix
ProbMort = 
  1 - 
  BodyCondition_Table[
    ,
    Survivorship
    ]

ProbMort_mtx = 
  matrix(
    rep(
      ProbMort,
      nrow(
        NODE_DATA
      )
    ),
    ncol = n_bins,
    byrow = TRUE
  )

# Natural forage decay rates by landcover type
# Convert this to a matrix with dimensions dependent upon the number
# of landcover types used
Shoreline_decay = 0.9998
Crops_decay = 0.997
WoodyWetlands_decay = 0.9965
HerbWetlands_decay = 0.991

# Establish inflection point--date on which importance of landscape
# variables switch directions
set.seed(7)
inflection = 
  floor(
    runif(
      1,
      (num_days / 2) - 15,
      (num_days / 2) + 15
    )
  )

# Generate sequence of weights
inflection_seq = 
  seq(
    0,
    1,
    length.out = num_days
  )

# Node departure probability

# Weights for WSI, body condition, and distance to breeding node for
# probability of departure functions
weight_WSI = 0.4
weight_DB = 0.2
weight_DG = 0.3
# weight_AD = 0.25

DepartScheme = 
  data.table(
    DayofNonBreedingPeriod = 
      1:num_days,
    WSIComponent = 
      (
        - (
          weight_WSI * 2
        ) * 
          (
            (
              inflection_seq - 
                (
                  inflection / 
                    num_days
                )
            ) ^ 2
          )
      ) + 
      weight_WSI,
    DBComponent = 
      (
        - (
          weight_DB * 2
        ) * 
          (
            (
              inflection_seq - 
                (
                  inflection / 
                    num_days
                )
            ) ^ 2
          )
      ) + 
      weight_DB,
    DGComponent = 
      (
        - (
          weight_DG * 2
        ) * 
          (
            (
              inflection_seq - 
                (
                  inflection / 
                    num_days
                )
            ) ^ 2
          )
      ) + 
      weight_DG,
    BCComponent =
      1 - (
        (
          (
            - (
              weight_WSI * 2
            ) * 
              (
                (
                  inflection_seq - 
                    (
                      inflection / 
                        num_days
                    )
                ) ^ 2
              )
          ) + 
            weight_WSI
        ) +
          (
            (
              - (
                weight_DB * 2
              ) * 
                (
                  (
                    inflection_seq - 
                      (
                        inflection / 
                          num_days
                      )
                  ) ^ 2
                )
            ) + 
              weight_DB
          ) +
          (
            (
              - (
                weight_DG * 2
              ) * 
                (
                  (
                    inflection_seq - 
                      (
                        inflection / 
                          num_days
                      )
                  ) ^ 2
                )
            ) + 
              weight_DG
          )
      )
  )

# DepartScheme_Time =
#   melt(
#     DepartScheme,
#     id.vars =
#       'DayofNonBreedingPeriod'
#   )
# 
# DepartScheme_Plot =
#   ggplot(
#     data = DepartScheme_Time,
#     aes(
#       x = DayofNonBreedingPeriod,
#       y = value,
#       linetype = variable
#     )
#   ) +
#   geom_line(
#     lwd = 1.25
#   ) +
#   scale_linetype_manual(
#     values =
#       c(
#         'solid',
#         'longdash',
#         'dotdash',
#         'dotted'
#       ),
#     guide = 'none'
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y =
#         1 - (
#           weight_WSI + weight_DB + weight_DG
#         ) - 0.03
#     ),
#     label = 'Body Condition',
#     size = 5
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = weight_WSI + 0.03
#     ),
#     label = 'WSI',
#     size = 5
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = weight_DG + 0.03
#     ),
#     label = 'Daily Gain',
#     size = 5
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = weight_DB + 0.03
#     ),
#     label = 'Breeding',
#     size = 5
#   ) +
#   ylim(
#     0,
#     1
#   ) +
#   labs(
#     x = 'Day of Nonbreeding Period',
#     y = 'Relative weight'
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text =
#       element_text(
#         face = 'bold',
#         size = 16
#       ),
#     axis.title =
#       element_text(
#         face = 'bold',
#         size = 20
#       )
#   )
# 
# ggsave(
#   'CobbDouglasWeights_Departure.jpeg',
#   DepartScheme_Plot,
#   width = 5,
#   height = 5
# )

# Node attraction probability
# Weights for Forage, Roosting habitat, Distance to nearest breeding 
# node, Air density, and default Gamma movement distance probability for
# probability of attraction functions
weight_F = 0.45
weight_AD = 0.3
weight_R = 0.2
weight_G = 0.05

AttractScheme = 
  data.table(
    DayofNonBreedingPeriod = 
      1:num_days,
    ForageAvailComponent = 
      (
        - (
          weight_F * 2
        ) * 
          (
            (
              inflection_seq - 
                (
                  inflection / 
                    num_days
                )
            ) ^ 2
          )
      ) + 
      weight_F,
    AirDensityComponent = 
      (
        - (
          weight_AD * 2
        ) * 
          (
            (
              inflection_seq - 
                (
                  inflection / 
                    num_days
                )
            ) ^ 2
          )
      ) + 
      weight_AD,
    RoostingHabComponent = 
      (
        - (
          weight_R * 2
        ) * 
          (
            (
              inflection_seq - 
                (
                  inflection / 
                    num_days
                )
            ) ^ 2
          )
      ) + 
      weight_R,
    BreedingDistComponent = 
      1 - (
        (
          (
            - (
              weight_F * 2
            ) * 
              (
                (
                  inflection_seq - 
                    (
                      inflection / 
                        num_days
                    )
                ) ^ 2
              )
          ) + 
            weight_F
        ) +
          (
            (
              - (
                weight_AD * 2
              ) * 
                (
                  (
                    inflection_seq - 
                      (
                        inflection / 
                          num_days
                      )
                  ) ^ 2
                )
            ) + 
              weight_AD
          ) +
          (
            (
              - (
                weight_R * 2
              ) * 
                (
                  (
                    inflection_seq - 
                      (
                        inflection / 
                          num_days
                      )
                  ) ^ 2
                )
            ) + 
              weight_R
          ) +
          (
            weight_G
          )
      ),
    GammaMoveComponent = 
      weight_G
  )

# AttractScheme_Time =
#   melt(
#     AttractScheme,
#     id.vars = 
#       'DayofNonBreedingPeriod'
#   )
# 
# AttractScheme_Plot =
#   ggplot(
#     data = AttractScheme_Time,
#     aes(
#       x = DayofNonBreedingPeriod,
#       y = value,
#       linetype = variable
#     )
#   ) +
#   geom_line(
#     lwd = 1.25
#   ) +
#   scale_linetype_manual(
#     values =
#       c(
#         'solid',
#         'longdash',
#         'dotdash',
#         'dotted',
#         'twodash'
#       ),
#     guide = 'none'
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = 
#         (
#           1 - 
#             sum(
#               weight_F,
#               weight_AD,
#               weight_R,
#               weight_G
#             )
#         ) - 0.03
#     ),
#     label = 'Breeding',
#     size = 5
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = weight_G + 0.05
#     ),
#     label = 'Gamma',
#     size = 5
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = weight_F + 0.03
#     ),
#     label = 'Forage',
#     size = 5
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = weight_R + 0.03
#     ),
#     label = 'Roosting',
#     size = 5
#   ) +
#   geom_text(
#     aes(
#       x = inflection,
#       y = weight_AD + 0.03
#     ),
#     label = 'Air Density',
#     size = 5
#   ) +
#   ylim(
#     (
#       1 - 
#         sum(
#           weight_F,
#           weight_AD,
#           weight_R,
#           weight_G
#         )
#     ) - 0.06,
#     1
#   ) +
#   labs(
#     x = 'Day of Nonbreeding Period',
#     y = 'Relative weight'
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text =
#       element_text(
#         face = 'bold',
#         size = 16
#       ),
#     axis.title =
#       element_text(
#         face = 'bold',
#         size = 20
#       )
#   )
# 
# ggsave(
#   'CobbDouglasWeights_Atraction.jpeg',
#   AttractScheme_Plot,
#   width = 5,
#   height = 5
# )


#
# Impose movement rules ---------------------------------------------------
gc()

rm(NODE_DATA)

deadbirds = vector()
daily_mort = list()
COMPLETE_DATA = NULL

Sys.time()
start_timer = Sys.time()

n_cores = detectCores() - 1
start_cluster = makeCluster(n_cores)
registerDoParallel(start_cluster)

foreach(
  j = 1:(
    length(
      year_from_record
    ) - 1
  ),
  .packages = 'data.table',
  .verbose = TRUE
) %dopar% {

  COMPLETE_DATA =
    fread(
      paste0(
        "NodeSpecifData_Temperature_WSI_AirDensity",
        year_from_record[j],
        ".txt"
      )
    )

  # Reduce forage material based on mass- and temperature-dependent
  # daily gain for each body mass/condition class and landcover-specific
  # decay rates, with place holders for decay rates
  for(i in 1:length(days)) {
    DayAbund =
      as.matrix(
        COMPLETE_DATA[
          ,
          paste0(
            "N_Abund_",
            i-1
          ),
          with = FALSE
          ]
      )

    DayTemper =
      as.matrix(
        COMPLETE_DATA[
          ,
          paste0(
            "Temperature_day_",
            i
          ),
          with = FALSE
          ]
      )

    DayAirDensity =
      as.matrix(
        COMPLETE_DATA[
          ,
          paste0(
            "AirDensity_day_",
            i
          ),
          with = FALSE
          ]
      )

    DailyGain =
      paste0(
        "DailyGain_day_",
        i
      )

    ReducedForage =
      paste0(
        "Forage_day_",
        i
      )

    BM_PopDist =
      as.vector(
        as.matrix(
          BodyCondition_Table[
            ,
            paste0(
              "BodyMass_dist_discrete_day",
              i-1
            ),
            with = FALSE
            ]
        )
      )

    # Daily gain is calculated with the inputs:
    #   HumanDisturb (Degree of human disturbance--urbanization)
    #   HarvestDisturb (Degree of human disturbance--hunter harvest)
    #   FDR (Fuel deposition rate, kJ ingested per day)
    #   TD_BMR (Temperature dependent BMR, covnert from kJ to kg)
    # DG = ((1 - HumanDisturb) * (1 - HarvestDisturb) * FDR) - TD_BMR

    # Calculate the number of birds in each body condition class
    # in each node on the day
    BC_Node_Abund =
      matrix(
        (
          DayAbund %o% BM_PopDist
        ),
        ncol = n_bins
      )

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++                   CHANGED           ++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Changed to make sure each of the below objects have the same dimensions.
    # Included WSI-related penalty for FDR.

    # Calculate the basal metabolic rate for each body
    # condition class in each node on the day
    BMR =
      BodyCondition_Table[
        ,
        matrix(
          rep(
            BodyMass_BMR,
            each = num_nodes
          ),
          ncol = n_bins
        )
        ]

    # Calculate the body mass for each body condition class
    # in each node on the day
    BC_BM_node =
      BodyCondition_Table[
        ,
        matrix(
          rep(
            BodyMass_kg,
            each = num_nodes
          ),
          ncol = n_bins
        )
        ]

    # Calculate the disturbance-dependent FDR (kg)
    WSI_mtx =
      matrix(
        rep(
          unlist(
            COMPLETE_DATA[
              ,
              paste0("WSI_day_",i),
              with = FALSE
              ]
          ),
          times = n_bins
        ),
        ncol = n_bins
      )

    Disturb_mtx =
      COMPLETE_DATA[
        ,
        matrix(
          rep(
            (
              (
                1 - HumanDisturb
              ) * (
                1 - HarvestDisturb
              )
            ),n_bins
          ),
          ncol = n_bins
        )
        ]

    # If WSI is greater than the threshold, FDR is set to 0; otherwise
    # it remains as defined.
    FDR_penalty =
      ifelse(
        WSI_mtx > wsi_cutoff,
        fdr_penalty_param,
        1
      )

    Disturb_FDR =
      Disturb_mtx *
      FDR_func(
        BC_BM_node
      ) *
      FDR_penalty

    rm(
      Disturb_mtx
    )

    # Calculate the temperature-dependent basal metabolic
    # rate for each body condition class in each node on the day
    LCT =
      matrix(
        rep(
          BodyCondition_Table[
            ,
            LowCritTemp
            ],
          each = num_nodes
        ),
        ncol = n_bins
      )

    DayTemper_BC =
      matrix(
        rep(
          DayTemper,
          n_bins
        ),
        ncol = n_bins
      )

    TD_BMR =
      TempDepMetRate(
        temperature = DayTemper_BC,
        basalmetabolicrate = BMR,
        LowCritTemp = LCT,
        body_mass = BC_BM_node
      ) /
      39700

    rm(
      DayTemper_BC
    )

    # Calculate the daily gain in each node on the day
    Node_Disturb_FDR =
      rowSums(
        Disturb_FDR
      )

    COMPLETE_DATA[
      ,
      (
        DailyGain
      ) :=
        Node_Disturb_FDR
      ]

    # Calculate the forage material reduction
    COMPLETE_DATA[
      ,
      (ReducedForage):=
        (
          (
            (
              ForageShoreline * (
                Shoreline_decay ^ i
              )
            ) +
              (
                ForageCrops * (
                  Crops_decay ^ i
                )
              ) +
              (
                ForageWoodyWetlands *
                  (
                    WoodyWetlands_decay ^ i
                  )
              ) +
              (
                ForageHerbWetlands *
                  (
                    HerbWetlands_decay ^ i
                  )
              )
          )
        ) -
        Node_Disturb_FDR
      ]

    NormReducedForage =
      paste0(
        "NormForage_day_",
        i
      )

    incols = ReducedForage

    COMPLETE_DATA[
      ,
      (NormReducedForage) :=
        lapply(
          .SD,
          norm_func
        ),
      .SDcols = incols
      ]

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++                   CHANGED           ++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Simplify BC transition function

    # Calculate kg gained per body condition to redistribute
    # individuals among body condition classes after foraging
    DG_BC_Node =
      Disturb_FDR - TD_BMR

    DG_BC_Node_Transition_ls =
      vector(
        'list'
      )

    for(
      k in 1:n_bins
    ) {
      BF_BC_delta =
        BodyCondition_Table[
          k,
          BodyFat_kg
        ] +
        DG_BC_Node[
          ,
          k
        ]

      DG_BC_Node_Transition_ls[[k]] =
        findInterval(
          BF_BC_delta,
          BodyCondition_Table[
            ,
            BodyFat_kg
            ]
        )
    }

    DG_BC_Node_Transition =
      matrix(
        unlist(
          DG_BC_Node_Transition_ls
        ),
        ncol = n_bins
      )

    DG_BC_Node_Transition_dt =
      data.table(
        node_num = 1:num_nodes,
        DG_BC_Node_Transition
      )

    setnames(
      DG_BC_Node_Transition_dt,
      names(
        DG_BC_Node_Transition_dt
      ),
      gsub(
        'V',
        '',
        names(
          DG_BC_Node_Transition_dt
        )
      )
    )

    DG_BC_Node_Transition_melt =
      data.table::melt(
        DG_BC_Node_Transition_dt,
        id.vars = 'node_num',
        value.name = 'BC_transition'
      )

    BC_Node_Propor =
      BC_Node_Abund / N_0

    BC_Node_Propor_dt =
      data.table(
        node_num = 1:num_nodes,
        BC_Node_Propor
      )

    setnames(
      BC_Node_Propor_dt,
      names(
        BC_Node_Propor_dt
      ),
      gsub(
        'V',
        '',
        names(
          BC_Node_Propor_dt
        )
      )
    )

    BC_Node_Propor_melt =
      data.table::melt(
        BC_Node_Propor_dt,
        id.vars = 'node_num',
        variable.name = 'BC_bin',
        value.name = 'BC_propor'
      )

    BC_Node_Propor_melt[
      ,
      BC_transition :=
        DG_BC_Node_Transition_melt[
          ,
          BC_transition
          ]
      ]

    BC_Node_Propor_melt[
      BC_transition == 0,
      BC_transition := 1
    ]

    BC_Node_Propor_melt[
      ,
      BC_bin_transition :=
        paste0(
          BC_bin,
          '_',
          BC_transition
        )
    ]

    DG_BC_Transition_dt =
      BC_Node_Propor_melt[
        ,
        sum(
          BC_propor
        ),
        by =
          .(
            BC_bin,
            BC_transition
          )
        ][
      ,
      sum(
        V1
      ),
      by = BC_transition
      ]

    DG_BC_Transition_dt[
      ,
      V1 :=
        V1 /
        sum(
          V1
        )
      ]

    BC_bin_missing =
      which(
        !1:n_bins %in%
          DG_BC_Transition_dt[
            ,
            unique(
              BC_transition
            )
            ]
      )

    DG_BC_Transition_dt =
      rbind(
        DG_BC_Transition_dt,
        cbind(
          BC_bin_missing,
          rep(
            0,
            length(
              BC_bin_missing
            )
          )
        ),
        use.names = FALSE
      )

    setkey(
      DG_BC_Transition_dt,
      BC_transition
    )

    BC_Propor_New =
      paste0(
        "BodyMass_dist_discrete_day",
        i
      )

    BodyCondition_Table[
      ,
      (
        BC_Propor_New
      ) :=
        DG_BC_Transition_dt[
          ,
          V1
          ]
      ]

    BC_Abund =
      matrix(
        (
          DayAbund %o% BM_PopDist
        ),
        ncol = n_bins
      )

    rm(
      DailyGain,
      ReducedForage,
      DG_BC_Node,
      DG_BC_Node_Transition_ls,
      DG_BC_Node_Transition,
      DG_BC_Node_Transition_dt,
      DG_BC_Node_Transition_melt,
      BC_Node_Propor,
      BC_Node_Propor_dt,
      BC_Node_Propor_melt,
      DG_BC_Transition_dt,
      BC_Propor_New
    )

    # Determine WSI, Distance to breeding grounds, and Body condition-
    # dependent probability of departure
    # Trying to implement air density-related departure does not seem
    # practical: no empirical relation of air pressure on relative
    # abundance of mallards. In fact, see:
    # Pearse, A. T. 2007. Design, evaluation and applications of an aerial
    # survey to estimate abundance of wintering waterfowl in Mississippi.
    # Dissertation, Mississippi State University, Mississippi State, USA.

    # Probability of departure component related to body condition
    ProbDepart_BC =
      c(
        0,
        (
          DepartScheme[
            i,
            BCComponent
          ] * (
            (
              c(
                1:n_bins
              ) ^ 8
            ) /
              (
                (
                  c(
                    1:n_bins
                  ) ^ 8
                ) +
                  (
                    BodyCondition_Table[
                      ,
                      which.max(
                        BodyMass_dist_discrete_day0
                        )
                      ] ^ 8
                  )
              )
          )
        )[-1]
      )

    ProbDepart_BC_mtx =
      matrix(
        rep(
          ProbDepart_BC,
          num_nodes
        ),
        ncol = n_bins,
        byrow = TRUE
      )

    # Probability of departure component related to distance to nearest
    # breeding node
    ProbDepart_DB =
      DepartScheme[
        i,
        DBComponent
        ] * (
        breedingdistances ^ 5
      ) /
      (
        (
          breedingdistances ^ 5
        ) +
          (
            (
              n_breedingdistances - (
                n_breedingdistances
              ) * 0.33
            )
          ) ^ 5
      )

    ProbDepart_DB_mtx =
      matrix(
        rep(
          ProbDepart_DB,
          n_bins
        ),
        ncol = n_bins,
        byrow = TRUE
      )

    # Probability of departure component related to weather severity index
    WSIProbDepart =
      as.matrix(
        COMPLETE_DATA[
          ,
          paste0("WSI_day_",i),
          with = FALSE
          ]
      )

    ProbDepart_WSI =
      as.vector(
        NAFunc(
          (
            DepartScheme[
              i,
              WSIComponent
              ] *
              (
                (
                  (
                    WSIProb(
                      WSIProbDepart
                    ) +
                      wsi_cutoff
                  ) ^ 3
                ) /
                  (
                    (
                      (
                        WSIProb(
                          WSIProbDepart
                        ) +
                          wsi_cutoff
                      ) ^ 3
                    ) +
                      (
                        wsi_cutoff ^ 3
                      )
                  )
              ) *
              (
                WSIProbDepart > wsi_cutoff
              )
          )
        )
      )

    ProbDepart_WSI_mtx =
      matrix(
        rep(
          ProbDepart_WSI,
          n_bins
        ),
        ncol = n_bins
      )

    # Probability of departure component related to daily gain 
    # (fuel deposition rate)
    ProbDepart_DG_mtx =
      DepartScheme[
        i,
        DGComponent
      ] *
      ZO_std(
        (
          1 - (
            (
              Disturb_FDR + 1
            ) / 5
          )
        )
      )

    ProbDepart =
      ProbDepart_BC_mtx +
      ProbDepart_WSI_mtx +
      ProbDepart_DG_mtx +
      ProbDepart_DB_mtx

    AbundDepart_vec =
      rowSums(
        matrix(
          (
            DayAbund %o% BM_PopDist
          ),
          ncol = n_bins
        ) * ProbDepart
      )

    rm(
      ProbDepart,
      ProbDepart_WSI,
      ProbDepart_WSI_mtx,
      ProbDepart_DG_mtx
    )

    # Use Cobb-Douglas function to calculate node-specific
    # attractiveness.
    NodeAttract =
      paste0(
        "NodeAttract_day_",
        i
      )

    NormAirDens =
      norm_func(
        as.matrix(
          COMPLETE_DATA[
            ,
            paste0(
              "AirDensity_day_",
              i
            ),
            with = FALSE
            ]
        )
      )

    COMPLETE_DATA[
      ,
      (NodeAttract) :=
        as.numeric(
          (
            WSIProbDepart <= wsi_cutoff
          ) *
            (
              Re(
                as.complex(
                  as.vector(
                    unlist(
                      COMPLETE_DATA[
                        ,
                        (
                          NormReducedForage
                        ),
                        with = FALSE
                        ]
                    )
                  )
                ) ^
                  AttractScheme[
                    i,
                    ForageAvailComponent
                    ]
              )
            ) *
            (
              NormRoostShoreline ^
                AttractScheme[
                  i,
                  RoostingHabComponent
                  ]
            ) *
            (
              NormDistBreed ^
                AttractScheme[
                  i,
                  BreedingDistComponent
                  ]
            ) *
            (
              NormAirDens ^
                AttractScheme[
                  i,
                  AirDensityComponent
                  ]
            ) *
            (
              GammaMoveProb ^
                AttractScheme[
                  i,
                  GammaMoveComponent
                  ]
            )
        )
      ]

    COMPLETE_DATA[
      ,
      (NodeAttract) :=
        lapply(
          .SD,
          NAFunc
        ),
      .SDcols = NodeAttract
      ]

    COMPLETE_DATA[
      ,
      (NodeAttract) :=
        lapply(
          .SD,
          norm_func
        ),
      .SDcols = NodeAttract
      ]

    rm(
      WSIProbDepart
    )

    # Determine next day's abundance based on probability of staying
    # and departing.
    NextAbund =
      paste0(
        "N_Abund_",
        i
      )

    NodeAttract =
      as.matrix(
        COMPLETE_DATA[
          ,
          paste0(
            "NodeAttract_day_",
            days[i]
          ),
          with=FALSE
          ]
      )

    gc()

    CurrAbund =
      (
        DayAbund - AbundDepart_vec
      ) +
      (
        rowSums(
          NodeAttract %o% AbundDepart_vec,
          na.rm = TRUE
        )
      )

    gc()

    # Redistribute departing population among body condition classes based
    # on movement among nodes
    BC_AbundDepart_vec =
      as.vector(
        colSums(
          AbundDepart_vec %o% BM_PopDist
        )
      )

    BC_AbundStay_vec =
      colSums(
        BC_Abund
      ) - BC_AbundDepart_vec

    BC_AbundDepart_mtx =
      matrix(
        rep(
          BC_AbundDepart_vec,
          each = n_bins
        ),
        ncol = n_bins
      )

    BC_AbundArrive_mtx =
      BC_Transition_mtx *
      BC_AbundDepart_mtx

    BC_AbundArrive_vec =
      rowSums(
        BC_AbundArrive_mtx
      )

    BC_Abund_Move =
      norm_func(
        BC_AbundStay_vec + BC_AbundArrive_vec
      )

    BM_PopDist =
      as.vector(
        as.matrix(
          BodyCondition_Table[
            ,
            paste0(
              "BodyMass_dist_discrete_day",
              i
            ),
            with = FALSE
            ]
        )
      )

    BC_Abund_Next =
      (
        BM_PopDist + BC_Abund_Move
      ) / 2

    BC_Propor_New =
      paste0(
        "BodyMass_dist_discrete_day",
        i
      )

    BodyCondition_Table[
      ,
      (
        BC_Propor_New
      ) :=
        BC_Abund_Next
      ]

    BM_PopDist =
      as.vector(
        as.matrix(
          BodyCondition_Table[
            ,
            paste0(
              "BodyMass_dist_discrete_day",
              i
            ),
            with=FALSE
            ]
        )
      )

    # Incur mortality
    ProporDie =
      rowSums(
        floor(
          matrix(
            (
              CurrAbund %o% BM_PopDist
            ),
            ncol = n_bins
          )
        ) * ProbMort_mtx
      )

    # Calulcate remaining abundance after mortality
    COMPLETE_DATA[
      ,
      (NextAbund) :=
        CurrAbund - ProporDie
      ]

    deadbirds[i] =
      sum(CurrAbund) -
      sum(
        COMPLETE_DATA[
          ,
          NextAbund,
          with = FALSE
          ]
      )

    # Remove objects from environment to optimize memory use
    rm(
      DayAbund,
      NextAbund,
      NodeAttract,
      AbundDepart_vec,
      CurrAbund,
      ProporDie,
      NormReducedForage,
      BC_Abund
    )

    print((i / num_days) * 100)
  }

  # Save COMPLETE_DATA and BodyCondition_Table elements and remove
  # to reduce memory usage
  fwrite(
    COMPLETE_DATA,
    file =
      paste0(
        year_from_record[j],
        "_Prediction.txt"
      )
  )

  fwrite(
    BodyCondition_Table,
    file =
      paste0(
        year_from_record[j],
        "_BCTable.txt"
      )
  )
  rm(COMPLETE_DATA)

  print(j)
}

stopCluster(start_cluster)

rm(
  ProbDepart_BC_mtx,
  ProbDepart_DB_mtx,
  ProbMort_mtx,
  DG_BC_Node_Transition,
  DG_BC_Node_Transition_mtx,
  breedingdistances,
  ProbDepart_DB,
  BC_DepartAbund_mtx,
  AttractScheme,
  BC_DepartAbundPropor_dist,
  BC_dist_seq_ls,
  BC_dist_seq_ls_maxln,
  BC_dist_seq_mtx,
  BC_distprob,
  BC_Abund_Gain,
  BC_DepartAbundDist,
  BC_DepartAbundDist_diags,
  BC_DepartAbundDist_mtx
)

Sys.time()
end_timer = Sys.time()
total_time = end_timer - start_timer
total_time


#
# Read in predictive migration dynamics for analysis ----------------------
availhabitat = MeanWSI_availhabitat = MedianWSI_availhabitat =
  MeanWSI_allhabitat = MedianWSI_allhabitat = daily_mort = vector('list')
annmean_availhabitat = annsd_availhabitat = annmin_availhabitat =
  median_forage = Minimum_annualtemp = vector()

# Particular years to evaluate: c(1957,1962, 1989, 1999, 2012)

for(
  i in 1:(
    length(
      year_from_record
    ) - 1
  )
) {
  COMPLETE_DATA =
    fread(
      paste0(
        getwd(),
        '/',
        year_from_record[i],
        "_Prediction.txt"
      )
    )

  # Daily mortality over non-breeding season
  mort_cols =
    names(
      COMPLETE_DATA
    )[
      names(
        COMPLETE_DATA
      ) %like%
        'N_Abund_'
      ]

  daily_mort[[i]] =
    as.vector(
      COMPLETE_DATA[
        ,
        colSums(
          .SD
        ),
        .SDcols =
          mort_cols
        ]
    )[
      1:(
        length(
          mort_cols
        ) - 1
      )
      ] -
    as.vector(
      COMPLETE_DATA[
        ,
        colSums(
          .SD
        ),
        .SDcols =
          mort_cols
        ]
    )[
      2:length(
        mort_cols
      )
      ]

  # Availability of habitat across non-breeding period (based on WSI)
  wsi_cols =
    names(
      COMPLETE_DATA
    )[
      which(
        names(
          COMPLETE_DATA
        ) %like% 'WSI_day_'
      )
      ]

  availhabitat[[i]] =
    COMPLETE_DATA[
      ,
      lapply(
        .SD,
        hospitable
      ),
      .SDcols = wsi_cols
      ]

  annmean_availhabitat[i] =
    mean(
      as.matrix(
        availhabitat[[i]]
      )
    )

  annsd_availhabitat[i] =
    sd(
      as.matrix(
        availhabitat[[i]]
      )
    )

  annmin_availhabitat[i] =
    min(
      as.matrix(
        availhabitat[[i]]
      )
    )

  # Mean WSI of available habitat
  MeanWSI_availhabitat[[i]] =
    as.matrix(
      COMPLETE_DATA[
        ,
        lapply(
          .SD,
          availmeanWSI
        ),
        .SDcols = wsi_cols
        ]
    )

  # Median WSI of available habitat
  MedianWSI_availhabitat[[i]] =
    as.matrix(
      COMPLETE_DATA[
        ,
        lapply(
          .SD,
          availmedianWSI
        ),
        .SDcols = wsi_cols
        ]
    )

  # Mean WSI of all habitat
  MeanWSI_allhabitat[[i]] =
    as.matrix(
      COMPLETE_DATA[
        ,
        lapply(
          .SD,
          ALLmeanWSI
        ),
        .SDcols = wsi_cols
        ]
    )

  # Median WSI of all habitat
  MedianWSI_allhabitat[[i]] =
    as.matrix(
      COMPLETE_DATA[
        ,
        lapply(
          .SD,
          ALLmedianWSI
        ),
        .SDcols = wsi_cols
        ]
    )

  # Median forage over year
  median_forage[i] =
    median(
      as.matrix(
        COMPLETE_DATA[
          ,
          names(
            COMPLETE_DATA
          )[
            names(
              COMPLETE_DATA
            ) %like%
              'NormForage_day_'
            ],
          with = FALSE
          ]
      )
    )

  # Minimum temperature over year
  Minimum_annualtemp[i] =
    min(
      COMPLETE_DATA[
        ,
        names(
          COMPLETE_DATA
        )[
          names(
            COMPLETE_DATA
          ) %like%
            'Temperature'
          ],
        with = FALSE
        ]
    )

  # Clean up environment
  rm(
    target_cols,
    COMPLETE_DATA
  )

  print(
    (
      i / (
        length(
          year_from_record
        ) - 1
      )
    ) * 100
  )
}

SummaryStats =
  data.table(
    Years =
      year_from_record[
        1:(
          length(
            year_from_record
          ) - 1
        )
        ],
    Mortality =
      unlist(
        lapply(
          daily_mort,
          sum
        )
      ),
    NormMortality =
      ZScore_func(
        unlist(
          lapply(
            daily_mort,
            sum
          )
        )
      ),
    AutumnMortality =
      unlist(
        lapply(
          daily_mort,
          function(
            x
          ) {
            sum(
              x[
                1:90
                ]
            )
          }
        )
      ),
    NormAutumnMortality =
      ZScore_func(
        unlist(
          lapply(
            daily_mort,
            function(
              x
            ) {
              sum(
                x[
                  1:90
                  ]
              )
            }
          )
        )
      ),
    WinterMortality =
      unlist(
        lapply(
          daily_mort,
          function(
            x
          ) {
            sum(
              x[
                91:180
                ]
            )
          }
        )
      ),
    NormWinterMortality =
      ZScore_func(
        unlist(
          lapply(
            daily_mort,
            function(
              x
            ) {
              sum(
                x[
                  91:180
                  ]
              )
            }
          )
        )
      ),
    SpringMortality =
      unlist(
        lapply(
          daily_mort,
          function(
            x
          ) {
            sum(
              x[
                181:num_days
                ]
            )
          }
        )
      ),
    NormSpringMortality = ZScore_func(
      unlist(
        lapply(
          daily_mort,
          function(
            x
          ) {
            sum(
              x[
                181:num_days
                ]
            )
          }
        )
      )
    ),
    Mean_AvailHab = annmean_availhabitat,
    StDev_AvailHab = annsd_availhabitat,
    NormMean_AvailHab =
      ZScore_func(
        annmean_availhabitat
      ),
    Min_AvailHab = annmin_availhabitat,
    NormAnnMinHab =
      ZScore_func(
        annmin_availhabitat
      ),
    MinTemp = Minimum_annualtemp,
    NormMinTemp =
      ZScore_func(
        Minimum_annualtemp
      ),
    MeanWSI_AvailHab =
      unlist(
        lapply(
          MeanWSI_availhabitat,
          mean
        )
      ),
    MedianWSI_AvailHab =
      unlist(
        lapply(
          MedianWSI_availhabitat,
          median
        )
      ),
    MeanWSI_AllHab =
      unlist(
        lapply(
          MeanWSI_allhabitat,
          mean
        )
      ),
    MedianWSI_AllHab =
      unlist(
        lapply(
          MedianWSI_allhabitat,
          median
        )
      ),
    NormAnnMeanWSI_AvailHab =
      ZScore_func(
        unlist(
          lapply(
            MeanWSI_availhabitat,
            mean
          )
        )
      ),
    NormAnnMedianWSI_AvailHab =
      ZScore_func(
        unlist(
          lapply(
            MedianWSI_availhabitat,
            median
          )
        )
      ),
    NormAnnMeanWSI_AllHab =
      ZScore_func(
        unlist(
          lapply(
            MeanWSI_allhabitat,
            mean
          )
        )
      ),
    NormAnnMedianWSI_AllHab =
      ZScore_func(
        unlist(
          lapply(
            MedianWSI_allhabitat,
            median
          )
        )
      )
  )

# Mean/Median mortality
SummaryStats[
  ,
  ann_mort :=
    (
      N_0 -
        Mortality
    ) /
      N_0
  ]

SummaryStats[
  ,
  mean(
    (
      N_0 -
        Mortality
    ) /
      N_0
  )
  ]

SummaryStats[
  ,
  range(
    (
      N_0 -
        Mortality
      ) /
      N_0
    )
  ]

SummaryStats[
  ,
  mean(
    (
      N_0 -
        Mortality
    ) /
      N_0
  ) +
    (
      (
        1.96 *
          sd(
            (
              N_0 -
                Mortality
            ) /
              N_0
          )
      )
    ) /
    nrow(
      SummaryStats
    )
  ]

SummaryStats[
  ,
  mean(
    (
      N_0 -
        Mortality
    ) /
      N_0
  ) -
    (
      (
        1.96 *
          sd(
            (
              N_0 -
                Mortality
            ) /
              N_0
          )
      )
    ) /
    nrow(
      SummaryStats
    )
  ]

SummaryStats[
  ,
  mean(
    (
      N_0 -
        AutumnMortality
    ) /
      N_0
  )
  ]

SummaryStats[
  ,
  median(
    (
      N_0 -
        AutumnMortality
    ) /
      N_0
  )
  ]

SummaryStats[
  ,
  range(
    (
      N_0 -
        AutumnMortality
    ) /
      N_0
  )
  ]

SummaryStats[
  ,
  mean(
    (
      (
        N_0 -
          AutumnMortality
      ) - WinterMortality
    ) /
      (
        N_0 -
          AutumnMortality
      )
  )
  ]

SummaryStats[
  ,
  median(
    (
      (
        N_0 -
          AutumnMortality
      ) -
        WinterMortality
    ) /
      (
        N_0 -
          AutumnMortality
      )
  )
  ]

SummaryStats[
  ,
  range(
    (
      (
        N_0 -
          AutumnMortality
      ) -
        WinterMortality
    ) /
      (
        N_0 -
          AutumnMortality
      )
  )
  ]

SummaryStats[
  ,
  mean(
    (
      (
        (
          N_0 -
            AutumnMortality
        ) -
          WinterMortality
      ) -
        SpringMortality
    ) /
      (
        (
          N_0 -
            AutumnMortality
        ) -
          WinterMortality
      )
  )
  ]

SummaryStats[
  ,
  median(
    (
      (
        (
          N_0 -
            AutumnMortality
        ) -
          WinterMortality
      ) -
        SpringMortality
    ) /
      (
        (
          N_0 -
            AutumnMortality
        ) -
          WinterMortality
      )
  )
  ]

SummaryStats[
  ,
  range(
    (
      (
        (
          N_0 -
            AutumnMortality
        ) -
          WinterMortality
      ) -
        SpringMortality
    ) /
      (
        (
          N_0 -
            AutumnMortality
        ) -
          WinterMortality
      )
  )
  ]

fwrite(
  SummaryStats,
  file =
    'FullRun_SummaryStats.txt'
)


#
# Visualize results graphically -------------------------------------------
SummaryStats = 
  fread(
    'FullRun_SummaryStats.txt'
  )

# NormMort_Plot =
#   ggplot(
#     data =
#       SummaryStats,
#     aes(
#       x = NormAnnMeanWSI_AllHab,
#       y = NormMortality
#     )
#   ) +
#   geom_point(
#     size = 2
#   ) +
#   scale_color_gradientn(
#     colours = 
#       rev(
#         viridis(
#           n = 5
#         )
#       )
#   ) +
#   xlab('Normalized WSI') +
#   ylab('Normalized Mortality') +
#   theme_minimal() +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18)
#   )
# 
# ggsave(
#   paste0(
#     "NormalizedMortality_v_NormalizedWSIAll.jpeg"
#   ),
#   NormMort_Plot,
#   width = 5,
#   height = 5
# )
# 
# NormAnnMinHab_Plot =
#   ggplot(
#     data = SummaryStats,
#     aes(
#       x = NormAnnMeanWSI_AllHab,
#       y = NormAnnMinHab,
#     )
#   ) +
#   geom_point(
#     size = 2
#   ) +
#   geom_smooth(
#     method = 'glm', 
#     color = 'black'
#   ) +
#   scale_color_gradientn(
#     colours = 
#       rev(
#         viridis(
#           n = 5
#         )
#       )
#   ) +
#   xlab('Nomalized WSI') +
#   ylab('Normalized Minimum Available Habitat') +
#   theme_minimal() +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18)
#   )
# 
# ggsave(
#   paste0(
#     'NormalizedAnnMinAvailHab_v_NormalizedWSIAll.jpeg'
#   ),
#   NormAnnMinHab_Plot,
#   width = 5,
#   height = 5
# )
# 
# Norm_Mort_MinHab_Plot =
#   ggplot(
#     data =
#       SummaryStats,
#     aes(
#       x = NormAnnMinHab,
#       y = NormMortality
#       )
#     ) +
#   geom_point(
#     size = 2
#     ) +
#   geom_smooth(
#     method = 'glm', 
#     color = 'black'
#     ) +
#   xlab(
#     'Normalized Minimum Available Habitat'
#   ) +
#   ylab(
#     'Normalized Mortality'
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18)
#   )
# 
# ggsave(
#   paste0(
#     'NormalizedMortality_v_NormalizedAnnMinAvailHab.jpeg'),
#   Norm_Mort_MinHab_Plot,
#   width = 5,
#   height = 5
# )
# 
# MeanWSIALL_Plot =
#   ggplot(
#     data = SummaryStats,
#     aes(
#       x = Years,
#       y = NormAnnMeanWSI_AllHab
#       )
#     ) +
#   geom_point(
#     size = 2
#     ) +
#   geom_smooth(
#     method = 'glm',
#     color = 'black'
#     ) +
#   xlab('Years') +
#   ylab('Normalized WSI of Landscape') +
#   theme_minimal() +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18)
#   )
# 
# ggsave(
#   paste0(
#     'NormMeanWSIAllHab.jpeg'
#   ),
#   MeanWSIALL_Plot,
#   width = 5,
#   height = 5
# )
# 
# NormMeanWSI_Plot =
#   ggplot(
#     data = SummaryStats,
#     aes(
#       x = Years,
#       y = NormAnnMeanWSI_AvailHab
#     )
#   ) +
#   geom_point(
#     size = 2
#   ) +
#   geom_smooth(
#     method = 'glm', 
#     color = 'black'
#   ) +
#   xlab('Years') +
#   ylab('Normalized WSI of Available Habitat') +
#   theme_minimal() +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18)
#   )
# 
# ggsave(
#   paste0(
#     'NormMeanWSIAvailHab.jpeg'
#   ),
#   NormMeanWSI_Plot,
#   width = 5,
#   height = 5
# )
# 
# Mort_HabWSI_Plot =
#   ggplot(
#     data = SummaryStats,
#     aes(
#       x = NormAnnMeanWSI_AvailHab,
#       y = NormMortality
#     )
#   ) +
#   geom_point(
#     size = 2
#   ) +
#   xlab('Normalized Mean WSI of Available Habitat') +
#   ylab('Normalized Mortality') +
#   theme_minimal()
# 
# ggsave(
#   paste0(
#     'Mortality_v_MeanWSIAvailHab.jpeg'
#   ),
#   Mort_HabWSI_Plot,
#   width = 5,
#   height = 5
# )
# 
# MeanHab_HabWSI_Plot =
#   ggplot(
#     data = SummaryStats,
#     aes(
#       x = NormAnnMeanWSI_AvailHab,
#       y = NormMean_AvailHab
#     )
#   ) +
#   geom_point(
#     size = 2
#   ) +
#   xlab('Normalized Mean WSI of Available Habitat') +
#   ylab('Normalized Mean Annual Available Habitat') +
#   theme_minimal()
# 
# ggsave(
#   paste0(
#     'ProporMeanAnnAvailHab_v_MeanWSIAvailHab.jpeg'
#   ),
#   MeanHab_HabWSI_Plot,
#   width = 5,
#   height = 5
# )

SDAvailHab_MeanDistAmong_Plot =
  ggplot(
    data = SummaryStats,
    aes(
      x = StDev_AvailHab,
      y = MeanDist_Among
    )
  ) +
  geom_point(
    size = 2
  ) +
  geom_smooth(
    method = 'lm',
    col = 'black'
  ) +
  xlab('Standard Deviation of Available Habitat') +
  ylab('Mean Distance Among Center-of-Mass Locations') +
  theme_minimal()

ggsave(
  paste0(
    'SDAvailHab_v_MeanDistAmongCoM.jpeg'
  ),
  SDAvailHab_MeanDistAmong_Plot,
  width = 5,
  height = 5
)

# Gather proportion of available habitat through time
annualavailhabitat_dt =
  data.table(
    cbind(
      c(1:num_days),
      do.call(
        cbind,
        lapply(
          availhabitat,
          function(
            x
          ) {
            as.vector(
              as.matrix(
                unlist(
                  x
                )
              )
            )
          }
        )
      )
    )
  )

setnames(
  annualavailhabitat_dt,
  names(
    annualavailhabitat_dt
  ),
  c(
    'MigDay',
    paste0(
      'AvailHab_',
      year_from_record[
        1:(
          length(
            year_from_record
          ) - 1
        )
        ]
    )
  )
)

annualavailhabitat_melt =
  melt(
    annualavailhabitat_dt,
    id.vars = 'MigDay',
    variable.name = 'Year',
    value = 'ProporAvailHab'
  )

annualavailhabitat_melt[
  ,
  Year :=
    gsub(
      'AvailHab_',
      '',
      Year
    )
  ]

setkey(
  annualavailhabitat_melt,
  MigDay,
  Year
)

daily_mortality =
  data.table(
    do.call(
      cbind,
      daily_mort
    )
  )

setnames(
  daily_mortality,
  names(
    daily_mortality
  ),
  paste0(
    'Mortality_',
    year_from_record[
      1:(
        length(
          year_from_record
        ) - 1
      )
      ]
  )
)

MigDay =
  1:length(
    days
  )

daily_mortality =
  cbind(
    MigDay,
    daily_mortality
  )

dailymortality_melt =
  melt(
    daily_mortality,
    id.vars = 'MigDay',
    variable.name = 'Year',
    value = 'DailyMortality'
  )

dailymortality_melt[
  ,
  Year :=
    gsub(
      'Mortality_',
      '',
      Year
    )
  ]

setkey(
  dailymortality_melt,
  MigDay,
  Year
)

annualavailhabitat_dailymortality_melt =
  dailymortality_melt[
    annualavailhabitat_melt
    ]

# # Plot of avarage daily mortality
# dailymortality_melt[
#   ,
#   mean_daily_mort :=
#     mean(
#       DailyMortality
#     ),
#   by = MigDay
# ]
# 
# mean_daily_mort_plot =
#   ggplot(
#     data = dailymortality_melt,
#     aes(
#       x = MigDay,
#       y = mean_daily_mort
#     )
#   ) +
#   geom_line() +
#   xlab('Day of Non-breeding Period') +
#   ylab('Number of birds dying') +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 14)
#   )
# 
# ggsave(
#   'mean_daily_mortality_plot.jpeg',
#   mean_daily_mort_plot,
#   width = 5,
#   height = 5
# )

# Evaluate mild versus severe years
SummaryStats[
  NormAnnMeanWSI_AllHab <
    quantile(
      NormAnnMeanWSI_AllHab,
      probs =
        seq(
          0,
          1,
          0.25
        )
    )[2],
  severity := 'mild'
]

SummaryStats[
  NormAnnMeanWSI_AllHab >
    quantile(
      NormAnnMeanWSI_AllHab,
      probs =
        seq(
          0,
          1,
          0.25
        )
    )[4],
  severity := 'severe'
]

severe_years = 
  SummaryStats[
    severity == 'severe',
    Years
  ]

mild_years = 
  SummaryStats[
    severity == 'mild',
    Years
  ]

annualavailhabitat_melt[
  Year %in% severe_years,
  sev_cat := 'Severe'
]

annualavailhabitat_melt[
  Year %in% mild_years,
  sev_cat := 'Mild'
]

mild_regress_est = 
  annualavailhabitat_melt[
    sev_cat == 'Mild',
    lm(
      ProporAvailHab ~ 
        MigDay + I(
          MigDay ^ 2
        )
    )
  ]

severe_regress_est = 
  annualavailhabitat_melt[
    sev_cat == 'Severe',
    lm(
      ProporAvailHab ~ 
        MigDay + I(
          MigDay ^ 2
        )
    )
  ]

mild_severe_comparison_dt = 
  data.table(
    Mig_day = 1:273,
    Severe = 
      severe_regress_est$coefficients[1] +
      (
        severe_regress_est$coefficients[2] *
          c(1:273)
      ) +
      (
        severe_regress_est$coefficients[3] *
          I(
            c(1:273) ^ 2
          )
      ),
    Mild = 
      mild_regress_est$coefficients[1] +
      (
        mild_regress_est$coefficients[2] *
          c(1:273)
      ) +
      (
        mild_regress_est$coefficients[3] *
          I(
            c(1:273) ^ 2
          )
      )
  )

mild_severe_comparison_melt = 
  data.table::melt(
    mild_severe_comparison_dt,
    id.vars = 'Mig_day',
    variable.name = 'Severity',
    value.name = 'Proportion_availhab'
  )

weather_severity_comparison =
  ggplot(
    data = mild_severe_comparison_melt,
    aes(
      x = Mig_day,
      y = Proportion_availhab,
      linetype = Severity,
      pch = Severity
    )
  ) +
  geom_smooth(
    method = 'lm',
    formula = y ~ poly(x, 2),
    col = 'black',
    fullrange = TRUE
  ) +
  ylim(0.1, 0.3) +
  xlim(75, 200) +
  xlab('Day of Non-breeding Period') +
  ylab('Proportion of Available Habitat') +
  labs(
    linetype = 'Weather Severity\nCategory',
    shape = 'Weather Severity\nCategory'
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(
  'wsi_comparison_plot.jpeg',
  weather_severity_comparison,
  width = 7,
  height = 5
)

# annualavailhabitat_dailymortality_melt[
#   Year %in% severe_years,
#   sev_cat :=
#     ifelse(
#       Year %in%
#         c(
#           1957,
#           1962,
#           1989
#         ),
#       'Freeze',
#       'Mild'
#     )
#   ]
# 
# annualavailhabitat_dailymortality_melt[
#   ,
#   DailyMortality:=
#     ZO_std(
#       DailyMortality
#     )
#   ]
# 
# annualavailhabitat_dailymortality_melt_dt =
#   melt(
#     annualavailhabitat_dailymortality_melt,
#     id.vars =
#       c(
#         'MigDay',
#         'Year',
#         'Extreme'
#       ),
#     measure.vars =
#       c(
#         'DailyMortality',
#         'ProporAvailHab'
#       ),
#     variable.name = 'Variable',
#     value.name = 'Normalized_Value'
#   )
# 
# Annual_AvailHab_Plot =
#   ggplot(
#     data = annualavailhabitat_melt,
#     aes(
#       x = MigDay,
#       y = ProporAvailHab,
#       group = Year
#     )
#   ) +
#   geom_line(
#     alpha = 0.1
#   ) +
#   scale_y_continuous(
#     breaks =
#       c(
#         0,
#         0.2,
#         0.4,
#         0.6,
#         0.8,
#         1
#       )
#   ) +
#   xlab('Day of Non-breeding Period') +
#   ylab('Normalized Mean Annual Available Habitat') +
#   theme_minimal()
# 
# ggsave(
#   paste0(
#     'ProporAnnAvailHab_OverTime.jpeg'
#   ),
#   Annual_AvailHab_Plot,
#   width = 7,
#   height = 7
# )
# 
# Annual_AvailHab_Mort_Plot =
#   ggplot(
#     data = annualavailhabitat_dailymortality_melt_dt,
#     aes(
#       x = MigDay,
#       y = Normalized_Value,
#       linetype = Variable
#     )
#   ) +
#   geom_line(
#     lwd = 1
#   ) +
#   facet_wrap( 
#     ~ Year
#   ) +
#   xlab('Day of Non-breeding Period') +
#   ylab('Normalized Mean Value (Available Habitat or (1 - Daily Mortality))') +
#   theme_minimal()
# 
# ggsave(
#   paste0(
#     'ProporAnnAvailHab_Mortality_OverTime_extremes.jpeg'
#   ),
#   Annual_AvailHab_Mort_Plot,
#   width = 10,
#   height = 5
# )
# 
# annualavailhabitat_dailymortality_melt_dt[
#   ,
#   Mean_Normalized_Extreme_Value :=
#     mean(
#       Normalized_Value
#     ),
#   by =
#     .(
#       MigDay,
#       Extreme
#     )
#   ]
# 
# Annual_AvailHab_Plot_smoothed =
#   ggplot(
#     data =
#       annualavailhabitat_dailymortality_melt_dt[
#         Variable %like% 'ProporAvailHab'
#         ],
#     aes(
#       x = MigDay,
#       y = Mean_Normalized_Extreme_Value,
#       linetype = Extreme
#     )
#   ) +
#   geom_smooth(
#     lwd = 1,
#     col = 'black'
#   ) +
#   scale_y_continuous(
#     breaks =
#       c(
#         0,
#         0.2,
#         0.4,
#         0.6,
#         0.8,
#         1
#       )
#   ) +
#   xlab('Day of Non-breeding Period') +
#   ylab('Normalized Mean Annual Available Habitat') +
#   theme_minimal()
# 
# ggsave(
#   paste0(
#     'ProporAnnAvailHab_OverTime_smoothed.jpeg'
#   ),
#   Annual_AvailHab_Plot_smoothed,
#   width = 5,
#   height = 5
# )
# 
# # With state lines
# na_countries_sts = 
#   c(
#     'USA',
#     'CAN'
#   )
# 
# northamerica_sts =
#   do.call(
#     'bind',
#     lapply(
#       na_countries_sts,
#       function(
#         x
#       ) {
#         getData(
#           'GADM',
#           country = x,
#           level = 1
#         )
#       }
#     )
#   )
# 
# # Remove Alaska and Hawaii
# CONUS_NA_sts =
#   northamerica_sts[
#     !(
#       toupper(
#         northamerica_sts$NAME_1
#       ) %like%
#         'ALASKA'|
#         toupper(
#           northamerica_sts$NAME_1
#         ) %like%
#         'HAWAII'
#     ),
#     ]
# 
# writeOGR(
#   obj = CONUS_NA_sts,
#   dsn = getwd(),
#   layer = 'states',
#   driver = 'ESRI Shapefile',
#   overwrite_layer = TRUE
# )
# 
CONUS_NA_sts =
  readOGR(
    'states',
    dsn = getwd()
  )

# Obtain spatial data for plotting
CoM_DailyDist_Between = CoM_DailyDist_Among = CoM_DailyDir =
  CenterofMass_coords = vector('list')

SumCoM_AdjDist = MeanCoM_AdjDist_Between = MeanCoM_AdjDist_Among =
  CoM_NS_Dist = Min_CoM_Lat_Day = Min_CoM_Lat = vector()

for(
  i in 1:(
    length(
      year_from_record
      ) - 1
    )
) {
  COMPLETE_DATA =
    fread(
      paste0(
        year_from_record[i],
        '_Prediction.txt'
      )
    )
  
  target_cols =
    names(
      COMPLETE_DATA
    )[
      which(
        names(
          COMPLETE_DATA
        ) %like% 'N_Abund_'
      )
      ][
        2:(
          num_days + 1
        )
        ]
  
  most_populace_node =
    as.vector(
      as.matrix(
        COMPLETE_DATA[
          ,
          lapply(
            .SD,
            node_populace
          ),
          .SDcols = target_cols
          ]
      )
    )
  
  COMPLETE_DATA[
    most_populace_node,
    within_populace := 1
    ]
  
  COMPLETE_DATA[
    ,
    within_populace :=
      NAFunc(
        within_populace
      )
    ]
  
  Mig_route = 
    COMPLETE_DATA[
      ,
      .(
        Longitude,
        Latitude,
        within_populace
      )
      ]
  
  coordinates(
    Mig_route
  ) = 
    c(
      'Longitude',
      'Latitude'
    )
  
  gridded(
    Mig_route
  ) = TRUE
  
  Mig_route_grd = 
    as(
      Mig_route, 
      'SpatialGridDataFrame'
    )
  
  crs(
    Mig_route_grd
  ) =
    '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
  +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
  
  Mig_route_rstr = 
    raster(
      Mig_route_grd
    )
  
  Mig_route_proj = 
    projectRaster(
      Mig_route_rstr,
      crs = geo_prj,
      res = 0.1
    )
  
  CenterofMass = 
    t(
      as.matrix(
        COMPLETE_DATA[
          ,
          lapply(
            .SD,
            Pop_CoM,
            x = Longitude,
            y = Latitude,
            z = NULL
          ),
          .SDcols = target_cols
          ][
            c(
              1,
              3
            )
            ]
      )
    )
  
  CenterofMass_pts = 
    SpatialPoints(
      CenterofMass
      )

  crs(
    CenterofMass_pts
    ) =
    '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
  +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
  
  CenterofMass_proj = 
    spTransform(
      CenterofMass_pts,
      CRS(
        geo_prj
      )
    )
  
  CenterofMass_coords[[i]] = 
    data.table(
      coordinates(
        CenterofMass_proj
      )
    )
  
  CenterofMass_coords[[i]][
    ,
    ColorIterate :=
      c(
        1:nrow(
          CenterofMass_coords[[i]]
        )
      )
    ]
  
  CoM_DailyDist_Between[[i]] =
    sapply(
      2:nrow(
        CenterofMass_coords[[i]]
      ),
      function(x) {
        distm(
          CenterofMass_coords[[i]][
            x - 1,
            .(
              coords.x1,
              coords.x2
            )
            ],
          CenterofMass_coords[[i]][
            x,
            .(
              coords.x1,
              coords.x2
            )
            ]
        ) /
          1000
      }
    )
  
  CoM_DailyDist_Among[[i]] =
    distm(
      CenterofMass_coords[[i]][
        ,
        .(
          coords.x1,
          coords.x2
        )
        ]
    ) / 1000
  
  CoM_DailyDist_Among[[i]] =
    ifelse(
      CoM_DailyDist_Among[[i]] == 0,
      NA,
      CoM_DailyDist_Among[[i]]
    )
  
  CoM_DailyDir[[i]] =
    CenterofMass_coords[[i]][
      ,
      coords.x2 -
        data.table::shift(
          coords.x2,
          n = 1L,
          fill = 0,
          type = 'lag'
        )
      ][
        -1
        ]
  
  Min_CoM_Lat_Day[i] = 
    CenterofMass_coords[[i]][
      ,
      which.min(
        coords.x2
      )
    ]
  
  Min_CoM_Lat[i] = 
    CenterofMass_coords[[i]][
      which.min(
        coords.x2
      ),
      coords.x2
    ]
  
  CoM_NS_Dist[i] =
    distm(
      CenterofMass_coords[[i]][
        which.min(
          coords.x2
        ),
        .(
          coords.x1,
          coords.x2
        )
      ],
      CenterofMass_coords[[i]][
        273,
        .(
          coords.x1,
          coords.x2
        )
      ]
    ) / 1000
      
  SumCoM_AdjDist[i] = 
    sum(
      CoM_DailyDist_Between[[i]]
    )
  
  MeanCoM_AdjDist_Between[i] = 
    mean(
      CoM_DailyDist_Between[[i]],
      na.rm = TRUE
    )
  
  MeanCoM_AdjDist_Among[i] = 
    mean(
      CoM_DailyDist_Among[[i]],
      na.rm = TRUE
    )
  
  col_pal = 
    colorRampPalette(
      c(
        'darkred',
        'deepskyblue4'
      )
    )
  
  CenterofMass_coords[[i]][
    ,
    ColorGradient :=
      col_pal(
        10
      )[
        as.numeric(
          cut(
            ColorIterate,
            breaks = 10
          )
        )
        ]
    ]
  
  tiff(
    file = 
      paste0(
        year_from_record[i],
        '_MigrationRoute_CoM_STPRVBRDR.tiff'
      ),
    bg = 'transparent',
    width = 1500,
    height = 900
  )
  
  par(
    mar =
      c(
        0,
        0,
        0,
        0
      ),
    oma =
      c(
        0,
        0,
        0,
        0
      )
  )
  
  image(
    Mig_route_proj,
    col =
      c(
        'white',
        'gray75'
      ),
    bty = 'n',
    xlab = NA,
    ylab = NA,
    xaxt = 'n',
    yaxt = 'n'
  )
  
  plot(
    CONUS_NA_sts,
    lty = 'dotted',
    border = 'gray',
    add = TRUE
  )
  
  points(
    CenterofMass_coords[[i]][
      ,
      .(
        coords.x1,
        coords.x2
      )
      ],
    pch = 16,
    cex = 0.5,
    col = 
      CenterofMass_coords[[i]][
        ,
        ColorGradient
        ]
  )
  
  dev.off()
  
  rm(COMPLETE_DATA)
  
  print(
    (
      i / (
        length(
          year_from_record
        ) - 1
      )
    ) * 100
  )
}

SummaryStats[
  ,
  `:=`(
    TotalDist_Traveled = SumCoM_AdjDist,
    MeanDist_Between = MeanCoM_AdjDist_Between,
    MeanDist_Among = MeanCoM_AdjDist_Among,
    NorthSouthDist_Traveled = CoM_NS_Dist
  )
  ]

fwrite(
  SummaryStats,
  file = 'FullRun_SummaryStats.txt'
)

SummaryStats =
  fread(
    file = 'FullRun_SummaryStats.txt'
  )

DistanceBetweenMetric =
  data.table(
    do.call(
      cbind,
      CoM_DailyDist_Between
    )
  )

setnames(
  DistanceBetweenMetric,
  names(
    DistanceBetweenMetric
  ),
  c(
    paste0(
      'DistanceBetween_',
      year_from_record[
        -which.max(
          year_from_record
        )
        ]
    )
  )
)

DistanceAmongMetric =
  matrix(
    rep(
      MeanCoM_AdjDist_Among,
      each = num_days - 1
    ),
    nrow = num_days - 1
  )

DistanceAmongMetric = 
  data.table(
    DistanceAmongMetric
  )

setnames(
  DistanceAmongMetric,
  names(
    DistanceAmongMetric
  ),
  c(
    paste0(
      'DistanceAmong_',
      year_from_record[
        -which.max(
          year_from_record
        )
        ]
    )
  )
)

DirectionMetric =
  data.table(
    do.call(
      cbind,
      CoM_DailyDir
    )
  )

setnames(
  DirectionMetric,
  names(
    DirectionMetric
  ),
  c(
    paste0(
      'Direction_',
      year_from_record[
        -which.max(
          year_from_record
        )
        ]
    )
  )
)

MovementMetrics = 
  cbind(
    DayofPeriod=days[
      -which.max(
        days
      )
      ],
    DistanceBetweenMetric,
    DistanceAmongMetric,
    DirectionMetric
  )

fwrite(
  MovementMetrics,
  file = 'MovementMetrics.txt'
)

MovementMetrics = 
  fread(
    'MovementMetrics.txt'
  )

for(
  i in 1:(
    length(
      year_from_record
      ) - 1
    )
) {
  CenterofMass_coords[[i]][,`:=`(
    Year = 
      rep(
        year_from_record[i],
        nrow(
          CenterofMass_coords[[i]]
          )
        ),
    MeanWSIAvail = 
      as.vector(
        MeanWSI_availhabitat[[i]]
        ),
    MedianWSIAvail = 
      as.vector(
        MedianWSI_availhabitat[[i]]
        ),
    MeanWSIAll = 
      as.vector(
        MeanWSI_allhabitat[[i]]
        ),
    MedianWSIAll = 
      as.vector(
        MedianWSI_allhabitat[[i]]
        ),
    Mean_AvailHab = 
      as.vector(
        as.matrix(
          availhabitat[[i]]
                  )
        ),
    Norm_AvailHab = 
      as.vector(
        ZScore_func(
          unlist(
            availhabitat[[i]]
            )
          )
        )
    )
    ]
}

CenterofMass_coords = 
  data.table(
    do.call(
      rbind,
      CenterofMass_coords
      )
    )

fwrite(
  CenterofMass_coords,
  file = 'CenterofMass_coords.txt'
  )

CenterofMass_coords = 
  fread(
    'CenterofMass_coords.txt'
  )

CenterofMass_coords[
  ,
  `:=`(
    meanDailyLat =
      mean( 
        coords.x2
      ),
    meanDailyWSIAvail =
      mean(
        MeanWSIAvail
      ),
    medianDailyWSIAvail =
      mean(
        MedianWSIAvail
      ),
    meanDailyWSIAll =
      mean(
        MeanWSIAll
      ),
    medianDailyWSIAll =
      mean(
        MedianWSIAll
      ),
    meanDailyAvailHab =
      mean(
        Mean_AvailHab
      )
  ),
  by = ColorIterate
  ]

CenterofMass_coords[
  ,
  `:=`(
    ColorIterate = NULL,
    ColorGradient = NULL
  )
  ]

CenterofMass_coords[
  ,
  `:=`(
    AnnualMeanWSIAvail = 
      mean(
        MeanWSIAvail
      ),
    AnnualMedianWSIAvail =
      mean(
        MedianWSIAvail
      ),
    AnnualMeanWSIAll =
      mean( 
        MeanWSIAll
      ),
    AnnualMedianWSIAll =
      mean(
        MedianWSIAll
      )
  ),
  by = Year
  ]

col_pal = 
  colorRampPalette(
    c(
      'red',
      'blue'
      )
    )

CenterofMass_coords[
  ,
  `:=`(
    ColorGradient =
      col_pal(
        57
      )[
        as.numeric(
          cut(
            MeanWSIAvail,
            breaks = 10
          )
        )
        ],
    NonbreedingDay = 
      rep(
        1:num_days,
        (
          length(
            year_from_record
          ) - 1
        )
      )
  )
  ]

scalemin = 
  CenterofMass_coords[
    ,
    trunc(
      min(
        AnnualMedianWSIAvail
      ) * 100
    ) / 100
    ]

scalemax =
  CenterofMass_coords[
    ,
    trunc(
      max(
        AnnualMedianWSIAvail
      ) * 100
    ) / 100
    ]

LatAbundTime_plot =
  ggplot(
    data =
      CenterofMass_coords,
    aes(
      x = NonbreedingDay,
      y = coords.x2,
      group = Year,
      col = AnnualMedianWSIAvail
    )
  ) +
  geom_line(
    stat = 'smooth',
    method = 'loess',
    se = FALSE,
    size = 1.2,
    alpha = 0.4
  ) +
  scale_color_continuous(
    name='Annual Median \nWeather Severity',
    breaks =
      c(
        seq(
          floor(scalemin),
          ceiling(scalemax),
          round(
            (
              ceiling(scalemax) - floor(scalemin)
            ) / 6,
            1
          )
        )
      ),
    low = 'gray90',
    high = 'black'
  ) +
  xlab('Day of Nonbreeding Period') +
  ylab('Latitude') +
  theme_minimal()

ggsave(
  'LatAbundTime.jpeg',
  LatAbundTime_plot,
  width = 5,
  height = 5
)

LatAbundWSI =
  CenterofMass_coords[
    ,
    .(
      min(
        coords.x2
      ),
      MedianWSIAvail[
        which.min(
          coords.x2
        )
        ]
    ),
    by = Year
    ]

setnames(
  LatAbundWSI,
  names(
    LatAbundWSI
  ),
  c(
    'Year',
    'Latitude',
    'WSI'
  )
)

LatAbundWSI_Plot =
  ggplot(
    data = LatAbundWSI,
    aes(
      x = WSI,
      y = Latitude
    )
  ) +
  geom_point() +
  geom_smooth(
    method = 'glm',
    color = 'black'
  ) +
  xlab('Weather Severity') +
  ylab('Latitude') +
  theme_minimal()

ggsave(
  'LatAbundWSI.jpeg',
  LatAbundWSI_Plot,
  width = 5,
  height = 5
)

#Regression of Most populace minimum latitude with WSI
MPML_WSI =
  CenterofMass_coords[
    ,
    .(
      min(
        coords.x2
      ),
      min(
        MeanWSIAvail
      ),
      mean(
        MeanWSIAvail
      ),
      min(
        MedianWSIAvail
      ),
      mean(
        MedianWSIAvail
      ),
      min(
        MeanWSIAll
      ),
      mean(
        MeanWSIAll
      ),
      min(
        MedianWSIAll
      ),
      mean(
        MedianWSIAll
      ),
      min(
        Mean_AvailHab
      ),
      mean(
        Mean_AvailHab
      )
    ),
    by = Year
    ]

setnames(
  MPML_WSI,
  names(
    MPML_WSI
  ),
  c(
    'Year',
    'Latitude',
    'Min_MeanWSIAvail',
    'Mean_MeanWSIAvail',
    'Min_MedianWSIAvail',
    'Mean_MedianWSIAvail',
    'Min_MeanWSIAll',
    'Mean_MeanWSIAll',
    'Min_MedianWSIAll',
    'Mean_MedianWSIAll',
    'Min_Mean_AvailHab',
    'AnnualMean_DailyMeanAvailableHabitat'
  )
)

fwrite(
  MPML_WSI,
  'most_populace_minimum_latitude.csv'
)

# # Plot body condition distribution over time
# Sys.time()
# start_timer = Sys.time()
# 
# n_cores = detectCores() - 1
# start_cluster = makeCluster(n_cores)
# registerDoParallel(start_cluster)
# 
# foreach(
#   j = 1:(
#     length(
#       year_from_record
#     ) - 1
#   ),
#   .packages = 'data.table',
#   .verbose = TRUE
# ) %dopar% {
#   BodyCondition_Table = 
#     fread(
#       file =
#         paste0(
#           year_from_record[j],
#           "_BCTable.txt"
#         )
#     )
#   
#   bc_dist_cols = 
#     names(
#       BodyCondition_Table
#     )[
#       names(
#         BodyCondition_Table
#         ) %like%
#         'BodyMass_dist_discrete'
#     ]
#   
#   BC_tbl_melt = 
#     data.table::melt(
#       BodyCondition_Table,
#       id.vars = 'BC_bins',
#       measure.vars = bc_dist_cols,
#       variable.name = 'Dist_day',
#       value.name = 'Proportion'
#     )
#   
#   BC_tbl_melt[
#     ,
#     Dist_day :=
#       gsub(
#         'BodyMass_dist_discrete_day',
#         '',
#         Dist_day
#       )
#   ]
#   
#   BC_plot_time = 
#     ggplot(
#       BC_tbl_melt[
#         Dist_day %in% c(0, 100, 200, 273)
#       ],
#       aes(
#         x = BC_bins,
#         y = Proportion,
#         col = Dist_day,
#         group = Dist_day
#       )
#     ) + 
#     geom_line() +
#     geom_area(
#       aes(
#         fill = Dist_day
#       ),
#       alpha = 0.2,
#       position = 'identity'
#     ) +
#     xlab('Body Condition Class') +
#     ylab('Abundance Proportion') +
#     labs(
#       color = 'Non-breeding\nPeriod Day'
#     ) +
#     guides(
#       fill = 'none'
#     ) +
#     theme_minimal() +
#     theme(
#       axis.title = element_text(size = 16),
#       axis.text = element_text(size = 14),
#       legend.text = element_text(size = 14),
#       legend.title = element_text(size = 16)
#     )
#   
#   ggsave(
#     'BodyCondition_overtime.png',
#     BC_plot_time,
#     width = 5,
#     height = 7
#   )
#   
#   print(j)
# }

stopCluster(start_cluster)

#
# Evaluate and visualize results ------------------------------------------
MPML_WSI =
  fread(
    'most_populace_minimum_latitude.csv'
  )

MPML_Year_plot =
  ggplot(
    data = MPML_WSI,
         aes(
           x = Year,
           y = Latitude
           )
    ) +
  geom_point() +
  geom_smooth(
    method = 'glm'
    ) +
  xlab('Year') +
  ylab('Latitude') +
  theme_minimal()

ggsave(
  'MPML_Year.jpeg',
  MPML_Year_plot,
  width = 5,
  height = 5
)

MPML_WSI_plot =
  ggplot(
    data =
      MPML_WSI,
    aes(
      x = Mean_MedianWSIAll,
      y = Latitude
    )
  ) +
  geom_point() +
  geom_smooth(
    method = 'glm'
  ) +
  xlab('WSI') +
  ylab('Latitude') +
  theme_minimal()

ggsave(
  'MPML_WSI.jpeg',
  MPML_WSI_plot,
  width = 5,
  height = 5
)

Year_WSI_plot =
  ggplot(
    data = MPML_WSI,
    aes(
      x = Year,
      y = Mean_MeanWSIAll
    )
  ) +
  geom_point() +
  geom_smooth(
    method = 'glm'
  ) +
  # scale_y_continuous(
  #   breaks =
  #     c(
  #       -9,
  #       -8.5,
  #       -8,
  #       -7.5,
  #       -7,
  #       -6.5,
  #       -6
  #     )
  # ) +
  xlab('Year') +
  ylab('WSI') +
  theme_minimal()

ggsave(
  'Year_WSI.jpeg',
  Year_WSI_plot,
  width = 5,
  height = 5
)

SummaryStats = 
  fread(
    'FullRun_SummaryStats.txt'
  )

Dist_Among =
  ggplot(
    data = SummaryStats,
    aes(
      x = Years,
      y = MeanDist_Among
    )
  ) +
  geom_point(
    size = 2
  ) + 
  geom_line( 
    lty = 'dotted'
  ) +
  xlab('Years') +
  ylab('Mean Distance (km) Among all Center-of-mass Points') +
  theme_minimal()

ggsave(
  paste0(
    'Among_Dist.jpeg'
  ),
  Dist_Among,
  width = 5,
  height = 5
)

Dist_Btwn =
  ggplot(
    data = SummaryStats,
    aes(
      x = Years,
      y = MeanDist_Between
    )
  ) +
  geom_point(
    size = 2
  ) + 
  geom_line(
    lty = 'dotted'
  ) +
  xlab('Years') +
  ylab('Mean Distance (km) Between Adjacent Center-of-mass Points') +
  theme_minimal()

ggsave(
  paste0(
    'Between_Dist.jpeg'
  ),
  Dist_Btwn,
  width = 5,
  height = 5
)

Dist_NS =
  ggplot(
    data = SummaryStats,
    aes(
      x = Years,
      y = NorthSouthDist_Traveled
    )
  ) +
  geom_point(
    size = 2
  ) + 
  geom_line(
    lty = 'dotted'
  ) +
  xlab('Years') +
  ylab('Mean Distance (km) Between Northern and Southern Extremes') +
  theme_minimal()

ggsave(
  paste0(
    'NS_Dist.jpeg'
  ),
  Dist_NS,
  width = 5,
  height = 5
)

# Abundance-weighted population center-of-mass across years
CenterofMass_data = data.table()

for(
  i in 1:(
    length(
      year_from_record
      ) - 1
    )
) {
  COMPLETE_DATA = 
    fread(
      paste0(
        year_from_record[i],
        '_Prediction.txt'
      )
    )
  
  target_cols =
    names(
      COMPLETE_DATA
    )[
      which(
        names(
          COMPLETE_DATA
        ) %like% 
          'N_Abund_'
      )
      ][
        2:(
          num_days + 1
        )
        ]
  
  CenterofMass = 
    t(
      as.matrix(
        COMPLETE_DATA[
          ,
          lapply(
            .SD,
            Pop_CoM,
            x = Longitude,
            y = Latitude,
            z = NULL
          ),
          .SDcols =
            target_cols
          ][
            c(
              1,
              3
            )
            ]
      )
    )
  
  CenterofMass_pts = 
    SpatialPoints(
      CenterofMass
    )
  
  crs(
    CenterofMass_pts
  ) =
    '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
  +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
  
  CenterofMass_proj = 
    spTransform(
      CenterofMass_pts,
      CRS(
        geo_prj
      )
    )
  
  CenterofMass_data_init = 
    data.table(
      Year = 
        rep(
          year_from_record[i],
          nrow(
            coordinates(
              CenterofMass_proj
            )
          )
        ),
      coordinates(
        CenterofMass_proj
      )
    )
  
  setnames(
    CenterofMass_data_init,
    names(
      CenterofMass_data_init
    ),
    c(
      'Year',
      'Longitude',
      'Latitude'
    )
  )
  
  CenterofMass_data = 
    rbind(
      CenterofMass_data, 
      CenterofMass_data_init
    )
  
  rm(COMPLETE_DATA)
  
  print( 
    (
      i / (
        length(
          year_from_record
        ) - 1
      )
    )  * 100
  )
}

fwrite(
  CenterofMass_data,
  file = 'CenterofMass_data.txt'
  )

CenterofMass_data = 
  fread(
    'CenterofMass_data.txt'
  )

CenterofMass_data[
  ,
  NonBreedingDay :=
    rep(
      days,
      (
        length(
          year_from_record
        ) - 1 
      )
    )
  ]

CenterofMass_data[
  , 
  `:=`(
    mean_Longitude = 
      mean(
        Longitude
      ),
    mean_Latitude =
      mean(
        Latitude
      )
  ),
  by = NonBreedingDay
  ]

# Obtain most populace nodes
COMPLETE_DATA = 
  fread(
    paste0(
      year_from_record[i],
      '_Prediction.txt'
    )
  )

# Migration route (based on most populous node per day)
target_cols =
  names(
    COMPLETE_DATA
  )[
    which(
      names(
        COMPLETE_DATA
      ) %like% 
        'N_Abund_'
    )
    ][
      2:(
        num_days + 1
      )
      ]

most_populace_node =
  as.vector(
    as.matrix(
      COMPLETE_DATA[
        ,
        lapply(
          .SD,
          node_populace
        ),
        .SDcols = target_cols
        ]
    )
  )

COMPLETE_DATA[
  most_populace_node,
  within_populace := 1
  ]

Mig_route = 
  COMPLETE_DATA[
    ,
    .(
      Longitude,
      Latitude,
      within_populace
      )
    ]

coordinates(
  Mig_route
  ) = 
  c(
    'Longitude',
    'Latitude'
    )

gridded(
  Mig_route
  ) = TRUE

Mig_route_grd = 
  as(
    Mig_route, 
    'SpatialGridDataFrame'
    )

crs(
  Mig_route_grd
  ) =
  '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

Mig_route_rstr = 
  raster(
    Mig_route_grd
  )

Mig_route_proj = 
  projectRaster(
    Mig_route_rstr,
    crs = geo_prj,
    res = 0.1
  )

image(
  Mig_route_proj,
  col =
    c(
      'white',
      'gray75'
    ),
  bty = 'n',
  xlab = NA,
  ylab = NA,
  xaxt = 'n',
  yaxt = 'n'
)

rm(COMPLETE_DATA)

for(
  i in 2:(
    length(
      year_from_record
    ) - 1
  )
) {
  COMPLETE_DATA = 
    fread(
      paste0(
        year_from_record[i],
        '_Prediction.txt'
      )
    )
  
  # Migration route (based on most populous node per day)
  target_cols =
    names(
      COMPLETE_DATA
    )[
      which(
        names(
          COMPLETE_DATA
        ) %like% 'N_Abund_'
      )
      ][
        2:(
          num_days + 1
        )
        ]
  
  most_populace_node =
    as.vector(
      as.matrix(
        COMPLETE_DATA[
          ,
          lapply(
            .SD,
            node_populace
          ),
          .SDcols = target_cols
          ]
      )
    )
  
  COMPLETE_DATA[
    most_populace_node,
    within_populace := 1
    ]

  Mig_route = 
    COMPLETE_DATA[
      ,
      .(
        Longitude,
        Latitude,
        within_populace
        )
      ]
  
  coordinates(
    Mig_route
  ) = 
    c(
      'Longitude',
      'Latitude'
    )
  
  gridded(
    Mig_route
  ) = TRUE
  
  Mig_route_grd = 
    as(
      Mig_route, 
      'SpatialGridDataFrame'
    )
  
  crs(
    Mig_route_grd
  ) =
    '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
  +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
  
  Mig_route_rstr = 
    raster(
      Mig_route_grd
    )
  
  Mig_route_proj = 
    projectRaster(
      Mig_route_rstr,
      crs = geo_prj,
      res = 0.1
    )
  
  image(
    Mig_route_proj,
    col = 
      c(
        'white',
        'gray75'
      ),
    bty = 'n',
    xlab = NA,
    ylab = NA,
    xaxt = 'n',
    yaxt = 'n',
    add = TRUE
  )

  rm(COMPLETE_DATA)
  print(i)
}

plot(
  CONUS_NA_sts,
  lty = 'dotted',
  border = 'gray',
  add = TRUE
)

for(
  i in 1:(
    length(
      year_from_record
    ) - 1
  )
) {
  lines(
    CenterofMass_data[
      Year %like% 
        year_from_record[i],
      Longitude
      ],
    CenterofMass_data[
      Year %like% 
        year_from_record[i],
      Latitude
      ],
    lwd = 1.5,
    col = 'gray40'
  )
}

lines(
  CenterofMass_data[
    ,
    unique(
      mean_Longitude
    )
    ],
  CenterofMass_data[
    ,
    unique(
      mean_Latitude
    )
    ],
  lwd = 3,
  col = 'black'
)


#
# Plot spatial variation in WSI -------------------------------------------
gc()

spatial_wsi_unique = list()
for(
  i in 1:
    (
      length(
        year_from_record
      ) - 1
    )
) {
  COMPLETE_DATA =
    fread(
      paste0(
        "NodeSpecifData_Temperature_WSI_AirDensity",
        year_from_record[i],
        ".txt"
      )
    )
  
  spatial_wsi_dt = 
    COMPLETE_DATA[
      ,
      c(
        'X_INDEX',
        'Y_INDEX',
        'Code',
        names(
          COMPLETE_DATA
        )[
          names(
            COMPLETE_DATA
          ) %like%
            'WSI_day_'
        ]
      ),
      with = FALSE
    ]
  
  spatial_wsi_melt = 
    data.table::melt(
      spatial_wsi_dt,
      id.vars = 
        c(
          'X_INDEX',
          'Y_INDEX',
          'Code'
        )
    )
  
  spatial_wsi_melt[
    ,
    wsi_sd :=
      mean(
        value
      ),
    by = Code
  ]
  
  spatial_wsi_unique[[i]] = 
    unique(
      spatial_wsi_melt[
        ,
        .(
          X_INDEX,
          Y_INDEX,
          Code,
          wsi_sd
        )
      ]
    )[
      ,
      wsi_sd
    ]
  
  rm(COMPLETE_DATA)
  
  print(
    (
      i / 
        length(
          year_from_record
        ) * 100
    )
  )
}

COMPLETE_DATA =
  fread(
    paste0(
      'NodeSpecifData_Temperature_WSI_AirDensity',
      year_from_record[1],
      '.txt'
    )
  )

wsi_var_dt = 
  cbind(
    COMPLETE_DATA[
      ,
      .(
        X_INDEX,
        Y_INDEX
      )
    ],
    do.call(
      cbind,
      spatial_wsi_unique
    )
  )

wsi_var_melt = 
  data.table::melt(
    wsi_var_dt,
    id.vars = 
      c(
        'X_INDEX',
        'Y_INDEX'
      )
  )

wsi_var_melt[
  ,
  wsi_sd :=
    sd(
      value
    ),
  by = 
    .(
      X_INDEX,
      Y_INDEX
    )
]

wsi_var_unique = 
  unique(
    wsi_var_melt[
      ,
      .(
        X_INDEX,
        Y_INDEX,
        wsi_sd
      )
    ]
  )

wsi_var_plot = 
  ggplot(
    data = wsi_var_unique, 
    aes(
      x = X_INDEX,
      y = Y_INDEX, 
      col = wsi_sd
    )
  ) + 
  geom_point() +
  labs(
    col = 'Standard Deviation\nof Mean Annual\nWeather Severity Index'
  ) +
  scale_colour_gradient(
    low = 'black',
    high = 'white'
  ) +
  theme_void()

ggsave(
  'WSI_variation_plot.jpeg',
  wsi_var_plot,
  width = 6,
  height = 5
)
    

#
# Animations to visualize output over time --------------------------------
# Create masking work-flow to improve efficiency
COMPLETE_DATA =
  fread(
    paste0(
      year_from_record[1],
      '_Prediction.txt'
    )
  )

DayAbund =
  ZO_std(
    as.matrix(
      COMPLETE_DATA[
        ,
        N_Abund_0
        ]
    )
  )

DayNodeAbund =
  cbind(
    COMPLETE_DATA[
      ,
      .(
        Longitude,
        Latitude
      )
      ],
    DayAbund
  )

coordinates(
  DayNodeAbund
) =
  c(
    'Longitude',
    'Latitude'
  )

gridded(
  DayNodeAbund
) = TRUE

DayNodeAbund_grd =
  as(
    DayNodeAbund,
    'SpatialGridDataFrame'
  )

crs(
  DayNodeAbund_grd
) =
  '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

DayNodeAbund_rstr =
  raster(
    DayNodeAbund_grd
  )

DayNodeAbund_proj =
  projectRaster(
    DayNodeAbund_rstr,
    crs = geo_prj,
    res = 0.1
  )

dummy_rstr =
  raster(
    nrows =
      nrow(
        DayNodeAbund_proj
      ),
    ncols =
      ncol(
        DayNodeAbund_proj
      ),
    xmn =
      xmin(
        DayNodeAbund_proj
      ),
    xmx =
      xmax(
        DayNodeAbund_proj
      ),
    ymn =
      ymin(
        DayNodeAbund_proj
      ),
    ymx =
      ymax(
        DayNodeAbund_proj
      )
  )

dummy_rstr[] = 1

CONUS_NA_sts =
  readOGR(
    'states',
    dsn = getwd()
  )

dummy_mask =
  mask(
    dummy_rstr,
    CONUS_NA_sts
  )

arr_dummy_mask =
  as.array(
    dummy_mask
  )[
    ,
    ,
    1
    ]

arr_dummy_mask[
  which(
    is.na(
      arr_dummy_mask
    )
  )
  ] = 0

arr_dummy_mask[
  which(
    arr_dummy_mask != 0
  )
  ] = 1

storage.mode(
  arr_dummy_mask
) = 'logical'

# Read in North America Borders
CONUS_NA_sts_fort =
  ggplot2::fortify(
    CONUS_NA_sts
  )

# Transform coordinates
COMPLETE_DATA_sp =
  COMPLETE_DATA[
    ,
    .(
      Longitude,
      Latitude
    )
    ]

coordinates(
  COMPLETE_DATA_sp
) =
  c(
    'Longitude',
    'Latitude'
  )

crs(
  COMPLETE_DATA_sp
) =
  '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

COMPLETE_DATA_sp =
  spTransform(
    COMPLETE_DATA_sp,
    CRS =
      CRS(
        geo_prj
      )
  )

trnsfrmd_coords =
  COMPLETE_DATA_sp@coords

# Create daily abundance transformed coordinate data tables per year
Sys.time()
start_timer = Sys.time()

n_cores = detectCores() - 1
start_cluster = makeCluster(n_cores)
registerDoParallel(start_cluster)

foreach(
  i = 1:(
    length(
      year_from_record
    ) - 1
  ),
  .packages = 
    c(
      'sp',
      'data.table',
      'raster',
      'grDevices',
      'fields',
      'SDMTools'
    ),
  .verbose = TRUE
) %dopar% {
  COMPLETE_DATA =
    fread(
      paste0(
        getwd(),
        '/',
        year_from_record[i],
        '_Prediction.txt'
      )
    )
  
  COMPLETE_DATA[
    ,
    `:=`(
      trnsfrmd_Longitude = trnsfrmd_coords[,1],
      trnsfrmd_Latitude = trnsfrmd_coords[,2]
    )
    ]
  
  DayNodeAbund_dt = list()
  for(
    j in 1:num_days
  ) {
    which_day =
      paste0(
        'N_Abund_',
        j - 1
      )
    
    DayAbund =
      ZO_std(
        as.matrix(
          COMPLETE_DATA[
            ,
            get(
              which_day
            )
            ]
        )
      )
    
    DayNodeAbund =
      cbind(
        COMPLETE_DATA[
          ,
          .(
            Longitude,
            Latitude
          )
          ],
        DayAbund
      )
    
    coordinates(
      DayNodeAbund
    ) =
      c(
        'Longitude',
        'Latitude'
      )
    
    gridded(
      DayNodeAbund
    ) = TRUE
    
    DayNodeAbund_grd =
      as(
        DayNodeAbund,
        'SpatialGridDataFrame'
      )
    
    crs(
      DayNodeAbund_grd
    ) =
      '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
    +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
    
    DayNodeAbund_rstr =
      raster(
        DayNodeAbund_grd
      )
    
    DayNodeAbund_proj =
      projectRaster(
        DayNodeAbund_rstr,
        crs = geo_prj,
        res = 0.1
      )
    
    crop_bounds =
      extent(
        CONUS_NA_sts
      )
    
    DayNodeAbund_crp =
      crop(
        DayNodeAbund_proj,
        crop_bounds
      )
    
    DayNodeAbund_array =
      as.array(
        DayNodeAbund_proj
      )[
        ,
        ,
        1
        ]
    
    DayNodeAbund_array_masked =
      fun.mask(
        DayNodeAbund_array
      )
    
    DayNodeAbund_msk =
      raster(
        DayNodeAbund_array_masked,
        xmn =
          xmin(
            DayNodeAbund_proj
          ),
        xmx =
          xmax(
            DayNodeAbund_proj
          ),
        ymn =
          ymin( DayNodeAbund_proj
          ),
        ymx =
          ymax(
            DayNodeAbund_proj
          )
      )
    
    DayNodeAbund_spdf =
      as(
        DayNodeAbund_msk,
        'SpatialPixelsDataFrame'
      )
    
    DayNodeAbund_dt[[j]] =
      data.table(
        as.data.frame(
          DayNodeAbund_spdf
        )
      )
    
    DayNodeAbund_dt[[j]][
      ,
      N_Abund_Day :=
        j - 1
      ]
    
    DayNodeAbund_dt[[j]][
      ,
      Date :=
        as.Date(
          as.numeric(
            N_Abund_Day
          ),
          origin = NonbreedingStart[i]
        )
      ]
    
    target_cols =
      names(
        COMPLETE_DATA
      )[
        which(
          names(
            COMPLETE_DATA
          ) %like% 'N_Abund_'
        )
        ][
          2:(
            num_days + 1
          )
          ]
    
    CenterofMass =
      t(
        as.matrix(
          COMPLETE_DATA[
            ,
            lapply(
              .SD,
              Pop_CoM,
              x = trnsfrmd_Longitude,
              y = trnsfrmd_Latitude,
              z = NULL
            ),
            .SDcols =
              target_cols
            ][
              c(
                1,
                3
              )
              ]
        )
      )
    
    CenterofMass_pts =
      SpatialPoints(
        CenterofMass
      )
    
    crs(
      CenterofMass_pts
    ) =
      geo_prj
    
    CenterofMass_coords =
      data.table(
        coordinates(
          CenterofMass_pts
        )
      )
    
    CenterofMass_coords[
      ,
      ColorIterate :=
        c(
          1:.N
        )
      ]
    
    col_pal =
      colorRampPalette(
        c(
          'white',
          'black'
        )
      )
    
    CenterofMass_coords[
      ,
      ColorGradient :=
        col_pal(
          .N
        )[
          as.numeric(
            cut(
              ColorIterate,
              breaks = .N
            )
          )
          ]
      ]
    
    DayNodeAbund_dt[[j]][
      which.min(
        rdist(
          CenterofMass_coords[
            j,
            .(
              coords.x1,
              coords.x2
            )
            ],
          DayNodeAbund_dt[[j]][
            N_Abund_Day == j - 1,
            .(
              x,
              y
            )
            ]
        )
      ),
      CoM_pt :=
        CenterofMass_coords[
          j,
          ColorGradient
          ]
      ]
    
    print(
      (
        j /
          num_days
      ) * 100
    )
  }
  
  DayNodeAbund_dt =
    rbindlist(
      DayNodeAbund_dt
    )
  
  fwrite(
    DayNodeAbund_dt,
    paste0(
      getwd(),
      '/',
      year_from_record[i],
      '_Prediction_transformed.txt'
    )
  )
  
  print(
    (
      i /
        (
          length(
            year_from_record
          ) - 1
        )
    ) * 100
  )
}

stopCluster(start_cluster)

Sys.time()
end_timer = Sys.time()
total_time = end_timer - start_timer
total_time

# Animation for daily node abundance (0-1 scale)
Sys.time()
start_timer = Sys.time()

n_cores = detectCores() - 1
start_cluster = makeCluster(n_cores)
registerDoParallel(start_cluster)

foreach(
  i = 13:(
    length(
      year_from_record
    ) - 1
  ),
  .packages = 
    c(
      'ggplot2',
      'data.table',
      'viridis',
      'animation'
    ),
  .verbose = TRUE
) %dopar% {
  DayNodeAbund_dt =
    fread(
      paste0(
        getwd(),
        '/',
        year_from_record[i],
        '_Prediction_transformed.txt'
      )
    )
  
  DailyAbundCoM_Animate_func =
    function() {
      for(
        j in 1:num_days
      ) {
        print(
          ggplot() +
            geom_tile(
              data =
                DayNodeAbund_dt[
                  N_Abund_Day == j - 1
                  ],
              aes(
                x = x,
                y = y,
                fill = layer
              ),
              alpha = 0.8
            ) +
            geom_polygon(
              data = CONUS_NA_sts_fort,
              aes(
                x = long,
                y = lat,
                group = group
              ),
              fill = NA,
              color = 'gray80',
              size = 0.25
            ) +
            geom_point(
              data =
                DayNodeAbund_dt[
                  N_Abund_Day >= 0 &
                    N_Abund_Day <= j - 1 &
                    CoM_pt != ''
                  ],
              aes(
                x = x,
                y = y,
                col = CoM_pt
              )
            ) +
            scale_color_manual(
              values =
                DayNodeAbund_dt[
                  N_Abund_Day >= 0 &
                    N_Abund_Day <= j - 1 &
                    CoM_pt != '',
                  CoM_pt
                  ]
            ) +
            guides(
              color = FALSE
            ) +
            scale_fill_viridis(
              name =
                paste0(
                  'Normalized Abundance: \n',
                  DayNodeAbund_dt[
                    N_Abund_Day == j - 1,
                    unique(
                      Date
                    )
                    ]
                )
            ) +
            theme(
              line = element_blank(),
              rect = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 16)
            )
        )
        
        print((j / num_days) * 100)
      }
    }
  
  ani.options(
    ani.width = 600,
    ani.height = 480
  )
  
  saveVideo(
    DailyAbundCoM_Animate_func(),
    interval = 0.1,
    video.name =
      paste0(
        'DailyAbundCoM_',
        year_from_record[i],
        '.mp4'
      ),
    ffmpeg = 'C:/ffmpeg/bin/ffmpeg.exe'
  )
  
  print(
    (
      i / (
        length(
          year_from_record
        ) - 1
      )
    ) * 100
  )
}

stopCluster(start_cluster)

Sys.time()
end_timer = Sys.time()
total_time = end_timer - start_timer
total_time

# Animation for daily body condition class proportional abundance
Sys.time()
start_timer = Sys.time()

n_cores = detectCores() - 1
start_cluster = makeCluster(n_cores)
registerDoParallel(start_cluster)

foreach(
  i = 1:(
    length(
      year_from_record
    ) - 1
  ),
  .packages = 
    c(
      'ggplot2',
      'data.table',
      'animation'
    ),
  .verbose = TRUE
) %dopar% {
  BC_Table = 
    fread(
      paste0(
        year_from_record[j],
        '_BCTable.txt'
      )
    )
  
  DailyBCPAAnimate_func = 
    function(
      
    ) {
      for(
        i in 1:length(
          days
        )
      ) {
        DayProportion =
          ZO_std(
            as.matrix(
              BC_Table[
                ,
                paste0(
                  'BodyMass_dist_discrete_day',
                  i
                ),
                with = FALSE
                ]
            )
          )
        
      print(
        ggplot(
          data = BC_Table,
               aes(
                 x = BC_bins,
                 y = DayProportion
                 )
          ) +
          geom_line(
            lwd = 1.2
            ) +
          xlab('Body Condition Class') +
          ylab('Normalized Proportional Abundance') +
          ggtitle(
            paste0(
              format(
                as.Date(
                  i,
                  origin = NonbreedingStart
                  ),
                                format = '%b %d'
                )
              )
            ) +
          theme_bw() +
          theme(
            axis.text.x = 
              element_text(
                size = 16,
                face = 'bold',
                colour = 'black'
                ),
            axis.text.y = 
              element_text(
                size = 16,
                           face = 'bold',
                           colour = 'black'
              ),
            axis.title = 
              element_text(
                size = 24,
                face = 'bold',
                colour = 'black'
                ),
            plot.title = 
              element_text(
              size = 24,
              face = 'bold',
              colour = 'black'
            )
            
          )
      )
      
      print(
        (
          i / 
            length(
              days
            )
        ) * 100
      )
    }
  }
  
  saveVideo(
    DailyBCPAAnimate_func(),
    interval = 0.1,
    video.name =
      paste0(
        'DailyBCPA',
        year_from_record[j],
        '.mp4'
      )
  )
  
  print(j)
}

stopCluster(start_cluster)

Sys.time()
end_timer = Sys.time()
total_time = end_timer - start_timer
total_time


#