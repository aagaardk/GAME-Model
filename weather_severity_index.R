########################################################################################
# WEATHER SEVERITY INDEX 
# for,
# FLYWAY LEVEL WATERFOWL MIGRATION MODEL
#---------------------------------------------------------------------------------------
#
# Created by Kevin Aagaard in collaboration with Eric Lonsdorf, Sarah Jacobi, 
# Wayne Thogmartin, and Tim Jones, and the Integrated Waterbird Management and 
# Monitoring Program
#
# Modified: 
Sys.time()
# By:        Kevin Aagaard
#
# **ENSURE THAT ALL FILES ARE IN THE 'R' FOLDER FOR PROPER WORKING DIRECTORY FUNCTION**
#
########################################################################################

# Set directory -----------------------------------------------------------
#Set working directory (should be universal)
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path)) 
getwd()

# Clear workspace and load packages and libraries -------------------------
rm(list=ls())

# # For use at UMESC when library path directory is not working (permissioning issue)
# .libPaths(.libPaths()[2])

# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
x = c("timeDate","lubridate","data.table","plyr","lattice","raster","rgeos",
      "maptools","foreign","xts","limma","PerformanceAnalytics","compiler",
      "fields","rpart","RNCEP","ggmap","rgdal","dplyr","tidyr","tmap", 
      "geosphere","tgp","animation","colorspace","colorRamps")
# install.packages(x) #THIS MAY TAKE A WHILE
lapply(x, library, character.only=TRUE)
rm(x)

# Set some global conditions ----------------------------------------------
enableJIT(3)
timeit = proc.time()
options(digits=16)

# Read in global data -----------------------------------------------------
NODE_DATA = 
  fread(
    'NA_NODE_LC2.txt'
  )

# # Columns in NODE_DATA:
# # 1:FID; 2:Y-index; 3:X-index; 4:code; 5:long; 6:lat; 7:shoreline; 8: open
# # water; 9: ice/snow; 10:open space 11:low-int devel; 12:med-intense
# # devel; 13:high-int develop; 14:barren; 15:Decid Forest; 16:Evergreen
# # forest; 17:mixed forest; 18:dwarf scrub; 19:Shrub/scrub 20:herb grass;
# # 21:sedge/herb; 22:moss; 23:pasture/hay; 24:crops; 25:woody wetlands;
# # 26:herb wetlands; 27:mean dist to crops; 28:std to crops; 29:mean woody
# # wet; 30:std wood wet; 31:mean herb wet; 32:std herb wet
# 
# # To create IWMM 20-mile node dataset:
# NODE_LOCATIONS = 
#   fread(
#     'node_locations.csv'
#   )
# 
# # Add lat/long decimal degrees to NODE_DATA
# setkey(
#   NODE_DATA,
#   Y_INDEX,
#   X_INDEX
# )
# 
# setkey(
#   NODE_LOCATIONS,
#   Y_INDEX,
#   X_INDEX
# )
# 
# Twenty_mi_nodes =
#   NODE_LOCATIONS[
#     NODE_DATA
#   ]
# 
# setnames(
#   Twenty_mi_nodes,
#   c(
#     'LATITUDE',
#     'LONGITUDE',
#     'Long',
#     'Lat'
#   ),
#   c(
#     'LatDegree',
#     'LongDegree',
#     'LongUTM',
#     'LatUTM'
#   )
# )
# 
# Twenty_mi_nodes[
#   ,
#   Y_INDEX :=
#     0 - (
#       Y_INDEX - 
#         max(
#           Y_INDEX
#         )
#     )
# ]
# 
# # Save Twenty-mile Nodes
# fwrite(
#   Twenty_mi_nodes,
#   file = 'Twenty_mile_node_locations.csv'
# )

Twenty_mi_nodes = 
  fread(
    'Twenty_mile_node_locations.csv'
  )

# Global variables --------------------------------------------------------
num_nodes = nrow(NODE_DATA)
node_move = 3000
hab_dist = matrix(0,num_nodes,3)
num_years = 1
num_flocks = 1
num_jumps = 8

min_year = 1957
max_year = 2019

# Distance effects --------------------------------------------------------
NODE_DATA = 
  data.frame(
    NODE_DATA
  )

hab_dist[,1] = 
  pnorm(
    node_move,
    NODE_DATA[,27],
    sapply(
      NODE_DATA[,28],
      function(
        x
      ) {
        max(
          x,
          30
        )
      }
    )
  )

hab_dist[,2] = 
  pnorm(
    node_move,
    NODE_DATA[,29],
    sapply(
      NODE_DATA[,30],
      function(x) {
        max(
          x,
          30
        )
      }
    )
  )

hab_dist[,3] = 
  pnorm(
    node_move,
    NODE_DATA[,31],
    sapply(
      NODE_DATA[,32],
      function(x) {
        max(
          x,
          30
        )
      }
    )
  )

hab_dist[which(NODE_DATA[,24]==0),1] = 0
hab_dist[which(NODE_DATA[,25]==0),2] = 0
hab_dist[which(NODE_DATA[,26]==0),3] = 0

coords = NODE_DATA[,2:3]
bird_hab = cbind(as.numeric(NODE_DATA[,7]),
                  NODE_DATA[,24:26]*hab_dist)

rm(hab_dist)
bird_hab = bird_hab/900

cals_temp = rowSums(bird_hab)*33945

# Load mallard distribution data (BPOP and NatureServe) -------------------
S = 
  read.dbf(
    'NorthAmerica_20mi_grid_wAK_BPOP_NSmallard_join.dbf'
  )

S.lat_long = 
  merge(
    S, 
    Twenty_mi_nodes
  )

# Download data -----------------------------------------------------------
# # The following code need only be run once to create the files. If the files
# # are lost, this will recreate them. Once they have been saved to the wd(),
# # there is no need to run this.
# 
# Temp.extent =
#   NCEP.gather(
#     variable = 'air.2m',
#     level = 'gaussian',
#     months.minmax = c(1,12),
#     years.minmax = c(min_year,max_year),
#     lat.southnorth =
#       c(
#         min(Twenty_mi_nodes$LatDegree),
#         max(Twenty_mi_nodes$LatDegree)
#       ),
#     lon.westeast =
#       c(
#         min(Twenty_mi_nodes$LongDegree),
#         max(Twenty_mi_nodes$LongDegree)
#       ),
#     reanalysis2 = FALSE,
#     return.units = TRUE
#   )
# 
SnowDepth.extent =
  NCEP.gather(
    variable = 'weasd.sfc',
    level = 'gaussian',
    months.minmax =
      c(
        1,
        12
      ),
    years.minmax =
      c(
        min_year,
        max_year
      ),
    lat.southnorth =
      c(
        min(Twenty_mi_nodes$LatDegree),
        max(Twenty_mi_nodes$LatDegree)
      ),
    lon.westeast =
      c(
        min(Twenty_mi_nodes$LongDegree),
        max(Twenty_mi_nodes$LongDegree)
      ),
    reanalysis2 = FALSE,
    return.units = TRUE
  )
# 
# AirDensity.extent =
#   NCEP.gather(
#     variable = 'pres.lcb',
#     level = 'gaussian',
#     months.minmax = c(1,12),
#     years.minmax = c(min_year,max_year),
#     lat.southnorth =
#       c(
#         min(Twenty_mi_nodes$LatDegree),
#         max(Twenty_mi_nodes$LatDegree)
#       ),
#     lon.westeast =
#       c(
#         min(Twenty_mi_nodes$LongDegree),
#         max(Twenty_mi_nodes$LongDegree)
#       ),
#     reanalysis2 = FALSE,
#     return.units = TRUE
#   )
# 
# # # Don't know how to incorporate wind velocity yet
# # EW_WindVelocity.extent =
# #   NCEP.gather(
# #     variable = 'uwnd',
# #     level = 'gaussian',
# #     months.minmax = c(1,12),
# #     years.minmax = c(min_year,max_year),
# #     lat.southnorth =
# #       c(
# #         min(Twenty_mi_nodes$LatDegree),
# #         max(Twenty_mi_nodes$LatDegree)
# #       ),
# #     lon.westeast =
# #       c(
# #         min(Twenty_mi_nodes$LongDegree),
# #         max(Twenty_mi_nodes$LongDegree)
# #       ),
# #     reanalysis2 = FALSE,
# #     return.units = TRUE
# #   )
# # 
# # NS_WindVelocity.extent =
# #   NCEP.gather(
# #     variable = 'vwnd',
# #     level = 'gaussian',
# #     months.minmax = c(1,12),
# #     years.minmax = c(min_year,max_year),
# #     lat.southnorth =
# #       c(
# #         min(Twenty_mi_nodes$LatDegree),
# #         max(Twenty_mi_nodes$LatDegree)
# #       ),
# #     lon.westeast =
# #       c(
# #         min(Twenty_mi_nodes$LongDegree),
# #         max(Twenty_mi_nodes$LongDegree)
# #       ),
# #     reanalysis2 = FALSE,
# #     return.units = TRUE
# #   )

# Restrict timeframe to non-breeding period -------------------------------
# Temp.extent =
#   NCEP.restrict(
#     Temp.extent,
#     months2remove = 6,
#     set2na = FALSE
#   )
# 
SnowDepth.extent =
  NCEP.restrict(
    SnowDepth.extent,
    months2remove = 6,
    set2na = FALSE
  )
# 
# AirDensity.extent =
#   NCEP.restrict(
#     AirDensity.extent,
#     months2remove = 6,
#     set2na = FALSE
#   )
# 
# Perform conversions -----------------------------------------------------
# # Convert air temperature from Kelvins to Celsius
# Temp.extent = Temp.extent - 273.15
#
# Convert snow data from snow-water-equivalence to depth
# Snow density ranges from 10 to 400 in our conditions. We're using 257
# to account for extreme conditions near arctic.
# http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings/parameters/
# snow_water_equivalent.html
density = 257

SnowDepth.extent = SnowDepth.extent / density

# Convert snow depth to centimeters from meters
SnowDepth.extent = SnowDepth.extent * 100

# # Convert air pressure in Pa to air density in kg / m^3
# # Specific gas constant for dry air is 287.058 (J / (kg*K))
# AirDensity.extent =
#   AirDensity.extent / (287.058 * (Temp.extent + 273.15))
# 
# Resample data onto regular grid -----------------------------------------
# Temp.resamp = array(data=NA,dim=dim(Temp.extent),dimnames(Temp.extent))
# pb = winProgressBar(title="Resampling and interpolating onto regular grid",
#                     label="0% done",min=0,max=100,initial=0)
# for(i in 1:dim(Temp.extent)[3]){
#   Temp.reg=interp.loess(x=rep(as.numeric(dimnames(Temp.extent)[[2]]),
#                               each=length(dimnames(Temp.extent)[[1]])),
#                         y=rep(as.numeric(dimnames(Temp.extent)[[1]]),
#                               length(dimnames(Temp.extent)[[2]])),
#                         z=as.vector(Temp.extent[,,i]), span=0.6,
#                         gridlen=c(length(dimnames(Temp.extent)[[2]]),
#                                   length(dimnames(Temp.extent)[[1]])))
# 
#   Temp.mat=matrix(data=t(Temp.reg$z),nrow=length(Temp.reg$y),ncol=length(Temp.reg$x))
#   Temp.mat=apply(Temp.mat, 2, rev)
# 
#   Temp.resamp[,,i]=Temp.mat
# 
#   info = sprintf("%d%% done", round((i/dim(Temp.extent)[3])*100))
#   setWinProgressBar(pb, i/(dim(Temp.extent)[3])*100, label=info)
# }
# close(pb)

SnowDepth.resamp = array(data=NA,dim=dim(SnowDepth.extent),dimnames(SnowDepth.extent))
pb = winProgressBar(title="Resampling and interpolating onto regular grid",
                    label="0% done",min=0,max=100,initial=0)
for(i in 1:dim(SnowDepth.extent)[3]){
  SnowDepth.reg=interp.loess(x=rep(as.numeric(dimnames(SnowDepth.extent)[[2]]),
                                   each=length(dimnames(SnowDepth.extent)[[1]])),
                             y=rep(as.numeric(dimnames(SnowDepth.extent)[[1]]),
                                   length(dimnames(SnowDepth.extent)[[2]])),
                             z=as.vector(SnowDepth.extent[,,i]), span=0.75,
                             gridlen=c(length(dimnames(SnowDepth.extent)[[2]]),
                                       length(dimnames(SnowDepth.extent)[[1]])))

  SnowDepth.mat=matrix(data=t(SnowDepth.reg$z),nrow=length(SnowDepth.reg$y),
                       ncol=length(SnowDepth.reg$x))
  SnowDepth.mat=apply(SnowDepth.mat, 2, rev)

  SnowDepth.resamp[,,i]=SnowDepth.mat

  info = sprintf("%d%% done", round((i/dim(SnowDepth.extent)[3])*100))
  setWinProgressBar(pb, i/(dim(SnowDepth.extent)[3])*100, label=info)
}
close(pb)

# AirDensity.resamp =
#   array(
#     data = NA,
#     dim = dim(AirDensity.extent),
#     dimnames(AirDensity.extent)
#   )
# pb =
#   winProgressBar(
#     title = "Resampling and interpolating onto regular grid",
#     label = "0% done",
#     min = 0,
#     max = 100,
#     initial = 0
#   )
# for(i in 1:dim(AirDensity.extent)[3]){
#   AirDensity.reg =
#     interp.loess(
#       x =
#         rep(
#           as.numeric(
#             dimnames(
#               AirDensity.extent
#             )[[2]]
#           ),
#           each =
#             length(
#               dimnames(
#                 AirDensity.extent
#               )[[1]]
#             )
#         ),
#       y =
#         rep(
#           as.numeric(
#             dimnames(
#               AirDensity.extent
#             )[[1]]
#           ),
#           length(
#             dimnames(
#               AirDensity.extent
#             )[[2]]
#           )
#         ),
#       z =
#         as.vector(
#           AirDensity.extent[,,i]
#         ),
#       span = 0.75,
#       gridlen =
#         c(
#           length(
#             dimnames(
#               AirDensity.extent
#             )[[2]]
#           ),
#           length(
#             dimnames(
#               AirDensity.extent
#             )[[1]]
#           )
#         )
#     )
# 
#   AirDensity.mat =
#     matrix(
#       data =
#         t(
#           AirDensity.reg$z
#         ),
#       nrow =
#         length(
#           AirDensity.reg$y
#         ),
#       ncol =
#         length(
#           AirDensity.reg$x
#         )
#     )
#   AirDensity.mat =
#     apply(
#       AirDensity.mat,
#       2,
#       rev
#     )
# 
#   AirDensity.resamp[,,i] = AirDensity.mat
# 
#   info =
#     sprintf(
#       "%d%% done",
#       round(
#         (
#           i / dim(AirDensity.extent)[3]
#         ) * 100
#       )
#     )
#   setWinProgressBar(
#     pb,
#     i / (dim(AirDensity.extent)[3]) * 100,
#     label = info
#   )
# }
# close(pb)
# 
# Convert to data.table via data.frame ------------------------------------
# Temp.df = NCEP.array2df(Temp.resamp)
# Temp.df.cols=colnames(Temp.df)
# 
SnowDepth.df = NCEP.array2df(SnowDepth.resamp)
SnowDepth.df.cols=colnames(SnowDepth.df)

# AirDensity.df = NCEP.array2df(AirDensity.resamp)
# AirDensity.df.cols = colnames(AirDensity.df)
# 
# Save tables -------------------------------------------------------------
work_dir = getwd()
# 
# # Air temperature
# filename = '/NA_air_temperature_data.txt'
# data_saver = file.path(paste(work_dir,filename,sep=''))
# 
# numrows = 1000000
# chunksize = floor(nrow(Temp.df)/numrows)
# chunks = 0:(chunksize-1)
# chunkseq = chunks*numrows
# pb = winProgressBar(title='Saving Air Temperature Data',
#                     label='0% done',min=0,max=100,initial=0)
# for (i in 1:chunksize){
#   Temp.df_chunk=Temp.df[(1+chunkseq[i]):((1+chunkseq[i])+(numrows-1)),]
#   write.table(Temp.df_chunk,file=data_saver,append=T,row.names=F,col.names=F)
# 
#     info = sprintf('%d%% done', round((i/chunksize)*100))
#   setWinProgressBar(pb, i/(chunksize)*100, label=info)
# }
# close(pb)
# 
# Temp.df_remainder=Temp.df[((chunksize*numrows)+1):nrow(Temp.df),]
# write.table(Temp.df_remainder,file=data_saver,append=T,row.names=F,col.names=F)
# 
# Snow depth
filename='/NA_snowdepth_data.txt'
data_saver=file.path(paste(work_dir,filename,sep=''))

numrows=1000000
chunksize=floor(nrow(SnowDepth.df)/numrows)
chunks=0:(chunksize-1)
chunkseq=chunks*numrows
pb = winProgressBar(title='Saving Snow Depth Data',
                    label='0% done',min=0,max=100,initial=0)
for (i in 1:chunksize){
  SnowDepth.df_chunk=SnowDepth.df[(1+chunkseq[i]):((1+chunkseq[i])+(numrows-1)),]
  write.table(SnowDepth.df_chunk,file=data_saver,append=T,row.names=F,col.names=F)

  info = sprintf('%d%% done', round((i/chunksize)*100))
  setWinProgressBar(pb, i/(chunksize)*100, label=info)
}
close(pb)

SnowDepth.df_remainder=SnowDepth.df[((chunksize*numrows)+1):nrow(SnowDepth.df),]
write.table(SnowDepth.df_remainder,file=data_saver,append=T,row.names=F,col.names=F)

# # Air density
# filename = '/NA_airdensity_data.txt'
# data_saver =
#   file.path(
#     paste(
#       work_dir,
#       filename,
#       sep = ""
#     )
#   )
# 
# numrows = 1000000
# chunksize =
#   floor(
#     nrow(
#       AirDensity.df
#     ) / numrows
#   )
# chunks = 0:(chunksize-1)
# chunkseq = chunks * numrows
# pb =
#   winProgressBar(
#     title = "Saving Air Density Data",
#     label = "0% done",
#     min = 0,
#     max = 100,
#     initial = 0
#   )
# for (i in 1:chunksize){
#   AirDensity.df_chunk =
#     AirDensity.df[
#       (1 + chunkseq[i]):((1 + chunkseq[i]) + (numrows - 1)),
#       ]
#   write.table(
#     AirDensity.df_chunk,
#     file = data_saver,
#     append = TRUE,
#     row.names = FALSE,
#     col.names = FALSE
#   )
# 
#   info =
#     sprintf(
#       "%d%% done",
#       round(
#         (i / chunksize) * 100
#       )
#     )
#   setWinProgressBar(
#     pb,
#     i / (chunksize) * 100,
#     label = info
#   )
# }
# close(pb)
# 
# AirDensity.df_remainder =
#   AirDensity.df[
#     ((chunksize * numrows) + 1):nrow(AirDensity.df),
#     ]
# write.table(
#   AirDensity.df_remainder,
#   file = data_saver,
#   append = TRUE,
#   row.names = FALSE,
#   col.names = FALSE
# )

# Create historical and averaged datasets ---------------------------------
# Convert downloaded data from frames to arrays ---------------------------
# Make room for large files
rm(list=ls())
gc()

# Read in templates
NODE_DATA = fread("NA_NODE_LC2.txt")

NODE_DATA[
  ,
  Y_INDEX :=
    0 - (
      Y_INDEX -
        max(
          Y_INDEX
        )
    )
]

Twenty_mi_nodes =
  fread(
    'Twenty_mile_node_locations.csv'
  )

gc()

# # Read in previously created data frames; this takes ~ 7 mins
# Temp.dt =
#   fread(
#     'NA_air_temperature_data.txt'
#   )
# 
Snow.dt =
  fread(
    'NA_snowdepth_data.txt'
  )

# AirDensity.dt =
#   fread(
#     'NA_airdensity_data.txt'
#   )
# 
# # Rename columns
# setnames(
#   Temp.dt,
#   names(
#     Temp.dt
#   ),
#   c(
#     'date',
#     'NCEP_LATITUDE',
#     'NCEP_LONGITUDE',
#     'daily_temp'
#   )
# )
# 
# Temp.dt[
#   ,
#   date :=
#     substr(
#       date,
#       1,
#       10
#     )
# ]
# 
# Rename columns
setnames(
  Snow.dt,
  names(
    Snow.dt
  ),
  c('date',
    'NCEP_LATITUDE','NCEP_LONGITUDE',
    'daily_snow'
  )
)
Snow.dt[
  ,
  date :=
    substr(
      date,
      1,
      10
    )
]

# Correct negative snow depth values
Snow.dt[
  ,
  daily_snow :=
    ifelse(
      daily_snow < 0,
      0,
      daily_snow
    )
]

# # Rename columns
# setnames(
#   AirDensity.dt,
#   names(
#     AirDensity.dt
#   ),
#   c(
#     'date',
#     'NCEP_LATITUDE',
#     'NCEP_LONGITUDE',
#     'daily_airdensity'
#   )
# )
# 
# AirDensity.dt[
#   ,
#   date :=
#     substr(
#       date,
#       1,
#       10
#     )
#   ]
# 
# Create historical data sets ---------------------------------------------
# # Air Temperature
# HistTemp.dt =
#   Temp.dt[
#     ,
#     round(
#       mean(
#         daily_temp
#       ),
#       digits = 3
#     ),
#     by =
#       .(
#         date,
#         NCEP_LATITUDE,
#         NCEP_LONGITUDE
#       )
#   ]
# 
# rm(
#   Temp.dt
# )
# 
# setnames(
#   HistTemp.dt,
#   'V1',
#   'daily_temp'
# )
# 
# HistTemp.dt[
#   ,
#   year :=
#     substr(
#       date,
#       1,
#       4
#     )
# ]
# 
# HistTemp.dt =
#   HistTemp.dt[
#     ,
#     .(
#       date,
#       year,
#       NCEP_LATITUDE,
#       NCEP_LONGITUDE,
#       daily_temp
#     )
#   ]
# 
# daily_mean_temperature =
#   dcast(
#     HistTemp.dt,
#     NCEP_LATITUDE + NCEP_LONGITUDE ~ date,
#     value.var = 'daily_temp'
#   )
# 
# setkey(
#   daily_mean_temperature,
#   NCEP_LATITUDE,
#   NCEP_LONGITUDE
# )
# 
# gc()
# 
# Snow Depth
HistSnow.dt =
  Snow.dt[
    ,
    round(
      mean(
        daily_snow
      ),
      digits = 3
    ),
    by =
      .(
        date,
        NCEP_LATITUDE,
        NCEP_LONGITUDE
      )
  ]

rm(
  Snow.dt
)

setnames(
  HistSnow.dt,
  'V1',
  'daily_snow'
)

HistSnow.dt[
  ,
  year :=
    substr(
      date,
      1,
      4
    )
]

HistSnow.dt =
  HistSnow.dt[
    ,
    .(
      date,
      year,
      NCEP_LATITUDE,
      NCEP_LONGITUDE,
      daily_snow
    )
  ]

daily_mean_snow =
  dcast(
    HistSnow.dt,
    NCEP_LATITUDE + NCEP_LONGITUDE ~ date,
    value.var = 'daily_snow'
  )

setkey(
  daily_mean_snow,
  NCEP_LATITUDE,
  NCEP_LONGITUDE
)

gc()

# # Air Density
# HistAirDensity.dt =
#   AirDensity.dt[
#     ,
#     round(
#       mean(
#         daily_airdensity
#       ),
#       digits = 3
#     ),
#     by =
#       .(
#         date,
#         NCEP_LATITUDE,
#         NCEP_LONGITUDE
#       )
#     ]
# 
# rm(
#   AirDensity.dt
# )
# 
# setnames(
#   HistAirDensity.dt,
#   'V1',
#   'daily_airdensity'
# )
# 
# HistAirDensity.dt[
#   ,
#   year :=
#     substr(
#       date,
#       1,
#       4
#     )
#   ]
# HistAirDensity.dt =
#   HistAirDensity.dt[
#     ,
#     .(
#       date,
#       year,
#       NCEP_LATITUDE,
#       NCEP_LONGITUDE,
#       daily_airdensity
#     )
#     ]
# 
# daily_mean_airdensity =
#   dcast(
#     HistAirDensity.dt,
#     NCEP_LATITUDE + NCEP_LONGITUDE ~ date,
#     value.var = 'daily_airdensity'
#   )
# setkey(
#   daily_mean_airdensity,
#   NCEP_LATITUDE,
#   NCEP_LONGITUDE
# )
# 
# gc()
# 
# Define nodes based on distance from node centroid to sampling location
Twenty_mi_nodes_strpd =
  Twenty_mi_nodes[
    ,
    .(
      Y_INDEX,
      X_INDEX,
      code,
      LatDegree,
      LongDegree
    )
    ]
setnames(
  Twenty_mi_nodes_strpd,
  old =
    colnames(Twenty_mi_nodes_strpd),
  new =
    c(
      "Y_INDEX",
      "X_INDEX",
      "NODE",
      "TMN_LATITUDE",
      "TMN_LONGITUDE"
    )
)
Twenty_mi_nodes_strpd[
  ,
  TMN_LONGITUDE :=
    TMN_LONGITUDE + 360
  ]

node_locations =
  Twenty_mi_nodes_strpd[
    ,
    .(
      TMN_LONGITUDE,
      TMN_LATITUDE
    )
    ]
node_locations =
  as.matrix(
    node_locations,
    nrow = nrow(Twenty_mi_nodes),
    ncol = 2
  )

# # Air Temperature
# samp_locations =
#   HistTemp.dt[
#     ,
#     .(
#       NCEP_LONGITUDE,
#       NCEP_LATITUDE
#     )
#   ]
# 
# samp_locations =
#   unique(
#     samp_locations
#   )
# 
# samp_locations =
#   as.matrix(
#     samp_locations,
#     nrow =
#       nrow(
#         HistTemp.dt
#       ),
#     ncol = 2
#   )
# 
# closest_samp =
#   matrix(
#     NA,
#     nrow =
#       nrow(
#         node_locations
#       ),
#     ncol = 4
#   )
# 
# closest_samp[
#   ,
#   1:2
# ] =
#   node_locations
# 
# for(
#   i in 1:nrow(
#     node_locations
#   )
# ) {
#   closest_samp[
#     i,
#     3:4
#   ]=
#     samp_locations[
#       which.min(
#         distHaversine(
#           node_locations[
#             i,
#           ],
#           samp_locations
#         )
#       ),
#     ]
# }
# 
# colnames(
#   closest_samp
# ) =
#   c(
#     'TMN_LONGITUDE',
#     'TMN_LATITUDE',
#     'NCEP_LONGITUDE',
#     'NCEP_LATITUDE'
#   )
# 
# closest_samp =
#   data.table(
#     closest_samp
#   )
# 
# TMN_daily_mean_temperature =
#   data.table(
#     matrix(
#       0,
#       nrow =
#         nrow(
#           closest_samp
#         ),
#       ncol =
#         ncol(
#           daily_mean_temperature
#         )
#     )
#   )
# 
# TMN_daily_mean_temperature[
#   ,
#   `:=`(
#     V1 =
#       closest_samp[
#         ,
#         NCEP_LATITUDE
#       ],
#     V2 =
#       closest_samp[
#         ,
#         NCEP_LONGITUDE
#       ]
#   )
# ]
# 
# setnames(
#   TMN_daily_mean_temperature,
#   names(
#     TMN_daily_mean_temperature
#   ),
#   colnames(
#     daily_mean_temperature
#   )
# )
# 
# setkey(
#   TMN_daily_mean_temperature,
#   NCEP_LATITUDE,
#   NCEP_LONGITUDE
# )
# 
# TMN_daily_mean_temperature =
#   daily_mean_temperature[
#     TMN_daily_mean_temperature
#   ][
#     ,
#     1:ncol(
#       daily_mean_temperature
#     ),
#     with = FALSE
#   ]
# 
# setkey(
#   closest_samp,
#   TMN_LONGITUDE,
#   TMN_LATITUDE
# )
# 
# setkey(
#   Twenty_mi_nodes_strpd,
#   TMN_LONGITUDE,
#   TMN_LATITUDE
# )
# 
# full_node_coords =
#   closest_samp[
#     Twenty_mi_nodes_strpd
#   ]
# 
# setkey(
#   TMN_daily_mean_temperature,
#   NCEP_LONGITUDE,
#   NCEP_LATITUDE
# )
# 
# setkey(
#   full_node_coords,
#   NCEP_LONGITUDE,
#   NCEP_LATITUDE
# )
# 
# TMN_daily_mean_temperature[
#   ,
#   `:=`(
#     NCEP_LATITUDE = NULL,
#     NCEP_LONGITUDE = NULL
#   )
# ]
# 
# TMN_daily_mean_temperature =
#   cbind(
#     full_node_coords[
#       ,
#       .(
#         Y_INDEX,
#         X_INDEX
#       )
#     ],
#     TMN_daily_mean_temperature
#   )
# 
# TMN_daily_mean_temperature =
#   as.matrix(
#     TMN_daily_mean_temperature
#   )
# 
# Snow Depth
samp_locations =
  HistSnow.dt[
    ,
    .(
      NCEP_LONGITUDE,
      NCEP_LATITUDE
    )
  ]

samp_locations =
  unique(
    samp_locations
  )

samp_locations =
  as.matrix(
    samp_locations,
    nrow =
      nrow(
        HistSnow.dt
      ),
    ncol = 2
  )

closest_samp =
  matrix(
    NA,
    nrow =
      nrow(
        node_locations
      ),
    ncol = 4
  )

closest_samp[
  ,
  1:2
] =
  node_locations

for(
  i in 1:nrow(
    node_locations
  )
) {
  closest_samp[
    i,
    3:4
  ] =
    samp_locations[
      which.min(
        distHaversine(
          node_locations[
            i,
          ],
          samp_locations
        )
      ),
    ]
}

colnames(
  closest_samp
) =
  c(
    'TMN_LONGITUDE',
    'TMN_LATITUDE',
    'NCEP_LONGITUDE',
    'NCEP_LATITUDE'
  )

closest_samp =
  data.table(
    closest_samp
  )

TMN_daily_mean_snow =
  data.table(
    matrix(
      0,
      nrow =
        nrow(
          closest_samp
        ),
      ncol =
        ncol(
          daily_mean_snow
        )
    )
  )

TMN_daily_mean_snow[
  ,
  `:=`(
    V1 =
      closest_samp[
        ,
        NCEP_LATITUDE
      ],
    V2 =
      closest_samp[
        ,
        NCEP_LONGITUDE
      ]
  )
]

setnames(
  TMN_daily_mean_snow,
  names(
    TMN_daily_mean_snow
  ),
  colnames(
    daily_mean_snow
  )
)

setkey(
  TMN_daily_mean_snow,
  NCEP_LATITUDE,
  NCEP_LONGITUDE
)

TMN_daily_mean_snow=
  daily_mean_snow[
    TMN_daily_mean_snow
  ][
    ,
    1:ncol(
      daily_mean_snow
    ),
    with = FALSE
  ]

setkey(
  closest_samp,
  TMN_LONGITUDE,
  TMN_LATITUDE
  )

setkey(
  Twenty_mi_nodes_strpd,
  TMN_LONGITUDE,
  TMN_LATITUDE
)

full_node_coords =
  closest_samp[
    Twenty_mi_nodes_strpd
  ]

setkey(
  TMN_daily_mean_snow,
  NCEP_LONGITUDE,
  NCEP_LATITUDE
)

setkey(
  full_node_coords,
  NCEP_LONGITUDE,
  NCEP_LATITUDE
)

TMN_daily_mean_snow[
  ,
  `:=`(
    NCEP_LATITUDE = NULL,
    NCEP_LONGITUDE = NULL
  )
]

TMN_daily_mean_snow =
  cbind(
    full_node_coords[
      ,
      .(
        Y_INDEX,
        X_INDEX
      )
    ],
    TMN_daily_mean_snow
  )

TMN_daily_mean_snow =
  as.matrix(
    TMN_daily_mean_snow
  )

# # Air Density
# samp_locations =
#   HistAirDensity.dt[
#     ,
#     .(
#       NCEP_LONGITUDE,
#       NCEP_LATITUDE
#     )
#     ]
# samp_locations =
#   unique(
#     samp_locations
#   )
# samp_locations =
#   as.matrix(
#     samp_locations,
#     nrow =
#       nrow(HistAirDensity.dt),
#     ncol = 2
#   )
# 
# closest_samp =
#   matrix(
#     NA,
#     nrow = nrow(node_locations),
#     ncol = 4
#   )
# closest_samp[,1:2] = node_locations
# for(i in 1:nrow(node_locations)){
#   closest_samp[i,3:4]=
#     samp_locations[
#       which.min(
#         distHaversine(
#           node_locations[i,],
#           samp_locations
#         )
#       ),
#       ]
# }
# 
# colnames(closest_samp) =
#   c(
#     "TMN_LONGITUDE",
#     "TMN_LATITUDE",
#     "NCEP_LONGITUDE",
#     "NCEP_LATITUDE"
#   )
# closest_samp =
#   data.table(
#     closest_samp
#   )
# 
# TMN_daily_mean_airdensity =
#   data.table(
#     matrix(
#       0,
#       nrow =
#         nrow(
#           closest_samp
#         ),
#       ncol =
#         ncol(
#           daily_mean_airdensity
#         )
#     )
#   )
# 
# TMN_daily_mean_airdensity[
#   ,
#   `:=`(
#     V1 =
#       closest_samp[,NCEP_LATITUDE],
#     V2 =
#       closest_samp[,NCEP_LONGITUDE]
#   )
#   ]
# 
# setnames(
#   TMN_daily_mean_airdensity,
#   names(TMN_daily_mean_airdensity),
#   colnames(daily_mean_airdensity)
# )
# setkey(
#   TMN_daily_mean_airdensity,
#   NCEP_LATITUDE,
#   NCEP_LONGITUDE
# )
# 
# TMN_daily_mean_airdensity =
#   daily_mean_airdensity[
#     TMN_daily_mean_airdensity
#     ][
#       ,
#       1:ncol(daily_mean_airdensity),
#       with = FALSE
#       ]
# 
# setkey(
#   closest_samp,
#   TMN_LONGITUDE,
#   TMN_LATITUDE
# )
# setkey(
#   Twenty_mi_nodes_strpd,
#   TMN_LONGITUDE,
#   TMN_LATITUDE
# )
# 
# full_node_coords =
#   closest_samp[Twenty_mi_nodes_strpd]
# 
# setkey(
#   TMN_daily_mean_airdensity,
#   NCEP_LONGITUDE,
#   NCEP_LATITUDE
# )
# setkey(
#   full_node_coords,
#   NCEP_LONGITUDE,
#   NCEP_LATITUDE
# )
# 
# TMN_daily_mean_airdensity[
#   ,
#   `:=`(
#     NCEP_LATITUDE = NULL,
#     NCEP_LONGITUDE = NULL
#   )
#   ]
# 
# TMN_daily_mean_airdensity =
#   cbind(
#     full_node_coords[
#       ,
#       .(
#         Y_INDEX,
#         X_INDEX
#       )
#       ],
#     TMN_daily_mean_airdensity
#   )
# 
# TMN_daily_mean_airdensity =
#   as.matrix(
#     TMN_daily_mean_airdensity
#   )
# 
# Save files --------------------------------------------------------------
# # Air Temperature
# write.table(TMN_daily_mean_temperature,
#             paste(getwd(), '/NorAm_historical_air_temperature_node_info.txt',
#                   sep=''), sep='\t', col.names=T, row.names=F)
# 
# Snow Depth
write.table(
  TMN_daily_mean_snow,
  paste(
    getwd(),
    '/NorAm_historical_snow_node_info.txt',
    sep=''
  ),
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE
)

# # Air Density
# write.table(
#   TMN_daily_mean_airdensity,
#   paste(
#     getwd(),
#     '/NorAm_historical_airdensity_node_info.txt',
#     sep = ''
#   ),
#   sep ='\t',
#   col.names = TRUE,
#   row.names = FALSE
# )
# 
# Create averaged data sets -----------------------------------------------
# # Make room for large files
# rm(list=ls())
# gc()
# 
# # Read in templates
# NODE_DATA = fread("NA_NODE_LC2.txt")
# NODE_DATA = data.frame(NODE_DATA)
# 
# Twenty_mi_nodes = fread("Twenty_mile_node_locations.csv")
# 
# # Read in previously created data frames; this takes ~ 7 mins
# Temp.dt = fread("NA_air_temperature_data.txt")
# Snow.dt = fread("NA_snowdepth_data.txt")
# AirDensity.dt = fread("NA_airdensity_data.txt")
# 
# # Rename columns
# setnames(Temp.dt,names(Temp.dt),c("date","NCEP_LATITUDE","NCEP_LONGITUDE",
#                                   "daily_temp"))
# Temp.dt[,date:=substr(date,1,10)]
# 
# # Rename columns
# setnames(Snow.dt,names(Snow.dt),c("date","NCEP_LATITUDE","NCEP_LONGITUDE",
#                                   "daily_snow"))
# Snow.dt[,date:=substr(date,1,10)]
# 
# # Correct negative snow depth values
# Snow.dt[,daily_snow:=ifelse(daily_snow<0,0,daily_snow)]
# 
# # Rename columns
# setnames(
#   AirDensity.dt,
#   names(AirDensity.dt),
#   c(
#     "date",
#     "NCEP_LATITUDE",
#     "NCEP_LONGITUDE",
#     "daily_airdensity"
#   )
# )
# AirDensity.dt[
#   ,
#   date :=
#     substr(
#       date,
#       1,
#       10
#     )
# ]
# 
# # Air Temperature
# Temp.dt[,date:=substr(date,6,10)]
# MeanTemp.dt=Temp.dt[,round(mean(daily_temp),digits=3),
#                     by=.(date,NCEP_LATITUDE,NCEP_LONGITUDE)]
# setnames(MeanTemp.dt,"V1","daily_temp")
# MeanTemp.dt=MeanTemp.dt[,.(date,NCEP_LATITUDE,NCEP_LONGITUDE,daily_temp)]
# 
# daily_mean_temperature=
#   dcast(MeanTemp.dt,NCEP_LATITUDE+NCEP_LONGITUDE~date,
#         value.var='daily_temp')
# setkey(daily_mean_temperature,NCEP_LATITUDE,NCEP_LONGITUDE)
# 
# # Snow Depth
# Snow.dt[,date:=substr(date,6,10)]
# MeanSnow.dt=Snow.dt[,round(mean(daily_snow),digits=3),
#                     by=.(date,NCEP_LATITUDE,NCEP_LONGITUDE)]
# setnames(MeanSnow.dt,"V1","daily_snow")
# MeanSnow.dt=MeanSnow.dt[,.(date,NCEP_LATITUDE,NCEP_LONGITUDE,daily_snow)]
# 
# daily_mean_snow=
#   dcast(MeanSnow.dt,NCEP_LATITUDE+NCEP_LONGITUDE~date,
#         value.var='daily_snow')
# setkey(daily_mean_snow,NCEP_LATITUDE,NCEP_LONGITUDE)
# 
# # Air Density
# AirDensity.dt[
#   ,
#   date :=
#     substr(
#       date,
#       6,
#       10
#     )
#   ]
# MeanAirDensity.dt =
#   AirDensity.dt[
#     ,
#     round(
#       mean(
#         daily_airdensity
#       ),
#       digits = 3
#     ),
#     by =
#       .(
#         date,
#         NCEP_LATITUDE,
#         NCEP_LONGITUDE
#       )
#     ]
# setnames(
#   MeanAirDensity.dt,
#   "V1",
#   "daily_airdensity"
# )
# MeanAirDensity.dt =
#   MeanAirDensity.dt[
#     ,
#     .(
#       date,
#       NCEP_LATITUDE,
#       NCEP_LONGITUDE,
#       daily_airdensity
#     )
#     ]
# 
# daily_mean_airdensity =
#   dcast(
#     MeanAirDensity.dt,
#     NCEP_LATITUDE + NCEP_LONGITUDE ~ date,
#     value.var = 'daily_airdensity'
#   )
# setkey(
#   daily_mean_airdensity,
#   NCEP_LATITUDE,
#   NCEP_LONGITUDE
# )
# 
# # Define nodes based on distance from node centroid to sampling location
# Twenty_mi_nodes_strpd =
#   Twenty_mi_nodes[
#     ,
#     .(
#       Y_INDEX,
#       X_INDEX,
#       code,
#       LatDegree,
#       LongDegree
#     )
#     ]
# setnames(
#   Twenty_mi_nodes_strpd,
#   old =
#     colnames(Twenty_mi_nodes_strpd),
#   new =
#     c(
#       "Y_INDEX",
#       "X_INDEX",
#       "NODE",
#       "TMN_LATITUDE",
#       "TMN_LONGITUDE"
#     )
# )
# Twenty_mi_nodes_strpd[
#   ,
#   TMN_LONGITUDE :=
#     TMN_LONGITUDE + 360
#   ]
# 
# node_locations =
#   Twenty_mi_nodes_strpd[
#     ,
#     .(
#       TMN_LONGITUDE,
#       TMN_LATITUDE
#     )
#     ]
# node_locations =
#   as.matrix(
#     node_locations,
#     nrow = nrow(Twenty_mi_nodes),
#     ncol = 2
#   )
# 
# #Define nodes based on distance from node centroid to sampling location
# # Air Temperature
# samp_locations = MeanTemp.dt[,.(NCEP_LONGITUDE,NCEP_LATITUDE)]
# samp_locations = unique(samp_locations)
# samp_locations=as.matrix(samp_locations,nrow=nrow(MeanTemp.dt),ncol=2)
# 
# closest_samp=matrix(NA,nrow=nrow(node_locations),ncol=4)
# closest_samp[,1:2]=node_locations
# for(i in 1:nrow(node_locations)){
#   closest_samp[i,3:4]=
#     samp_locations[which.min(distHaversine(node_locations[i,],
#                                            samp_locations)),]
# }
# 
# colnames(closest_samp)=c("TMN_LONGITUDE","TMN_LATITUDE","NCEP_LONGITUDE",
#                          "NCEP_LATITUDE")
# closest_samp=data.table(closest_samp)
# 
# TMN_daily_mean_temperature = 
#   data.table(
#     matrix(
#       0,
#       nrow =
#         nrow(
#           closest_samp
#         ),
#       ncol =
#         ncol(
#           daily_mean_temperature
#         )
#     )
#   )
# 
# TMN_daily_mean_temperature[,`:=`(V1=closest_samp[,NCEP_LATITUDE],
#                                  V2=closest_samp[,NCEP_LONGITUDE])]
# 
# setnames(TMN_daily_mean_temperature,names(TMN_daily_mean_temperature),
#          colnames(daily_mean_temperature))
# setkey(TMN_daily_mean_temperature,NCEP_LATITUDE,NCEP_LONGITUDE)
# 
# TMN_daily_mean_temperature=
#   daily_mean_temperature[TMN_daily_mean_temperature][
#     ,1:ncol(daily_mean_temperature),with=F]
# 
# setkey(closest_samp,TMN_LONGITUDE,TMN_LATITUDE)
# setkey(Twenty_mi_nodes_strpd,TMN_LONGITUDE,TMN_LATITUDE)
# 
# full_node_coords=closest_samp[Twenty_mi_nodes_strpd]
# 
# setkey(TMN_daily_mean_temperature,NCEP_LONGITUDE,NCEP_LATITUDE)
# setkey(full_node_coords,NCEP_LONGITUDE,NCEP_LATITUDE)
# 
# TMN_daily_mean_temperature[,`:=`(NCEP_LATITUDE=NULL,NCEP_LONGITUDE=NULL)]
# 
# TMN_daily_mean_temperature=cbind(full_node_coords[,.(Y_INDEX,X_INDEX)],
#                                  TMN_daily_mean_temperature)
# 
# TMN_daily_mean_temperature=as.matrix(TMN_daily_mean_temperature)
# 
# # Snow Depth
# samp_locations = MeanSnow.dt[,.(NCEP_LONGITUDE,NCEP_LATITUDE)]
# samp_locations = unique(samp_locations)
# samp_locations=as.matrix(samp_locations,nrow=nrow(MeanSnow.dt),ncol=2)
# 
# closest_samp=matrix(NA,nrow=nrow(node_locations),ncol=4)
# closest_samp[,1:2]=node_locations
# for(i in 1:nrow(node_locations)){
#   closest_samp[i,3:4]=
#     samp_locations[which.min(distHaversine(node_locations[i,],
#                                            samp_locations)),]
# }
# 
# colnames(closest_samp)=c("TMN_LONGITUDE","TMN_LATITUDE","NCEP_LONGITUDE",
#                          "NCEP_LATITUDE")
# closest_samp=data.table(closest_samp)
# 
# TMN_daily_mean_snow = 
#   data.table(
#     matrix(
#       0,
#       nrow =
#         nrow(
#           closest_samp
#         ),
#       ncol =
#         ncol(
#           daily_mean_snow
#         )
#     )
#   )
# 
# TMN_daily_mean_snow = data.table(TMN_daily_mean_snow)
# 
# TMN_daily_mean_snow[,`:=`(V1=closest_samp[,NCEP_LATITUDE],
#                           V2=closest_samp[,NCEP_LONGITUDE])]
# 
# setnames(TMN_daily_mean_snow,names(TMN_daily_mean_snow),
#          colnames(daily_mean_snow))
# setkey(TMN_daily_mean_snow,NCEP_LATITUDE,NCEP_LONGITUDE)
# 
# TMN_daily_mean_snow=
#   daily_mean_snow[TMN_daily_mean_snow][
#     ,1:ncol(daily_mean_snow),with=F]
# 
# setkey(closest_samp,TMN_LONGITUDE,TMN_LATITUDE)
# setkey(Twenty_mi_nodes_strpd,TMN_LONGITUDE,TMN_LATITUDE)
# 
# full_node_coords=closest_samp[Twenty_mi_nodes_strpd]
# 
# setkey(TMN_daily_mean_snow,NCEP_LONGITUDE,NCEP_LATITUDE)
# setkey(full_node_coords,NCEP_LONGITUDE,NCEP_LATITUDE)
# 
# TMN_daily_mean_snow[,`:=`(NCEP_LATITUDE=NULL,NCEP_LONGITUDE=NULL)]
# 
# TMN_daily_mean_snow=cbind(full_node_coords[,.(Y_INDEX,X_INDEX)],
#                           TMN_daily_mean_snow)
# 
# TMN_daily_mean_snow=as.matrix(TMN_daily_mean_snow)
# 
# # Air Density
# samp_locations =
#   MeanAirDensity.dt[
#     ,
#     .(
#       NCEP_LONGITUDE,
#       NCEP_LATITUDE
#     )
#     ]
# samp_locations =
#   unique(
#     samp_locations
#   )
# samp_locations =
#   as.matrix(
#     samp_locations,
#     nrow =
#       nrow(
#         MeanAirDensity.dt
#       ),
#     ncol = 2
#   )
# 
# closest_samp =
#   matrix(
#     NA,
#     nrow =
#       nrow(
#         node_locations
#       ),
#     ncol = 4
#   )
# closest_samp[,1:2] =
#   node_locations
# for(i in 1:nrow(node_locations)){
#   closest_samp[i,3:4] =
#     samp_locations[
#       which.min(
#         distHaversine(
#           node_locations[i,],
#           samp_locations
#         )
#       ),
#       ]
# }
# 
# colnames(closest_samp) =
#   c(
#     "TMN_LONGITUDE",
#     "TMN_LATITUDE",
#     "NCEP_LONGITUDE",
#     "NCEP_LATITUDE"
#   )
# closest_samp =
#   data.table(
#     closest_samp
#   )
# 
# TMN_daily_mean_airdensity =
#   data.table(
#     matrix(
#       0,
#       nrow = nrow(closest_samp),
#       ncol = ncol(daily_mean_airdensity)
#     )
#   )
# 
# TMN_daily_mean_airdensity[
#   ,
#   `:=`(
#     V1 =
#       closest_samp[,NCEP_LATITUDE],
#     V2 =
#       closest_samp[,NCEP_LONGITUDE]
#   )
#   ]
# 
# setnames(
#   TMN_daily_mean_airdensity,
#   names(TMN_daily_mean_airdensity),
#   colnames(daily_mean_airdensity)
# )
# setkey(
#   TMN_daily_mean_airdensity,
#   NCEP_LATITUDE,
#   NCEP_LONGITUDE
# )
# 
# TMN_daily_mean_airdensity =
#   daily_mean_airdensity[
#     TMN_daily_mean_airdensity][
#       ,
#       1:ncol(daily_mean_airdensity),
#       with = FALSE
#       ]
# 
# setkey(
#   closest_samp,
#   TMN_LONGITUDE,
#   TMN_LATITUDE
# )
# setkey(
#   Twenty_mi_nodes_strpd,
#   TMN_LONGITUDE,
#   TMN_LATITUDE
# )
# 
# full_node_coords =
#   closest_samp[
#     Twenty_mi_nodes_strpd
#     ]
# 
# setkey(
#   TMN_daily_mean_airdensity,
#   NCEP_LONGITUDE,
#   NCEP_LATITUDE
# )
# setkey(
#   full_node_coords,
#   NCEP_LONGITUDE,
#   NCEP_LATITUDE
# )
# 
# TMN_daily_mean_airdensity[
#   ,
#   `:=`(
#     NCEP_LATITUDE = NULL,
#     NCEP_LONGITUDE = NULL
#   )
#   ]
# 
# TMN_daily_mean_airdensity =
#   cbind(
#     full_node_coords[
#       ,
#       .(
#         Y_INDEX,
#         X_INDEX
#       )
#       ],
#     TMN_daily_mean_airdensity
#   )
# 
# TMN_daily_mean_airdensity =
#   as.matrix(
#     TMN_daily_mean_airdensity
#   )

# Save files --------------------------------------------------------------
# # Air Temperature
# write.table(TMN_daily_mean_temperature,
#             paste(getwd(), "/NorAm_mean_air_temperature_node_info.txt",
#                   sep=""), sep="\t", col.names=T, row.names=F)
# 
# # Snow Depth
# write.table(TMN_daily_mean_snow,
#             paste(getwd(), "/NorAm_mean_snow_node_info.txt",
#                   sep=""), sep="\t", col.names=T, row.names=F)
# 
# # Air Density
# write.table(
#   TMN_daily_mean_airdensity,
#   paste(
#     getwd(),
#     "/NorAm_mean_airdensity_node_info.txt",
#     sep = ""
#   ),
#   sep = "\t",
#   col.names = TRUE,
#   row.names = FALSE
# )

# Read in files for diagnostics for wrapper -------------------------------
{
  # # The model can run with a number of options:
  # #
  # # 1) AVERAGE: 'mean_air_temperature_file' represents the average temperature 
  # # in a given node on a given day, using data from 1957 to 2020. This can be 
  # # used for a static representation of temprature
  # #
  # # 2) HISTORICAL: Any individual year (1957 to 2020) can be used to set the 
  # # node-surface temperatures.
  # #
  # # 3) SIMULATION: Using the mean and the standard deviation of temperature 
  # # per node per day, a simulation can be run using randomly selected 
  # # temperatures from the distribution defined by these two metrics.
  # #
  # # Ultimately, I'd like to put this in a function giving the user the option
  # # to select (1), (2), or (3). This would then call the appropriate script to 
  # # apply the model to averaged data, historical data, or simulation data.
  # # 
  # # # Air Temperature
  # # # AVERAGE
  # # mean_temperature_node_info = 
  # #   fread(paste(getwd(), "/NorAm_mean_air_temperature_node_info.txt", sep=""), 
  # #         header=T)
  # # 
  # # day_of_interest = "12_01"  #must be a character
  # # target_day=
  # #   names(mean_temperature_node_info)[
  # #     which(grepl(day_of_interest,names(mean_temperature_node_info)))]
  # # 
  # # mean_temperature_node_matrix=
  # #   as.matrix(dcast(mean_temperature_node_info,X_INDEX~Y_INDEX,
  # #                   value.var=target_day))
  # # rownames(mean_temperature_node_matrix)=
  # #   mean_temperature_node_matrix[,1]
  # # mean_temperature_node_matrix=mean_temperature_node_matrix[,-1]
  # # 
  # # mat_rot = function(x) apply(t(x), 2, rev)
  # # image.plot(t(mat_rot(mean_temperature_node_matrix)),horizontal=T)
  # 
  # # # HISTORICAL -- just change 'year_of_interest' parameter
  # # historical_temperature_node_data = 
  # #   fread(paste(getwd(),"/NorAm_historical_air_temperature_node_info.txt",sep=""))
  # # 
  # # historical_temperature_node_frame = 
  # #   data.frame(historical_temperature_node_data)
  # # data_columns = colnames(historical_temperature_node_frame)
  # # 
  # # year_of_interest = "2020" #must be a character
  # # day_of_interest = "12_01"  #must be a character
  # # 
  # # data_subset = data_columns[grepl(year_of_interest, data_columns)]
  # # historical_temperature_node_info = 
  # #   historical_temperature_node_frame[,c("Y_INDEX","X_INDEX",data_subset)]
  # # 
  # # target_day=
  # #   names(historical_temperature_node_info)[
  # #     which(grepl(day_of_interest,names(historical_temperature_node_info)))]
  # # 
  # # historical_temperature_node_matrix=
  # #   as.matrix(dcast(historical_temperature_node_info,X_INDEX~Y_INDEX,
  # #                   value.var=target_day))
  # # rownames(historical_temperature_node_matrix)=
  # #   historical_temperature_node_matrix[,1]
  # # historical_temperature_node_matrix=historical_temperature_node_matrix[,-1]
  # # 
  # # mat_rot = function(x) apply(t(x), 2, rev)
  # # image.plot(t(mat_rot(historical_temperature_node_matrix)),horizontal=T)
  # 
  # # # SIMULATION
  # # mean_temperature_node_info = 
  # #   fread(paste(getwd(), "/NorAm_mean_air_temperature_node_info.txt", sep=""), 
  # #         header=T)
  # # 
  # # sd_temperature_node_info = 
  # #   fread(paste(getwd(), "/NorAm_sd_air_temperature_node_info.txt", sep=""), 
  # #         header=T)
  # # 
  # # temperature_node_info = matrix(NA, 25723, 96)
  # # temperature_node_info[,1:4] = cbind(mean_temperature_node_info$Y_INDEX, 
  # #                                     mean_temperature_node_info$X_INDEX, 
  # #                                     mean_temperature_node_info$TMN_LATITUDE,
  # #                                     mean_temperature_node_info$TMN_LONGITUDE)
  # # 
  # # for(i in 5:ncol(temperature_node_info)) {
  # #   temperature_node_info[,i] = 
  # #     rnorm(length(as.data.frame(mean_temperature_node_info)[,i]),
  # #           mean=as.numeric(as.data.frame(mean_temperature_node_info)[,i]), 
  # #           sd=as.numeric(as.data.frame(sd_temperature_node_info)[,i]))  
  # # }
  # # 
  # # sim_temperature_node_info = data.table(temperature_node_info)
  # # setnames(sim_temperature_node_info, old = colnames(sim_temperature_node_info),
  # #          new = colnames(sd_temperature_node_info))
  # 
  # # # Snow Depth
  # # # AVERAGE
  # # mean_snow_node_info = 
  # #   fread(paste(getwd(), "/NorAm_mean_snow_node_info.txt", sep=""), 
  # #         header=T)
  # # 
  # # day_of_interest = "12_01"  #must be a character
  # # target_day=
  # #   names(mean_snow_node_info)[
  # #     which(grepl(day_of_interest,names(mean_snow_node_info)))]
  # # 
  # # mean_snow_node_matrix=
  # #   as.matrix(dcast(mean_snow_node_info,X_INDEX~Y_INDEX,
  # #                   value.var=target_day))
  # # rownames(mean_snow_node_matrix)=
  # #   mean_snow_node_matrix[,1]
  # # mean_snow_node_matrix=mean_snow_node_matrix[,-1]
  # # 
  # # mat_rot = function(x) apply(t(x), 2, rev)
  # # image.plot(t(mat_rot(mean_snow_node_matrix)),horizontal=T)
  #  
  # # # HISTORICAL -- just change 'year_of_interest' parameter
  # # historical_snow_node_data = 
  # #   fread(paste(getwd(),"/NorAm_historical_snow_node_info.txt",sep=""))
  # # 
  # # historical_snow_node_frame = 
  # #   data.frame(historical_snow_node_data)
  # # data_columns = colnames(historical_snow_node_frame)
  # # 
  # # year_of_interest = "1959" #must be a character
  # # day_of_interest = "12_01"  #must be a character
  # # 
  # # data_subset = data_columns[grepl(year_of_interest, data_columns)]
  # # historical_snow_node_info = 
  # #   historical_snow_node_frame[,c("Y_INDEX","X_INDEX",data_subset)]
  # # 
  # # target_day=
  # #   names(historical_snow_node_info)[
  # #     which(grepl(day_of_interest,names(historical_snow_node_info)))]
  # # 
  # # historical_snow_node_matrix=
  # #   as.matrix(dcast(historical_snow_node_info,X_INDEX~Y_INDEX,
  # #                   value.var=target_day))
  # # rownames(historical_snow_node_matrix)=
  # #   historical_snow_node_matrix[,1]
  # # historical_snow_node_matrix=historical_snow_node_matrix[,-1]
  # # 
  # # mat_rot = function(x) apply(t(x), 2, rev)
  # # image.plot(t(mat_rot(historical_snow_node_matrix)),horizontal=T)
  #  
  # # # SIMULATION
  # # mean_snow_node_info = 
  # #   fread(paste(getwd(), "/NorAm_mean_snow_node_info.txt", sep=""), 
  # #         header=T)
  # # 
  # # sd_snow_node_info = 
  # #   fread(paste(getwd(), "/NorAm_sd_snow_node_info.txt", sep=""), 
  # #         header=T)
  # # 
  # # snow_node_info = matrix(NA, 25723, 96)
  # # snow_node_info[,1:4] = cbind(mean_snow_node_info$Y_INDEX, 
  # #                                     mean_snow_node_info$X_INDEX, 
  # #                                     mean_snow_node_info$TMN_LATITUDE,
  # #                                     mean_snow_node_info$TMN_LONGITUDE)
  # # 
  # # for(i in 5:ncol(snow_node_info)) {
  # #   snow_node_info[,i] = 
  # #     rnorm(length(as.data.frame(mean_snow_node_info)[,i]),
  # #           mean=as.numeric(as.data.frame(mean_snow_node_info)[,i]), 
  # #           sd=as.numeric(as.data.frame(sd_snow_node_info)[,i]))  
  # # }
  # # 
  # # sim_snow_node_info = data.table(snow_node_info)
  # # setnames(sim_snow_node_info, old = colnames(sim_snow_node_info),
  # #          new = colnames(sd_snow_node_info))  
}

# Assemble Weather Severity Index files -----------------------------------
# Combine temperature and snow depth into weather table -------------------
Temp.dt =
  fread(
    'NA_air_temperature_data.txt'
    )

setnames(
  Temp.dt,
  names(
    Temp.dt
  ),
  c(
    'date',
    'NCEP_LATITUDE',
    'NCEP_LONGITUDE',
    'daily_temp'
  )
)

Temp.dt[
  ,
  date :=
    substr(
      date,
      1,
      10
    )
]

Temp.dt =
  Temp.dt[
    ,
    round(
      mean(
        daily_temp
      ),
      digits = 3
    ),
    by =
      .(
        date,
        NCEP_LATITUDE,
        NCEP_LONGITUDE
      )
  ]

setnames(
  Temp.dt,
  'V1',
  'daily_temp'
)

Snow.dt =
  fread(
    'NA_snowdepth_data.txt'
  )

setnames(
  Snow.dt,
  names(
    Snow.dt
  ),
  c(
    'date',
    'NCEP_LATITUDE',
    'NCEP_LONGITUDE',
    'daily_snow'
  )
)

Snow.dt[
  ,
  date :=
    substr(
      date,
      1,
      10
    )
]

Snow.dt =
  Snow.dt[
    ,
    round(
      mean(
        daily_snow
      ),
      digits = 3
    ),
    by =
      .(
        date,
        NCEP_LATITUDE,
        NCEP_LONGITUDE
      )
  ]

setnames(
  Snow.dt,
  'V1',
  'daily_snow'
)

Snow.dt[
  ,
  daily_snow :=
    ifelse(
      daily_snow < 0,
      0,
      daily_snow
    )
]

Weather.dt =
  Temp.dt[
    ,
    .(
      date,
      NCEP_LATITUDE,
      NCEP_LONGITUDE
    )
  ]

Weather.dt[
  ,
  `:=`(
    daily_snow_score =
      (
        Snow.dt[
          ,
          (
            (
              daily_snow * 257
            ) /
              400
          )
        ]
      ) *
      0.394,
    daily_temp_score =
      -(
        Temp.dt[
          ,
          daily_temp
        ]
      ),
    year =
      substr(
        date,
        1,
        4
      )
  )
]

# Calculate temp_day and snow_day, from Schummer et al. 2010
Weather.dt[
  ,
  `:=`(
    cold =
      ifelse(
        Temp.dt[
          ,
          daily_temp
        ] < 0,
        1,
        0
      ),
    deep =
      ifelse(
        Snow.dt[
          ,
          daily_snow
        ] >=
          2.54,
        1,
        0
      )
  )
]

Weather.dt[
  ,
  `:=`(
    temp_day =
      sequence(
        rle(
          cold
        )$lengths
      ) * cold,
    snow_day =
      sequence(
        rle(
          deep
        )$lengths
      ) * deep
  ),
  by = year
]

# Calculate WSI 
WSI_data =
  Weather.dt[
    ,
    .(
      daily_temp_score,
      temp_day,
      daily_snow_score,
      snow_day
    )
  ]

WSI_data[
  ,
  `:=`(
    WSI = 
      rowSums(
        WSI_data
      )
  )
]

Weather.dt[
  ,
  `:=`(
    WSI = 
      WSI_data[
        ,
        WSI
      ]
  )
]

Weather.dt = 
  Weather.dt[
    ,
    .(
      date,
      NCEP_LATITUDE,
      NCEP_LONGITUDE,
      daily_temp_score,
      daily_snow_score,
      temp_day,
      snow_day,
      WSI
    )
  ]

WSI.dt = 
  Weather.dt[
    ,
    .(
      date,
      NCEP_LATITUDE,
      NCEP_LONGITUDE,
      WSI
    )
  ]

# Historical WSI per node
HistWSI.dt = WSI.dt

setnames(
  HistWSI.dt,
  'WSI',
  'daily_wsi'
)

daily_mean_wsi =
  dcast(
    HistWSI.dt,
    NCEP_LATITUDE + NCEP_LONGITUDE ~ date,
    value.var = 'daily_wsi'
  )

setkey(
  daily_mean_wsi,
  NCEP_LATITUDE,
  NCEP_LONGITUDE
)

#Define nodes by wsi based on distance from node centroid to
#wsi sampling location
Twenty_mi_nodes_strpd =
  Twenty_mi_nodes[
    ,
    .(
      Y_INDEX,
      X_INDEX,
      code,
      LatDegree,
      LongDegree
    )
  ]

setnames(
  Twenty_mi_nodes_strpd, 
  old =
    colnames(
      Twenty_mi_nodes_strpd
    ),
  new = 
    c(
      'Y_INDEX',
      'X_INDEX',
      'NODE',
      'TMN_LATITUDE',
      'TMN_LONGITUDE'
    )
)

Twenty_mi_nodes_strpd[
  ,
  TMN_LONGITUDE :=
    TMN_LONGITUDE + 360
]

node_locations = 
  Twenty_mi_nodes_strpd[
    ,
    .(
      TMN_LONGITUDE,
      TMN_LATITUDE
    )
  ]

samp_locations = 
  HistWSI.dt[
    ,
    .(
      NCEP_LONGITUDE,
      NCEP_LATITUDE
    )
  ]

samp_locations = 
  unique(
    samp_locations
  )

node_locations =
  as.matrix(
    node_locations,
    nrow =
      nrow(
        Twenty_mi_nodes
      ),
    ncol = 2
  )

samp_locations =
  as.matrix(
    samp_locations,
    nrow =
      nrow(
        HistWSI.dt
      ),
    ncol = 2
  )

closest_samp =
  matrix(
    NA,
    nrow = 
      nrow(
        node_locations
      ),
    ncol = 4
  )

closest_samp[
  ,
  1:2
] =
  node_locations

for(
  i in 1:nrow(
    node_locations
  )
) {
  closest_samp[
    i,
    3:4
  ]=
    samp_locations[
      which.min(
        distHaversine(
          node_locations[
            i,
          ],
          samp_locations
        )
      ),
    ]
}

colnames(
  closest_samp
) = 
  c(
    'TMN_LONGITUDE',
    'TMN_LATITUDE',
    'NCEP_LONGITUDE',
    'NCEP_LATITUDE'
  )

closest_samp =
  data.table(
    closest_samp
  )

TMN_daily_mean_wsi =
  data.table(
    matrix(
      0,
      nrow =
        nrow(
          closest_samp
        ),
      ncol =
        ncol(
          daily_mean_wsi
        )
    )
  )

TMN_daily_mean_wsi[
  ,
  `:=`(
    V1 =
      closest_samp[
        ,
        NCEP_LATITUDE
      ],
    V2 =
      closest_samp[
        ,
        NCEP_LONGITUDE
      ]
  )
]

setnames(
  TMN_daily_mean_wsi,
  names(
    TMN_daily_mean_wsi
  ),
  colnames(
    daily_mean_wsi
  )
)

setkey(
  TMN_daily_mean_wsi,
  NCEP_LATITUDE,
  NCEP_LONGITUDE
)

TMN_daily_mean_wsi=
  daily_mean_wsi[
    TMN_daily_mean_wsi
  ][
    ,
    1:ncol(
      daily_mean_wsi
    ),
    with = FALSE
  ]

setkey(
  closest_samp,
  TMN_LONGITUDE,
  TMN_LATITUDE
)

setkey(
  Twenty_mi_nodes_strpd,
  TMN_LONGITUDE,
  TMN_LATITUDE
)

full_node_coords =
  closest_samp[
    Twenty_mi_nodes_strpd
  ]

setkey(
  TMN_daily_mean_wsi,
  NCEP_LONGITUDE,
  NCEP_LATITUDE
)

setkey(
  full_node_coords,
  NCEP_LONGITUDE,
  NCEP_LATITUDE
)

TMN_daily_mean_wsi[
  ,
  `:=`(
    NCEP_LATITUDE = NULL,
    NCEP_LONGITUDE = NULL
  )
]

TMN_daily_mean_wsi =
  cbind(
    full_node_coords[
      ,
      .(
        Y_INDEX,
        X_INDEX
      )
    ],
    TMN_daily_mean_wsi
  )

# TMN_daily_mean_wsi =
#   as.matrix(
#     TMN_daily_mean_wsi
#   )

fwrite(
  TMN_daily_mean_wsi,
  'NorAm_historical_wsi_node_info.txt',
)

# Average WSI per node
WSI.dt[
  ,
  date :=
    substr(
      date,
      6,
      10
    )
]

MeanWSI.dt =
  WSI.dt[
    ,
    round(
      mean(
        daily_wsi
      ),
      digits = 3
    ),
    by =
      .(
        date,
        NCEP_LATITUDE,
        NCEP_LONGITUDE
      )
  ]

setnames(
  MeanWSI.dt,
  'V1',
  'daily_wsi'
)

MeanWSI.dt =
  MeanWSI.dt[
    ,
    .(
      date,
      NCEP_LATITUDE,
      NCEP_LONGITUDE,
      daily_wsi
    )
  ]

daily_mean_wsi=
  dcast(
    MeanWSI.dt,
    NCEP_LATITUDE + NCEP_LONGITUDE ~ date,
    value.var = 'daily_wsi'
  )

setkey(
  daily_mean_wsi,
  NCEP_LATITUDE,
  NCEP_LONGITUDE
)

#Define nodes by wsi based on distance from node centroid to
#wsi sampling location
setnames(
  Twenty_mi_nodes_strpd, 
  old =
    colnames(
      Twenty_mi_nodes_strpd
    ),
  new = 
    c(
      'Y_INDEX',
      'X_INDEX',
      'NODE',
      'TMN_LATITUDE',
      'TMN_LONGITUDE'
    )
)

node_locations = 
  Twenty_mi_nodes_strpd[
    ,
    .(
      TMN_LONGITUDE,
      TMN_LATITUDE
    )
  ]

samp_locations = 
  MeanWSI.dt[
    ,
    .(
      NCEP_LONGITUDE,
      NCEP_LATITUDE
    )
  ]

samp_locations = 
  unique(
    samp_locations
  )

node_locations =
  as.matrix(
    node_locations,
    nrow =
      nrow(
        Twenty_mi_nodes
      ),
    ncol = 2
  )

samp_locations =
  as.matrix(
    samp_locations,
    nrow = 
      nrow(
        MeanWSI.dt
      ),
    ncol = 2
  )

closest_samp =
  matrix(
    NA,
    nrow = 
      nrow(
        node_locations
      ),
    ncol = 4
  )

closest_samp[
  ,
  1:2
] = 
  node_locations

for(
  i in 1:nrow(
    node_locations
  )
) {
  closest_samp[
    i,
    3:4
  ] =
    samp_locations[
      which.min(
        distHaversine(
          node_locations[
            i,
          ],
          samp_locations
        )
      )
      ,
    ]
}

colnames(
  closest_samp
) =
  c(
    'TMN_LONGITUDE',
    'TMN_LATITUDE',
    'NCEP_LONGITUDE',
    'NCEP_LATITUDE'
  )

closest_samp =
  data.table(
    closest_samp
  )

TMN_daily_mean_wsi =
  data.table(
    matrix(
      0,
      nrow =
        nrow(
          closest_samp
        ),
      ncol =
        ncol(
          daily_mean_wsi
        )
    )
  )

TMN_daily_mean_wsi[,`:=`(V1=closest_samp[,NCEP_LATITUDE],
                          V2=closest_samp[,NCEP_LONGITUDE])]

setnames(TMN_daily_mean_wsi,names(TMN_daily_mean_wsi),
         colnames(daily_mean_wsi))
setkey(TMN_daily_mean_wsi,NCEP_LATITUDE,NCEP_LONGITUDE)

TMN_daily_mean_wsi=
  daily_mean_wsi[TMN_daily_mean_wsi][
    ,1:ncol(daily_mean_wsi),with=F]

setkey(closest_samp,TMN_LONGITUDE,TMN_LATITUDE)
setkey(Twenty_mi_nodes_strpd,TMN_LONGITUDE,TMN_LATITUDE)

full_node_coords=closest_samp[Twenty_mi_nodes_strpd]

setkey(TMN_daily_mean_wsi,NCEP_LONGITUDE,NCEP_LATITUDE)
setkey(full_node_coords,NCEP_LONGITUDE,NCEP_LATITUDE)

TMN_daily_mean_wsi[,`:=`(NCEP_LATITUDE=NULL,NCEP_LONGITUDE=NULL)]

TMN_daily_mean_wsi=cbind(full_node_coords[,.(Y_INDEX,X_INDEX)],
                          TMN_daily_mean_wsi)

TMN_daily_mean_wsi=as.matrix(TMN_daily_mean_wsi)

write.table(TMN_daily_mean_wsi,
            paste(getwd(), '/NorAm_mean_wsi_node_info.txt',
                  sep=''), sep='\t', col.names=T, row.names=F)

# Read in files for diagnostics for wrapper -------------------------------
{
  # The model can run with a number of options:
  #
  # 1) AVERAGE: 'mean_wsi_file' represents the average wsi 
  # in a given node on a given day, using data form 1957 to 2020. This can be 
  # used for a static representation of wsi
  #
  # 2) HISTORICAL: Any individual year (1957 to 2020) can be used to set the 
  # node-surface wsis.
  #
  # 3) SIMULATION: Using the mean and the standard deviation of wsi 
  # per node per day, a simulation can be run using randomly selected 
  # wsis from the distribution defined by these two metrics.
  #
  # Ultimately, I'd like to put this in a function giving the user the option
  # to select (1), (2), or (3). This would then call the appropriate script to 
  # apply the model to averaged data, historical data, or simulation data.
  # 
  # # AVERAGE
  # mean_wsi_node_info = 
  #   fread(paste(getwd(), "/NorAm_mean_wsi_node_info.txt", sep=""), 
  #         header=T)
  # 
  # day_of_interest = "12_01"  #must be a character
  # target_day=
  #   names(mean_wsi_node_info)[
  #     which(grepl(day_of_interest,names(mean_wsi_node_info)))]
  # 
  # mean_wsi_node_matrix=
  #   as.matrix(dcast(mean_wsi_node_info,X_INDEX~Y_INDEX,
  #                   value.var=target_day))
  # rownames(mean_wsi_node_matrix)=
  #   mean_wsi_node_matrix[,1]
  # mean_wsi_node_matrix=mean_wsi_node_matrix[,-1]
  # 
  # mat_rot = function(x) apply(t(x), 2, rev)
  # image.plot(t(mat_rot(mean_wsi_node_matrix)),horizontal=T)
  
  # # HISTORICAL -- just change 'year_of_interest' parameter
  # historical_wsi_node_data = 
  #   fread(paste(getwd(),"/NorAm_historical_wsi_node_info.txt",sep=""))
  # 
  # historical_wsi_node_frame = 
  #   data.frame(historical_wsi_node_data)
  # data_columns = colnames(historical_wsi_node_frame)
  # 
  # year_of_interest = "2011" #must be a character
  # day_of_interest = "12_31"  #must be a character
  # 
  # data_subset = data_columns[grepl(year_of_interest, data_columns)]
  # historical_wsi_node_info = 
  #   historical_wsi_node_frame[,c("Y_INDEX","X_INDEX",data_subset)]
  # 
  # target_day=
  #   names(historical_wsi_node_info)[
  #     which(grepl(day_of_interest,names(historical_wsi_node_info)))]
  # 
  # historical_wsi_node_matrix=
  #   as.matrix(dcast(historical_wsi_node_info,X_INDEX~Y_INDEX,
  #                   value.var=target_day))
  # rownames(historical_wsi_node_matrix)=
  #   historical_wsi_node_matrix[,1]
  # historical_wsi_node_matrix=historical_wsi_node_matrix[,-1]
  # 
  # mat_rot = function(x) apply(t(x), 2, rev)
  # 
  # filename = paste("/NA_Weather_Severity_Index_Charts_",day_of_interest,
  #                  year_of_interest,".jpg",sep="")
  # WSIChartpath = file.path(paste(getwd(),filename,sep=""))
  # jpeg(file = WSIChartpath, width = 1500, height = 1000)
  # # windows(width=10, height=8)
  # # par(bty='n')
  # image.plot(t(mat_rot(historical_wsi_node_matrix)),horizontal=T,
  #            legend.lab=paste(day_of_interest,year_of_interest))
  # dev.off()
  # 
  # filename = paste("/US_Weather_Severity_Index_Charts_",day_of_interest,
  #                  year_of_interest,".jpg",sep="")
  # WSIChartpath = file.path(paste(getwd(),filename,sep=""))
  # jpeg(file = WSIChartpath, width = 1500, height = 1000)
  # # windows(width=10, height=8)
  # # par(bty='n')
  # image.plot(t(mat_rot(historical_wsi_node_matrix))[,1:75],horizontal=T,
  #            legend.lab=paste(day_of_interest,year_of_interest))
  # dev.off()
  
  # # SIMULATION
  # mean_wsi_node_info = 
  #   fread(paste(getwd(), "/NorAm_mean_wsi_node_info.txt", sep=""), 
  #         header=T)
  # 
  # sd_wsi_node_info = 
  #   fread(paste(getwd(), "/NorAm_sd_wsi_node_info.txt", sep=""), 
  #         header=T)
  # 
  # wsi_node_info = matrix(NA, 25723, 96)
  # wsi_node_info[,1:4] = cbind(mean_wsi_node_info$Y_INDEX, 
  #                                     mean_wsi_node_info$X_INDEX, 
  #                                     mean_wsi_node_info$TMN_LATITUDE,
  #                                     mean_wsi_node_info$TMN_LONGITUDE)
  # 
  # for(i in 5:ncol(wsi_node_info)) {
  #   wsi_node_info[,i] = 
  #     rnorm(length(as.data.frame(mean_wsi_node_info)[,i]),
  #           mean=as.numeric(as.data.frame(mean_wsi_node_info)[,i]), 
  #           sd=as.numeric(as.data.frame(sd_wsi_node_info)[,i]))  
  # }
  # 
  # sim_wsi_node_info = data.table(wsi_node_info)
  # setnames(sim_wsi_node_info, old = colnames(sim_wsi_node_info),
  #          new = colnames(sd_wsi_node_info))  
}

# Assemble layers for animation -------------------------------------------
# historical_wsi_node_data =
#   fread(paste(getwd(),"/NorAm_historical_wsi_node_info.txt",sep=""))
# 
# historical_wsi_node_frame = data.frame(historical_wsi_node_data)
# colnames(historical_wsi_node_frame) =
#   gsub("X1","1",colnames(historical_wsi_node_frame))
# colnames(historical_wsi_node_frame) =
#   gsub("X2","2",colnames(historical_wsi_node_frame))
# data_columns = colnames(historical_wsi_node_frame)
# 
# date_range=data_columns[3:length(data_columns)]
# date_range=date_range[!grepl("02_29",date_range)]
# 
# years_of_interest = unique(substr(date_range,1,4))
# days_of_interest = unique(substr(date_range,5,10))
# 
# mat_rot = function(x) apply(t(x), 2, rev)
# 
# pb1 = winProgressBar(title="Years",label="0% done",min=0,max=100,
#                      initial=0)
# for(i in 1:length(years_of_interest)){
#   pb2 = winProgressBar(title="Days",label="0% done",min=0,max=100,
#                        initial=0)
#   for(j in 1:length(days_of_interest)){
#     data_subset = data_columns[grepl(years_of_interest[i], data_columns)]
#     historical_wsi_node_info =
#       historical_wsi_node_frame[,c("Y_INDEX","X_INDEX",data_subset)]
# 
#     target_day=
#       names(historical_wsi_node_info)[
#         which(grepl(days_of_interest[j],names(historical_wsi_node_info)))]
# 
#     historical_wsi_node_matrix=
#       as.matrix(dcast(historical_wsi_node_info,X_INDEX~Y_INDEX,
#                       value.var=target_day))
# 
#     rownames(historical_wsi_node_matrix)=historical_wsi_node_matrix[,1]
# 
#     historical_wsi_node_matrix=historical_wsi_node_matrix[,-1]
# 
#     historical_wsi_node_matrix=mat_rot(historical_wsi_node_matrix)
# 
#     write.table(historical_wsi_node_matrix,
#                 paste(getwd(),"/WSIGrids","/NorAm_comp_wsi_node_info",
#                       paste("_",years_of_interest[i],days_of_interest[j],
#                             sep=""),".txt",sep=""),
#                 sep="\t",col.names=T,row.names=T)
# 
#     info=sprintf("%d%% done",round((j/length(days_of_interest))*100))
#     setWinProgressBar(pb2,j/(length(days_of_interest))*100,label=info)
# 
#   }
# 
#   close(pb2)
#   info=sprintf("%d%% done",round((i/length(years_of_interest))*100))
#   setWinProgressBar(pb1,i/(length(years_of_interest))*100,label=info)
# 
# }
# close(pb1)

# Animate data ------------------------------------------------------------
historical_wsi_node_data = 
  fread(paste(getwd(),"/NorAm_historical_wsi_node_info.txt",sep=""))
historical_wsi_node_frame = data.frame(historical_wsi_node_data)

colnames(historical_wsi_node_frame) = 
  gsub("X1","1",colnames(historical_wsi_node_frame))
colnames(historical_wsi_node_frame) = 
  gsub("X2","2",colnames(historical_wsi_node_frame))
data_columns = colnames(historical_wsi_node_frame)

date_range=data_columns[3:length(data_columns)]
date_range=date_range[!grepl("02_29",date_range)]

# lowest observed WSI = -29.042
# highest observed WSI = 488931.0078
# WSI threshold for mallard departure = 7.2 (Schummer et al. 2010)

brks = c(seq(-30,-7.2,1),1500)
plotcols = c(rev(tim.colors(length(brks)-1)))

wsi_vid = function() {
  for (TS in date_range) {
    options(warn=-1)

    split.screen(rbind(c(0,.8,0,1), c(.8,1,0,1)))
        
    wsi_file=paste0(getwd(),"/WSIGrids","/NorAm_comp_wsi_node_info_",TS,
                    ".txt")
    wsi_grid=fread(wsi_file)
    
    wsi_grid_cols=unlist(strsplit(readLines(wsi_file,1),"\"\t\""))
    wsi_grid_cols[1]=substr(wsi_grid_cols[1],nchar(wsi_grid_cols[1])-1,
                            nchar(wsi_grid_cols[1]))
    wsi_grid_cols[length(wsi_grid_cols)]=
      substr(wsi_grid_cols[length(wsi_grid_cols)],1,
             nchar(wsi_grid_cols[length(wsi_grid_cols)])-1)
    
    colnames(wsi_grid)=c("rows",wsi_grid_cols)
    
    wsi_grid_mat=as.matrix(wsi_grid)
    rownames(wsi_grid_mat)=wsi_grid_mat[,1]
    wsi_grid_mat=wsi_grid_mat[,-1]
    class(wsi_grid_mat)="numeric"
    
#     image(t(wsi_grid_mat),axes=F,breaks=brks,col=plotcols)
#     title(main=paste0(format(as.Date(TS,"%Y_%m_%d"),format="%d %B %Y")))
    
    screen(1)
    image(t(wsi_grid_mat),axes=F,breaks=brks,col=plotcols)
    title(main=paste0(format(as.Date(TS,"%Y_%m_%d"),format="%d %B %Y")))
    
    screen(2)
    image.plot(zlim=c(-30,10),legend.only=TRUE,col=plotcols,
               legend.lab="Cumulative WSI (Schummer et al. 2010)",
               smallplot=c(.1,.2,.3,.7))
    close.screen(all.screens=T)
  }
}

saveVideo(wsi_vid(),interval=0.1,video.name="Daily_NA_WSI_extendedcut.mp4")


