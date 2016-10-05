#######################################################################
#
# Main Migration Wrapper
#---------------------------------------------------------------------------------------
#
# Created by Sarah Jacobi, Eric Lonsdorf, Kevin Aagaard, Wayne Thogmartin
# Tim Jones and Vicky Hunt in collaboration with the Integrated Waterbird Management 
# and Monitoring Program
#
# Modified:  
Sys.time()
# By: Sarah Jacobi - adding new speed ups and cleaning up code. Trying to push.
#
# Create necessary inputs used in the movement script
######################################################################

# This script is used to gather the necessary data for the flyway model.
# Vicky will now work on this version to implement data table improvements.
# Figured out how to get pull and push buttons back.

#Clear the workspace
rm(list=ls())

#Set working directory
# #For Sarah:
# setwd("~/IDriveSync/IWMMGit/inputs")

#For Kevin
.libPaths("C:/R-3.3.1/library")
work_dir = file.path("//Igsarfebfslacr3","Users","kaagaard","My Documents",
                    "IWMM Efforts","Flyway Model Project","Modelling",
                    "Scripts and data")
setwd(work_dir)
setwd(paste(getwd(),'/R',sep=''))

#Include necessary libraries
library(compiler)
library(fields)
library(maptools)
library(data.table)
library(plyr)
library(foreign)
library(dplyr)

#Enable the just in time compiler to speed up analysis
enableJIT(3)

#Grab current time
timeit<-proc.time()

#Set the significant digits to 15
options(digits=16)

# First, read in the node data from NA_NODE_LC2.csv
NODE_DATA <- fread("NA_NODE_LC2.txt")
NODE_DATA = data.frame(NODE_DATA)

# These are the columns in NODE_DATA:
# 1:FID; 2:Y-index; 3:X-index; 4:code; 5:long; 6:lat; 7:shoreline; 8: open
# water; 9: ice/snow; 10:open space 11:low-int devel; 12:med-intense
# devel; 13:high-int develop; 14:barren; 15:Decid Forest; 16:Evergreen
# forest; 17:mixed forest; 18:dwarf scrub; 19:Shrub/scrub 20:herb grass;
# 21:sedge/herb; 22:moss; 23:pasture/hay; 24:crops; 25:woody wetlands;
# 26:herb wetlands; 27:mean dist to crops; 28:std to crops; 29:mean woody
# wet; 30:std wood wet; 31:mean herb wet; 32:std herb wet

# Determine number of nodes
num_nodes <- nrow(NODE_DATA)

#Rename columns to make it easier to reference them
setnames(NODE_DATA, old= c("VALUE_1", "VALUE_82" ,"VALUE_90","VALUE_95", "MEAN_82", "MEAN_90", "MEAN_95", "STD_82", "STD_90" ,"STD_95"),
         new=c("Shoreline", "Crops", "Woody_wet","Herb_wet", "Mean_crops","Mean_woody_wet","Mean_herb_wet", "STD_Crops","STD_woody_wet", "STD_herb_wet"))


#Set node size in miles
node_size = 20 

# set other variables
node_move <- 3000 #How far birds can move within nodes in meters
num_days <<- 237 #number of days in migration
num_land_cover_class= 3 #number of land cover classes currently analyzed

###############################################################################
# This section is used to determine the distribution of starting birds and 
# forage/roosting data
###############################################################################

#set up matrix to hold habitat data
hab_dist <- matrix(0,num_nodes,num_land_cover_class)

# Incorporate distance effects. Create normal distribution using mean and stdev dist to cover
hab_dist[,1]<-pnorm(node_move,NODE_DATA[,"Mean_crops"],sapply(NODE_DATA[,"STD_Crops"], function(x) max(x,30)))
hab_dist[,2]<-pnorm(node_move,NODE_DATA[,"Mean_woody_wet"],sapply(NODE_DATA[,"STD_woody_wet"], function(x) max(x,30)))
hab_dist[,3]<-pnorm(node_move,NODE_DATA[,"Mean_herb_wet"],sapply(NODE_DATA[,"STD_herb_wet"], function(x) max(x,30)))
hab_dist[which(NODE_DATA[,"Crops"]==0),1]<-0
hab_dist[which(NODE_DATA[,"Woody_wet"]==0),2]<-0
hab_dist[which(NODE_DATA[,"Herb_wet"]==0),3]<-0

#Grab Y_INDEX and X_INDEX
coords <- NODE_DATA[,2:3]

#Grab FID, Y_INDEX, X_INDEX, Shoreline, crops, woody wet and herb wet. Mult by cum. norm value
bird_hab <- cbind(NODE_DATA[1:6],as.numeric(NODE_DATA[,"Shoreline"]),NODE_DATA[,names(NODE_DATA) %in% c("Crops","Woody_wet","Herb_wet")]*hab_dist) 

#Rename columns
setnames(bird_hab,c("as.numeric(NODE_DATA[, \"Shoreline\"])"),c("Shoreline"))

rm(hab_dist) #remove unneeded data
bird_hab[,7:10] <- bird_hab[,7:10]/900 #Data is in square meters.  Covert to pixels
#convert area into pixels and each pixel is a 30x30 cell (900 meters squared)

#This calculation is used to distribute the birds among breeding nodes (shoreline + herb grass)
#33945 kcal = 143k kJ, based on energy in breeding habitat

cals_by_cover <- bird_hab[,7:10] #Use this to keep land cover types separate going forward

#fixed on August 8, 2016
NODE_DATA$cals_sum<-rowSums(cals_by_cover)*33945 #Use this total number of kcal to distribute the birds

# Load data on mallard distribution according to BPOP and NatureServe
# We should really just change the root file name to remove the extra 'e'.
# Also, the double periods before the 'dbf' throw and error for me. Also
# probably a file name issue.
S<- read.dbf("NorthAmerica_20mi_grid_wAK_BPOP_NSmallard_join.dbf") #much faster to read in

#Put Shape data into temp variable.  Note that Natureserve = 3 is non-breeding only; 2 is breeding only
tempS <- data.frame(S$ID, S$Y_INDEX, S$X_INDEX, S$NORMPOP, S$COUNT, S$Natureserv)

#Grab number of columns in NODE_DATA
c<-ncol(NODE_DATA)

#Grab the NORMPOP for only breeding nodes
atest <- S$NORMPOP*(S$Natureserv==2)
atest[is.na(atest)]<-0 #Set NA to zero

bird_start <- 1000000*sum(atest) #Sum up NORMPOP and mult. by 1M to get 20M birds to start

tempS$conc <- paste0(tempS[,2], tempS[,3]) #Grab x,y coords from naturserve data
NODE_DATA$conc <- paste0(NODE_DATA[,2],NODE_DATA[,3]) #Grab x,y coords from NOde data
BB1 <- tempS[tempS[,"conc"] %in% NODE_DATA[,"conc"],1] #Find overlapping node coords
BB <- BB1[1:num_nodes] #unsure why we get rid of last node

#Node definition - use nature serve range maps to select nodes
NODE_DATA$newcol1 <- tempS[BB,4] #Normpop
NODE_DATA$newcol2 <- tempS[BB,6] #Naturserv category

#Grab number of nodes 
num_bird_nodes <- nrow(NODE_DATA)

#Identify breeding and wintering nodes

#Subset NODE_DATA by breeding nodes
NODE_DATA_breed<-data.table(NODE_DATA)

#Calculate normalized BPOP normpop over breeding nodes
NODE_DATA_breed=NODE_DATA_breed[which(NODE_DATA_breed[,newcol2]==2),BPOP_wt:=newcol1/sum(newcol1,na.rm=TRUE)]

#Calculate normalized calories over breeding nodes
NODE_DATA_breed = NODE_DATA_breed[which(NODE_DATA_breed[,newcol2]==2),cal_wt:=cals_sum/sum(cals_sum,na.rm=TRUE)]

#Use the normalized calories and BPOPs to determine normalized wt for breeding nodes
NODE_DATA_breed = NODE_DATA_breed[which(NODE_DATA_breed[,newcol2]==2),breed_wt:=BPOP_wt*cal_wt/sum(BPOP_wt*cal_wt,na.rm=TRUE)]

#Distribute the birds based on the breeding wt
NODE_DATA_breed = NODE_DATA_breed[which(NODE_DATA_breed[,newcol2]==2),breed_pop:=round(bird_start*breed_wt)]


#Clean up temp variables
rm(S,tempS,BB,BB1,atest)
gc()

#Grab number of new subset of rows
num_nodes <<- nrow(NODE_DATA)

#Grab node id, long, lat
mtxNode <- data.frame(NODE_DATA[,1],NODE_DATA[,5],NODE_DATA[,6])

#subset of mtxNode that are breeding nodes
mtxBreedNodes <- NODE_DATA_breed[which(NODE_DATA_breed[,newcol2]==2),.(FID,Long,Lat)]

#start_nodes <- breeding_nodes
init_birds <- NODE_DATA_breed[which(NODE_DATA_breed[,newcol2]==2),breed_pop]

#Determine number of start_nodes
n_start <- length(init_birds)

#Determine proportion of starting birds in each starting node
start_birds <- init_birds/sum(init_birds)

#We assume that max flight distance is fixed at 1400 miles
#Nodes are 20 miles wide x 20 miles tall

#There are two ways to determine window radius.  1) Fix the max flight distance. 2) grab distance from parameter file.
Max_flight_distance= 1400 #miles
#flight_dist<- parameters[rr,10] #Make sure this is the same value as is used to create window size
#flight_dist = 1400/.6 #Convert to km
flight_dist = Max_flight_distance
window_radius = Max_flight_distance/node_size

#Create a moving window matrix
moving_window = matrix(0,window_radius*2+1,window_radius*2+1)
rownames(moving_window)<-seq(-window_radius,window_radius,by=1)
colnames(moving_window)<-seq(-window_radius,window_radius,by=1)

#Fill in moving_window with the distance between the central node to each radial node
for (i in 0:window_radius) {
  for (j in 0:window_radius) {
    rowstring = toString(i)
    colstring = toString(j)
    negrowstring=toString(-i)
    negcolstring = toString(-j)
    
    if (rowstring=="0") {
      moving_window[rowstring,colstring] = node_size * j
      moving_window[rowstring,negcolstring] = node_size * j
    }    
    else if (colstring == "0") {
      moving_window[rowstring,colstring] = node_size * i
      moving_window[negrowstring,colstring] = node_size * i
    }
    else {
      moving_window[rowstring,colstring] = sqrt((node_size * i)^2 + (node_size * j)^2)
      moving_window[negrowstring,colstring] = sqrt((node_size * i)^2 + (node_size * j)^2)
      moving_window[rowstring,negcolstring] = sqrt((node_size * i)^2 + (node_size * j)^2)
      moving_window[negrowstring,negcolstring] = sqrt((node_size * i)^2 + (node_size * j)^2)
    }
  }
}

moving_window[moving_window>Max_flight_distance]<-0 #Zero out any nodes beyond max flight distance
moving_window=moving_window*1.609344 #convert to km

#Calculate the distance from each node to each breeding node
distance_to_breed_nodes <- rdist(mtxNode[,2:3],mtxBreedNodes[,2:3,with=FALSE])


#Remove unneeded variables
rm(mtxNode,mtxBreedNodes)

#Determine distance to nearest breeding node
distance_to_breed <-as.matrix(apply(distance_to_breed_nodes,1,min))

#Garbage collector
gc()

#Sum up all the habitat data contained in columns 7 to 26 of NODE_DATA
tot_hab <- rowSums(NODE_DATA[,7:26])

#Get percent of total habitat that is shoreline (column 7)
percent_shore <- as.matrix(NODE_DATA[,"Shoreline"]/tot_hab) # dabbler roost

#Grab node id, y_index and x_index
mtxNode <- NODE_DATA[,1:3]

#Disturbance will affect the movement/expenditures to forage
#May affect leave/stay decisions
#Disturbance will decrease daily gain, maybe enough to force them to leave

#Create different distrubance weights for each cover type (open space, low intense, med intense, 
#high intense devel.)
disturb_weights <- c(0.4, 0.6, 0.8, 1)

#Disturbance is calculated as sum of open space, low/med/high intensity development divided by 
#total habitat weight by disturbance weights
disturbance <-as.matrix(NODE_DATA[,10:13])%*%(disturb_weights)/tot_hab

#Daily gain was found by 100*(1-sum(disturbance))
#This disturbance only includes "human non-hunting" effects.  Need to add hunting and temp
#inpacts on daily gain
#DG = (Beta0i*di*ingestion)-BMR
#Beta0i = current measure of human disturbance
#di= disturbance coefficient describing harvest
#ingestion=describe kJ consumed per day (maybe per habitat perhaps)
#BMR= basal metabolic rate that is temp dependent

#current DG is in kcals I believe. body condition is in grams.  Each body condition bin differs by 12.5 grams.
#For now, place holder moves birds one body condition bin per day.  This is the same as converting 100 kcals to 12.5 grams
DG <<- 100*rep(1,num_nodes)*(1-disturbance)/8 #8 is the conversion factor to change 100 kcal to 12.5 grams


# Parameter scenarios
# c1 param ID; c2 shoreline cals; c3 crop cals; c4 herb wet
# cals; c5 woody wet cals; c6 tank cals; c7 flight cost; c8 flight speed;
# c9 flight fuel efficiency (c8/c7); c10 flight range (c9*c6)

param_scenarios <- read.delim("param_scenarios_5_8_14.txt",header= TRUE)
gc()

rr=1 #This is the parameter row

cals=matrix(0,length(bird_hab),4) #Create empty matrix to hold calorie data for 
#shoreline, crops, woody wet, herb wet

#Keep cals separate by cover type
cals.shore = as.matrix(as.numeric(bird_hab[,"Shoreline"])*as.matrix(param_scenarios[rr,"shoreline"]))
cals.crops = as.matrix(as.numeric(bird_hab[,"Crops"])*as.matrix(param_scenarios[rr,"crops"]))
cals.woody_wet = as.matrix(as.numeric(bird_hab[,"Woody_wet"])*as.matrix(param_scenarios[rr,"woody_wet"]))
cals.herb_wet = as.matrix(as.numeric(bird_hab[,"Herb_wet"])*as.matrix(param_scenarios[rr,"herbaceous_wet"]))

#Add node id data to input data
forage=cbind(NODE_DATA[,1:6],cals.shore,cals.crops,cals.woody_wet,cals.herb_wet)
roosting = cbind(NODE_DATA[,1:6],percent_shore)
distance_to_breed = cbind(NODE_DATA[,1:6],distance_to_breed)
DG = cbind(NODE_DATA[,1:6],DG)

###################################################################################
# Calculate distances within moving window using gamma distribution
###################################################################################

range_limit <- 0.8 # can travel 80% of total
gamma_a <- 1 # used to calculate distance prob.

parameters= param_scenarios[rr,] #Grab parameter data 

# Calculate distance birds can fly based on starting with a full tank

# Determine probability of going from central node to each node based on gamma
# distribution
gamma_x = ((flight_dist-moving_window)*(moving_window>0))
gamma_b = flight_dist-flight_dist*range_limit

t1 <- gamma_x^(gamma_a-1)
t2 <- gamma_x/gamma_b
t3 <- gamma(gamma_a)

distance_to_arrival<- ((t1)*(exp(1)^(-t2)/(gamma_b^gamma_a*t3)))
distance_to_arrival<-distance_to_arrival/sum(sum(distance_to_arrival))

#create matrix of grid nodes in moving window
xgridmin = min(NODE_DATA[,"X_INDEX"])
xgridmax = max(NODE_DATA[,"X_INDEX"])

ygridmin = min(NODE_DATA[,"Y_INDEX"])
ygridmax= max(NODE_DATA[,"Y_INDEX"])

grid_nodes = matrix (0, xgridmax-xgridmin+1, ygridmax- ygridmin+1)

rownames(grid_nodes) = seq(xgridmin,xgridmax,by=1)
colnames(grid_nodes) = seq(ygridmin, ygridmax, by=1)

#Convert the ordering of the Y_INDEX field to go from small to large
forage=data.table(forage)
forage[,`:=`(New_Y_INDEX=max(Y_INDEX)+1-Y_INDEX)]
roosting=data.table(roosting)
roosting[,`:=`(New_Y_INDEX=max(Y_INDEX)+1-Y_INDEX)]
distance_to_breed=data.table(distance_to_breed)
distance_to_breed[,`:=`(New_Y_INDEX=max(Y_INDEX)+1-Y_INDEX)]
DG = data.table(DG)
DG[,':='(New_Y_INDEX=max(Y_INDEX)+1-Y_INDEX)]

#put data into correct format to match with moving_window
NODE_DATA_breed[,':='(NEW_Y_INDEX=max(Y_INDEX)+1-Y_INDEX)]
init_birds.matrix= as.matrix(dcast(NODE_DATA_breed,X_INDEX~NEW_Y_INDEX, value.var="breed_pop"))
rownames(init_birds.matrix) = init_birds.matrix[,1]
init_birds.matrix = init_birds.matrix[,-1]

forage.shore.matrix = as.matrix(dcast(forage, X_INDEX~New_Y_INDEX,value.var="cals.shore"))
forage.crops.matrix = as.matrix(dcast(forage, X_INDEX~New_Y_INDEX,value.var="cals.crops"))
forage.woody_wet.matrix = as.matrix(dcast(forage, X_INDEX~New_Y_INDEX,value.var="cals.woody_wet"))
forage.herb_wet.matrix = as.matrix(dcast(forage, X_INDEX~New_Y_INDEX,value.var="cals.herb_wet"))
rownames(forage.shore.matrix)=forage.shore.matrix[,1]
rownames(forage.crops.matrix)=forage.crops.matrix[,1]
rownames(forage.woody_wet.matrix)=forage.woody_wet.matrix[,1]
rownames(forage.herb_wet.matrix)=forage.herb_wet.matrix[,1]
forage.shore.matrix=forage.shore.matrix[,-1]
forage.crops.matrix=forage.crops.matrix[,-1]
forage.woody_wet.matrix=forage.woody_wet.matrix[,-1]
forage.herb_wet.matrix=forage.herb_wet.matrix[,-1]

roosting.matrix=as.matrix(dcast(roosting,X_INDEX~New_Y_INDEX,value.var="percent_shore"))
rownames(roosting.matrix)=roosting.matrix[,1] #Set row names
roosting.matrix=roosting.matrix[,-1] #remove column of row names

distance_to_breed.matrix=as.matrix(dcast(distance_to_breed,X_INDEX~New_Y_INDEX,value.var="distance_to_breed"))
rownames(distance_to_breed.matrix)=distance_to_breed.matrix[,1]
distance_to_breed.matrix=distance_to_breed.matrix[,-1]

DG.matrix = as.matrix(dcast(DG,X_INDEX~New_Y_INDEX,value.var="DG"))
rownames(DG.matrix)=DG.matrix[,1]
DG.matrix = DG.matrix[,-1]

#Plot images of input layers
image.plot(forage.shore.matrix)
image.plot(forage.crops.matrix)
image.plot(forage.woody_wet.matrix)
image.plot(forage.herb_wet.matrix)
image.plot(roosting.matrix)
image.plot(distance_to_breed.matrix)
image.plot(distance_to_arrival)
image.plot(DG.matrix)

#Augment the matrices by adding zeros the size of the radius in all directions
min_data_node = window_radius +1 #This is the node id that starts holding data for the actual landscape.  Holds for x and y

num_rows_orig = nrow(forage.shore.matrix) #size of unpadded landscape
num_cols_orig= ncol(forage.shore.matrix)

num_rows = nrow(forage.shore.matrix)+2*window_radius #add padding
num_cols = ncol(forage.shore.matrix)+2*window_radius
min_x_node = window_radius+1
min_y_node = window_radius+1
max_x_node = window_radius+nrow(forage.shore.matrix) 
max_y_node=window_radius+ncol(forage.shore.matrix) 

#add NA's of padding for easier moving window analysis
forage.shore.matrix.full = matrix(data=NA,nrow=num_rows, ncol =num_cols)
forage.crops.matrix.full = matrix(data=NA,nrow=num_rows, ncol =num_cols)
forage.woody_wet.matrix.full = matrix(data=NA,nrow=num_rows, ncol =num_cols)
forage.herb_wet.matrix.full = matrix(data=NA,nrow=num_rows, ncol =num_cols)
roosting.matrix.full = matrix(data=NA,nrow=num_rows, ncol =num_cols)
distance_to_breed.matrix.full = matrix(data=NA,nrow=num_rows, ncol =num_cols)
DG.matrix.full = matrix(data=NA,nrow=num_rows, ncol =num_cols)

forage.shore.matrix.full[min_x_node:max_x_node, min_y_node:max_y_node] = forage.shore.matrix
forage.crops.matrix.full[min_x_node:max_x_node, min_y_node:max_y_node] = forage.crops.matrix
forage.woody_wet.matrix.full[min_x_node:max_x_node, min_y_node:max_y_node] = forage.woody_wet.matrix
forage.herb_wet.matrix.full[min_x_node:max_x_node, min_y_node:max_y_node] = forage.herb_wet.matrix
roosting.matrix.full[min_x_node:max_x_node, min_y_node:max_y_node] = roosting.matrix
distance_to_breed.matrix.full[min_x_node:max_x_node, min_y_node:max_y_node] = distance_to_breed.matrix
DG.matrix.full[min_x_node:max_x_node, min_y_node:max_y_node] = DG.matrix

#Incorporate decay by land cover - this can be preprocessed
forage.by.day = array(data=NA, dim=c(num_rows,num_cols,num_days))

#Place holder for decay rate
decay = as.matrix(c(0.9998, 0.997, 0.9965, 0.991),nrows = 4, ncol =1)
rownames(decay)=c("shoreline","crops","woody_wet","herb_wet")

for (i in 1:num_days) {
  forage.by.day[,,i]=forage.shore.matrix.full*decay["shoreline",1]^(i-1)+forage.crops.matrix.full*decay["crops",1]^(i-1)+
    forage.woody_wet.matrix.full*decay["woody_wet",1]^(i-1)+forage.herb_wet.matrix.full*decay["herb_wet",1]^(i-1)
}


#Create body condition bins. Here we assume 20
body_condition_bins = 20

#Include daily survival for each body condition class

#Bottom bin has 99.1% daily survivorship
min_survival = .991
#Top bin has 99.997% daily survivorship
max_survival = .99997

daily_survival = matrix(0,body_condition_bins)

#Divide range equally
for (i in 1:body_condition_bins) {
  daily_survival[i,1] = (max_survival - min_survival)/body_condition_bins*i + min_survival
}

#Assume that the energy used in migration is from lipid reserves exclusively (Whyte and Bolen 1988)
#Note: this assumption is not entirely true.  Hanson et al. (1990) show carbs used in Sept. and lipids Nov/Dec.

#Divide expected max lipid levels into 20 bins (see Whyte and Bolen 1988 for similar construct)
#Mallards have low end of 20g and high end of 270 g.

#Create a vector of body conditions from the lowest fat reserves to highest

body_low = 20 #grams
body_high = 270 #grams

body_condition_bin_size = (body_high-body_low)/body_condition_bins
body_conditions_low = seq(body_low,body_high-body_condition_bin_size,by=body_condition_bin_size)
body_conditions_high = seq(body_low+body_condition_bin_size,body_high,by=body_condition_bin_size)
body_conditions = cbind(body_conditions_low, body_conditions_high)

#add mean
body_conditions=cbind(body_conditions,rowMeans(body_conditions))

#For now, create two 4-D arrays.  One for am and one for pm
#Create a 4-D array to use for birds. Dim = days, x by y nodes (including padding), body condition bin

bird_pop.am = array(data=NA, dim=c(num_days,num_rows-2*window_radius,num_cols-2*window_radius,body_condition_bins))
bird_pop.pm = bird_pop.am

#Place all the initial birds into the top body condition on start of day 1
bird_pop.am[1,,,body_condition_bins] = init_birds.matrix
bird_pop.am[is.na(bird_pop.am)] =0


#Create a logical array that descirbes the possibility of moving from each body condition to each other
#for all combos of nodes in the window.  This is not temp dependent

body_change_logical = array(data=NA,dim=c(body_condition_bins,body_condition_bins,dim(distance_to_arrival)))

#Empirically determined cost of flight per unit distance = HBMR*12/Flight speed.
#For this, we assume HBMR in thermal-neutral zone = 18.8KJ and flight speed = 72 km/hr
#Convert energy in KJ to g of lipids: 37 kJ/g
convert_to_g = (18.8*12/72)/37 #units = g/km
moving_window.grams = moving_window*convert_to_g#This stores the cost in grams of moving from focal node to all nodes in window

#Create "pages" of logical node maps describing which from-to body conditions are possible
for (i in 1:body_condition_bins){
  new_body = (body_conditions[i,3]-moving_window.grams)*(moving_window.grams!=0)
  new_body = new_body*(new_body>0) #remove any unattainable nodes
  for (j in 1:body_condition_bins) {
    body_change_logical[i,j,,]=(new_body>=body_conditions[j,1])&(new_body<body_conditions[j,2])  
  }
}


######################################################################################################
#
#Starting from the first node, create a window around it. Normalize forage, roosting, distance based on
#WSI mask.
#
#
######################################################################################################

#read in WSI data from text file
WSI_DATA <- fread("NorAm_historical_wsi_node_info.txt")
WSI_DATA_frame = data.frame(WSI_DATA)
data_columns = colnames(WSI_DATA_frame)

#subset out a particular migration period (don't want to base this on calendar year, but
#rather on annual migratory cycle; i.e., fall one year to spring the next)
target_year = 2002
nonbreeding_year = c(as.character(target_year),as.character(target_year+1))
target_dates=c(paste0(nonbreeding_year[1],c("_09_","_10_","_11_","_12_")),
               paste0(nonbreeding_year[2],c("_01_","_02_","_03_","_04_")))

nonbreeding_period = data_columns[grepl(target_dates[1],data_columns)]
for(i in 2:length(target_dates)){
  nonbreeding_period = 
    c(nonbreeding_period,data_columns[grepl(target_dates[i],data_columns)]) 
}

WSI_DATA.nonbreeding_period = WSI_DATA_frame[,c('X_INDEX','Y_INDEX',nonbreeding_period)]
#Convert to moving_window format
WSI=data.table(WSI_DATA.nonbreeding_period)
WSI[,`:=`(New_Y_INDEX=max(Y_INDEX)+1-Y_INDEX)]

#Create empty matrix with padding to store WSI data
WSI.array = array(data=NA,dim=c(dim(roosting.matrix.full),ncol(WSI-2)))

#Create "pages" by day of migration 
#create counter to switch from column/date names to days
j=1
for (i in colnames(WSI)) { #Loop through the columns, each is a day
  if (grepl("INDEX",i)==FALSE){ #First 2 columns are location data
    A = as.matrix(dcast(WSI,X_INDEX~New_Y_INDEX,value.var=i)) #Covert to moving_window "raster"
    rownames(A)=A[,1] #Move row names
    A=A[,-1] #Remove row names from data
    WSI.array[min_x_node:max_x_node, min_y_node:max_y_node,j] = A # fill in the non-padding portion of array with data
    j=j+1
    
  }
}

#Place holder for leave/stay rule
prob_depart.WSI = (WSI.array+10)^3/((WSI.array+10)^3+10^3)*(WSI.array>0)

#prob_depart.DG = linear decreasing function. In Chicago meeting in April 2016, 1-(DG+1)/5. DG = {-1,0,1,2,3,
weight.WSI = 0.5
weight.BC = 0.5

#Place holder for weights
W1 = .25
W2 = .25
W3 = .25
W4 = .25

#Limit loops by only grabbing nodes with data
data_nodes = which(!is.na(roosting.matrix.full),TRUE)

#remove items not needed to reduce memory

rm(list= ls()[!(ls() %in% c('init_birds.matrix','num_rows_orig', 'num_cols_orig','num_days','roosting.matrix','data_nodes','window_radius','bird_pop.am','forage.by.day','roosting.matrix.full',
                            'distance_to_breed.matrix.full', 'WSI.array', 'distance_to_arrival','W1', 'W2', 'W3', 'W4','body_condition_bins',
                            'DG.matrix.full', 'weight.WSI', 'prob_depart.WSI', 'daily_survival', 'weight.BC', 'bird_pop.pm','body_conditions','body_change_logical'))])
gc()

#Pre-determine the windows around each possible focal node.

#Grad x any y indices of any node with data (not NA)
focal_node_x =data_nodes[,1]
focal_node_y=data_nodes[,2]

#Used to convert to correct indices based on padding of zeros
x = focal_node_x-window_radius
y=focal_node_y-window_radius

#Determine lower and upper bound on window
lower_bound_x = focal_node_x-window_radius 
upper_bound_x = focal_node_x+window_radius

lower_bound_y = focal_node_y-window_radius
upper_bound_y = focal_node_y+window_radius

#These are lower and upper bounds used to determine how the smaller array that stores number of birds
#maps to the full padded arrays that store data

lb_x = ifelse(x-window_radius<1,1,x-window_radius)
up_x = ifelse(x+window_radius>num_rows_orig, num_rows_orig, x+window_radius)

lb_y = ifelse(y-window_radius<1,1,y-window_radius)
up_y = ifelse(y+window_radius>num_cols_orig, num_cols_orig, y+window_radius)

#used to convert data indices to bird pop indices since they aren't the same size
lb_x_small = ifelse((x-window_radius)>=1,1,((window_radius*2+1)-(up_x-lb_x)))
up_x_small = ifelse((x+window_radius)<=num_rows_orig, (window_radius*2+1), (up_x-lb_x+1))

lb_y_small = ifelse((y-window_radius)>=1, 1, (window_radius*2+1)-(up_y-lb_y))
up_y_small = ifelse((y+window_radius)<=num_cols_orig,(window_radius*2+1), (up_y-lb_y+1))

WSI.cutoff = 7.5
#######################################################################
#Define necessary functions used in for loop

#######################################################################
#Function Name: new_body_condition
#Purpose: determines the resulting body condition after foraging
#
#Inputs: body_conditions (matrix holding body condition data)
#        new.bc (the resulting body_condition after daily gain)
#
#Outputs: new.bc.bin (discrete bin in which resulting body condition falls)
#        
######################################################################
new_body_condition <- function(body_conditions, bcvector) {
  if (new.bc > max(body_conditions[,2])) {
    new.bc.bin = nrow(body_conditions)
  } else if (new.bc < min(body_conditions[1])) {
    new.bc.bin = 0 #don't think this can happen
  } else if (new.bc == min(body_conditions[1])) {
    new.bc.bin = 1
  } else {
    new.bc.bin = which(new.bc>body_conditions[,1] & body_conditions[,2]>=new.bc)
  }
  return(new.bc.bin)
} 

new_bc_vector <- function(body_conditions, bcvector) {
  which(bcvector>max(body_conditions[,2]))
}

getNewBC <- function(bcvector, body_conditions) {
  tmp <- body_conditions %>%
    filter(body_conditions[,1]<= bcvector, bcvector >= body_conditions[,2])
  return (row_number())
}
#######################################################################
#Function Name: bc_to_loop
#Purpose: used to loop over "to" body conditions. Calculates how many
#         birds move from focal node to each node in window
#
#Inputs: pm.by.movements (product of window_movements and bird_pop.pm
#           describes number of birds in each body condition that move to 
#           each node within window based on movement probabilities)
#        body_change_logical (logical array describing whether or not birds
#           can move between nodes and body conditions based on tank)
#        new.bc.bin - body condition after foraging
#        bird_pop.am - number of birds in window in am
#        start_x, end_x, start_y, end_y - size of window in bird_pop array
#        bc_to - "to" body conditions
#
#Outputs: BP.am.out = number of birds in each body condition and node in window
#                     on next day of migrations (am birds)
#        
######################################################################

bc_to_loop <- function(pm.by.movements,body_change_logical, new.bc.bin, bp.am, 
                       start_x, end_x, start_y, end_y,bc_to){
  
  #birds that end up in each node
  Added_birds= pm.by.movements*body_change_logical[new.bc.bin,bc_to,,]
  
  #add these birds to any birds already in node.  Only include part of Added_bird
  #within the smaller bird_pop window
  
  BP.am.out =bp.am[,,bc_to]+
    Added_birds[start_x:end_x, start_y:end_y]
  
  return(BP.am.out)
}

#######################################################################
#Function Name: window_setup
#Purpose: used to create window of input layes
#
#Inputs: A = forage.by.day
#        B = roosting.matrix.full
#        C = distance_to_breed.matrix.full
#        D = distance_to_arrival
#        E = WSI.array
#
#Outputs: window_movements = probability of moving to each node within
#         the window
#        
######################################################################
window_setup <- function(A, B, C, D, E) {
  #create a window around focal node
  forage.window = A[lower_bound_x[i]:upper_bound_x[i], lower_bound_y[i]: upper_bound_y[i],d]
  roosting.window=B[lower_bound_x[i]:upper_bound_x[i],lower_bound_y[i]:upper_bound_y[i]]
  distance_to_breed.window=C[lower_bound_x[i]:upper_bound_x[i],lower_bound_y[i]:upper_bound_y[i]]
  
  #Normalize what is in the window
  forage.window.norm=forage.window/sum(sum(forage.window,na.rm=TRUE))
  roosting.window.norm=roosting.window/sum(sum(roosting.window,na.rm=TRUE))
  
  #Distance to breed has opposite direction of imporance (closer is better)
  max.distance = max(max(distance_to_breed.window, na.rm=TRUE))
  distance_to_breed.norm =(max.distance-distance_to_breed.window)/sum(sum((max.distance-distance_to_breed.window),na.rm=TRUE))
  
  #WSI should be mask to remove "poor weather" nodes
  #Using a WSI cutoff of 7.5 for mallards from the Schummer et al. 2010 paper
  
  WSI.array.window = E[lower_bound_x[i]:upper_bound_x[i], lower_bound_y[i]:upper_bound_y[i],]
  WSI.mask = (WSI.array.window[,,d]<=WSI.cutoff)*1
  
  #Use Cobb-Douglas function to combine components of "attractiveness"
  window_movements = WSI.mask*((forage.window.norm)^W1 * (roosting.window.norm)^W2 *
                                 (distance_to_breed.norm)^W3 * (D)^W4)
  
  #Normalize probability over window
  window_movements = window_movements/sum(sum(window_movements, na.rm=TRUE),na.rm=TRUE)
  
  return(window_movements)
}

#From body condition loop function
bc_from_loop <- function(bird_cond) {
  
  #Check that there are birds in the current body condition
  if (bird_pop.am[d,x[i],y[i],body_cond]>0){
    
    #Assume the birds at the beginning of the day forage for the day. Then decide to stay or leave at night
    
    #Determine new body condition based on daily gain
    new.bc = body_conditions[body_cond,3]+DG.matrix.full[focal_node_x[i],focal_node_y[i]]
    
    #Figure out which bin the new body condition is in
    new.bc.bin=new_body_condition(body_conditions, new.bc)
    
    #Compute landscape energetics required to fuel birds:
    #Convert grams to kJ using 37 kJ/g lipids
    energy_needed = DG.matrix.full[focal_node_x[i],focal_node_y[i]]*
      bird_pop.am[d,x[i],y[i],body_cond]*37
    
    #keep summing energy_need over body conditions
    energy_needed.total = energy_needed.total+energy_needed
    
    #Surviving birds at end of the day of foraging.
    #mult am birds by energy_needed/available forage and daily survival
    #if available forage exceeds energy needed, birds removed only due to daily survival
    #otherwise, only proportion of birds that can use available forage stay
    
    bird_pop.pm[d,x[i],y[i],new.bc.bin]=bird_pop.am[d,x[i],y[i],body_cond]*
      ifelse(forage.by.day[focal_node_x[i],focal_node_y[i],d]==0,0, 
             max(1,energy_needed/forage.by.day[focal_node_x[i],focal_node_y[i],d], na.rm=TRUE))*daily_survival[body_cond,1]
    
    #Reduce forage for next body conditions of birds
    forage.by.day[focal_node_x[i],focal_node_y[i],d] = max(0,(forage.by.day[focal_node_x[i],focal_node_y[i],d]-energy_needed))
    
    #Compute probabily of departing
    prob_depart = weight.WSI*prob_depart.WSI[focal_node_x[i],focal_node_y[i],d]+weight.BC*(body_cond^8/(body_cond^8+17^8))
    
    #Generate a random 0-1 and see if it is below the prob_depart.  If so, birds leave.  Else, birds stay
    if (runif(1)<prob_depart) { 
      #Birds leave
      #Move birds according to window_movements and body_change_logical
      
      #Compute pm.by.movements as product of pm birds and window_movments
      pm.by.movements=bird_pop.pm[d,x[i],y[i],new.bc.bin]*window_movements
      
      #use sapply to compute number of birds that end up in each "to" body condition:
      bird_pop.am[d+1,lb_x[i]:up_x[i], lb_y[i]:up_y[i], ]= 
        sapply(seq_len(body_condition_bins),bc_to_loop, pm.by.movements=pm.by.movements,body_change_logical=
                 body_change_logical,new.bc.bin=new.bc.bin, bp.am=bird_pop.am[d+1,lb_x[i]:up_x[i],lb_y[i]:up_y[i],], 
               start_x=lb_x_small[i], end_x=up_x_small[i], start_y=lb_y_small[i], end_y=up_y_small[i], simplify="array")
      
      
    } else {# if loop
      
      #Move all pm birds to next day am
      bird_pop.am[d+1,x[i],y[i],new.bc.bin]=bird_pop.pm[d,x[i],y[i],new.bc.bin]
    }
    
    #Any other case in which birds leave?  In prototype, they leave if they are in highest body condition class
    
    
  } #else loop
  return(list(energy_needed.total, bird_pop.pm, forage.by.day, bird_pop.am))
  
}

#####################################################################
# BIRD MOVEMENT FOR LOOP
#####################################################################

#Loop over days and nodes
Rprof("ForLoop.out",memory.profiling=FALSE)
start.time = Sys.time()
pb1=winProgressBar(title="Days",label="0% done",min=0,max=100,initial=0)
target=2
#for (d in 1:num_days) {  #Number of days of WSI data
for (d in 1:target) {
  timeit<-proc.time()
  start.time=Sys.time()
  pb2=winProgressBar(title="Nodes",label="0% done",min=0,max=100,
                     initial=0)
  for (i in 1:nrow(data_nodes) ) { #Loop through data nodes
    
    #Check if there are birds in the node at the start of the day  
    if (sum(bird_pop.am[d,focal_node_x[i]-window_radius,
                        focal_node_y[i]-window_radius,],
            na.rm=TRUE) >0 ) {
      #There are birds in the node. 
      window_movements = window_setup(forage.by.day, 
                                      roosting.matrix.full, 
                                      distance_to_breed.matrix.full, 
                                      distance_to_arrival, WSI.array)
      
      energy_needed.total=0 #clear out energy needed. add up as loop 
                            #goes
      
     #compute new body condition due to daily gain
      bcvector<- (body_conditions[,3]+DG.matrix.full[focal_node_x[i],
                                                     focal_node_y[i]])
      
      #determine which bin the new body condition falls in
      new.bc.bin = findInterval(bcvector, body_conditions[,1])
      
      #Compute landscape energetics required to fuel birds:
      #Convert grams to kJ using 37 kJ/g lipids
      energy_needed = 
        DG.matrix.full[focal_node_x[i],focal_node_y[i]]*
        bird_pop.am[d,x[i],y[i],]*37
      
      #keep summing energy_need over body conditions
      energy_needed.total = energy_needed.total+sum(energy_needed)
      
      #Surviving birds at end of the day of foraging.
      #mult am birds by energy_needed/available forage and daily 
      #survival if available forage exceeds energy needed, birds 
      #removed only due to daily survival otherwise, only proportion 
      #of birds that can use available forage stay
      
      PMPop = bird_pop.am[d,x[i],y[i],] * 
        ifelse(forage.by.day[focal_node_x[i],
                             focal_node_y[i],d]==0,0,
               max(1,energy_needed/forage.by.day[focal_node_x[i],
                                                 focal_node_y[i],d], 
                   na.rm=TRUE))*daily_survival[,1]
      BPAgg = aggregate(cbind(PMPop,new.bc.bin),by=list(new.bc.bin), 
                        FUN=sum)
      bird_pop.pm[d,x[i],y[i],BPAgg[,1]] = BPAgg[,2]

      #Reduce forage for next body conditions of birds
      forage.by.day[focal_node_x[i],focal_node_y[i],d] = 
        max(0,(forage.by.day[focal_node_x[i],focal_node_y[i],d]-
                 energy_needed))
      
      #Compute probabily of departing
      prob_depart = 
        weight.WSI*prob_depart.WSI[focal_node_x[i],focal_node_y[i],d]+
        weight.BC*(seq(1:body_condition_bins)^8/
                     (seq(1:body_condition_bins)^8+17^8))
      
      #Generate a random 0-1 and see if it is below the prob_depart. 
      #If so, birds leave.  Else, birds stay leave is a vector of 
      #indices describing which body condition classes will leave
      leave = which(runif(body_condition_bins)<prob_depart)  
      stay = setdiff(seq(1:body_condition_bins),leave)
      
      #For each body condition class that leaves
      #Move birds according to window_movements and body_change_logical
      for (bc.to in leave) {
        pm.by.movements =
          bird_pop.pm[d,x[i],y[i],bc.to]*window_movements
        bird_pop.am[d+1,lb_x[i]:up_x[i],lb_y[i]:up_y[i],]= 
          bird_pop.am[d+1,lb_x[i]:up_x[i],lb_y[i]:up_y[i],]+
          sapply(seq_len(body_condition_bins),bc_to_loop,
                 pm.by.movements=pm.by.movements,
                 body_change_logical=body_change_logical,
                 new.bc.bin=bc.to,
                 bp.am=bird_pop.am[d+1,lb_x[i]:up_x[i],
                                   lb_y[i]:up_y[i],],
                 start_x=lb_x_small[i],end_x=up_x_small[i],
                 start_y=lb_y_small[i],end_y=up_y_small[i],
                 simplify="array")
      }
      
      #If birds didn't leave, move all pm birds to next day am
      bird_pop.am[d+1,x[i],y[i],stay]=
        bird_pop.am[d+1,x[i],y[i],stay]+bird_pop.pm[d,x[i],y[i],stay]
      
    } #else loop
    info=sprintf("%d%% done",round((i/nrow(data_nodes))*100))
    setWinProgressBar(pb2,i/(nrow(data_nodes))*100,label=info)
  } #for loop
  close(pb2)
  #reduce the forage on the next day by amount of energy_needed
  forage.by.day[focal_node_x[i],focal_node_y[i],d+1]=
    max(0,(forage.by.day[focal_node_x[i],focal_node_y[i],d]-
             energy_needed.total))
time.take=Sys.time()-start.time
time.take
print(d)

info=sprintf("%d%% done",round((d/target)*100))
setWinProgressBar(pb1,(d/target)*100,label=info)

} #d loop
#} # d loop
close(pb1)
summaryRprof("ForLoop.out")

bird_pop.am.sum=apply(bird_pop.am,1:3,sum, na.rm=TRUE)
bird_pop.pm.sum=apply(bird_pop.pm,1:3,sum, na.rm=TRUE)

for (d in 1:100){
  image.plot(bird_pop.am.sum[d,,])
  image.plot(bird_pop.pm.sum[d,,])
}


#######################################################################
#My attempt at creating c++ function
######################################################################
#library(Rcpp)

#cppFunction('NumericVector bc_to_loop_c(NumericVector pmbymovements, NumericMatrix body_change_logical_c, int newbcbin) {
#             NumericVector out;
#             for (int i=0; i<20; i++) {
#             NumericMatrix Added_birds = pmbymovements*body_change_logical_c(newbcbin,i-1,_,_);
#             out[i]=Added_birds;
#             }
#             return out;
#             }')            

#             NumericMatrix BPamout = bird_pam + Added_birds(Range(start_x,end_x), Range(start_y,end_y));
#             return BPamout;
#             }'
#             
#)
