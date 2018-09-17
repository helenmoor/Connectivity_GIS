## Connectivity_Example.R
# Helen Moor 290618
# Calculating connectivity of focal points from age, volume and distance of forest, 
# following method in Mair et al. 2017, Ecology & Evolution:
# sum up spruce volume in pixels of age > 100 yrs, weighted by negative exponential distance decay from focal site.
# Example below for Finnish NFI data.


# INPUT:
#- point data focal sites
#- raster stand age (yrs)
#- raster spruce volume (m3/ha)


## PROCEDURE:
# - select pixels of age >= 100 yrs
# - sum spruce volume in these pixels; aggregate to 100 m pixel (using average) ->  spruce vol m3/ha of old forest
# - find distances from focal cell to all 100 m cells within 20 km buffer 
# - calculate connectivity as S-i = sum(Sprucevol_j * exp(-alpha * dist_ij)) for different values of alpha 
# - alpha values: alphas=c(2,1,0.5, 0.2,0.1)  # for 0.5, 1, 2 ,5, 10 km mean dispersal 


## NOTES:
# Here we used Finnish NFI data from 2009:
# Fin NFI 2009 data:
# - crs= EPSG 3067: Proj4js.defs["EPSG:3067"] = "+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs";
# - original resolution 20 x 20 m
# - Nodata values: -32768, 32766, 32767



### libraries:

library(raster)
library(gdalUtils)
library(rgdal)


### set your workdir:
setwd("~")

#start.time <- Sys.time() 

# set coordinate system: crsdef: EPSG:3067
mycrs <- crs( "+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs")

# site coordinates:
# get site point data:
sites <- read.csv(file="FungiSites_all180621/FIN553Sites_EPSG3067.csv", stringsAsFactors = F)
sites <- subset(sites, select = c( "StandID"  ,  "X_EPSG3067"  , "Y_EPSG3067" )) # keep Stand ID and coordinates

# make dataframe for output:
sitedata <- sites 

# make sites a SpatialPointsDataFrame
coordinates(sites) <- c("X_EPSG3067"  , "Y_EPSG3067" ) 
proj4string(sites) <- mycrs  # set coordinate system for sites


#####################################################################################
### RASTER DATA PREPARATION:----------------------------------------------------------

# merged from multiple tiles, edited and stored as brick
# first merge tiles onto one raster:

### STANDAGE:--------

#Build list of all raster files you want to join (in your current working directory).
setwd("Finland_LUKE/2009/Standage/2009/")
all_my_rasts <- list.files(pattern = "\\.tif$") # 13 files

# Make a template raster file to build onto.
# required extent was determined manually in QGis
# this is covering the whole area: (38592 x  30000 pixels)
e <- extent(116000, 733472, 6666000, 7146000) #xmin, xmax, ymin, ymax
template <- raster(e)
projection(template) <- mycrs
writeRaster(template, file="Standage2009.tif", format="GTiff") # 2.3 Gb single raster

# Merge all raster tiles into one big raster.
mosaic_rasters(gdalfile=all_my_rasts,dst_dataset="Standage2009.tif",of="GTiff")
gdalinfo("Standage2009.tif")
rm(e,template, all_my_rasts)


### SpruceVolume:---------

#Build list of all raster files you want to join (in your current working directory).
setwd("Finland_LUKE/2009/Sprucevol/2009/")
all_my_rasts <- list.files(pattern = "\\.tif$")

# Make a template raster file to build onto.
# required extent was determined manually in QGis
e <- extent(116000, 733472, 6666000, 7146000) #xmin, xmax, ymin, ymax
template <- raster(e)
projection(template) <- mycrs
writeRaster(template, file="Sprucevol2009.tif", format="GTiff")

# Merge all raster tiles into one big raster.
mosaic_rasters(gdalfile=all_my_rasts,dst_dataset="Sprucevol2009.tif",of="GTiff")
gdalinfo("Sprucevol2009.tif")
rm(e,template, all_my_rasts)
gc()


### COMBINE STAND AGE AND SPRUCE VOL:----------
# Load both rasters into r:
#setwd("~")

age <- raster("Standage/2009/Standage2009.tif")
spruce <- raster("Sprucevol/2009/Sprucevol2009.tif")

# set Nodata values to NA: ( -32768, 32766, 32767)  (use reclassify for better memory management)
age <- reclassify(age, cbind(c(-32768, 32766, 32767), NA))
spruce <- reclassify(spruce, cbind(c(-32768, 32766, 32767), NA))


## Store merged rasters as brick  
nfi2009 <- brick(age,spruce)
names(nfi2009) <- c("age", "sprucevol")

# # store the brick for later use:
# setwd("~")
# writeRaster(nfi2009, filename="NFI2009.grd", format="raster")
rm(age, spruce)
gc()

# Load the brick:
# setwd("~")
# nfi2013 <- brick('NFI2009.grd')
#  - OR continue with loaded brick from work environment.


#####################################################################################
#### Connectivity calculation: ---------------------------------------------

# PROCEDURE:
# 1. get aggregated spruce volume in old pixels:
# - us age as mask to find pixels in spruce with age >= 100 yrs:
# - aggregate to 100m pixels  using average

#start.time <- Sys.time()

# first select only relevant pixels (old cells, age >= 100yrs), set rest to 0
spruceagg <- nfi2009$sprucevol
spruceagg[nfi2009$age < 100] <- 0;

rm(nfi2009) # not required any longer
gc()

## aggregate resulting spruce vol in old pixels to 100m resolution (mean)
spruceagg <- aggregate(spruceagg, fact=5, fun=mean, expand=TRUE, na.rm=TRUE)

# # this may takes quite a while
# end.time <- Sys.time()
# time.taken <- end.time - start.time


# spruceagg now has low resolution average vol of spruce in old cells:
# save to spare the time for aggregation:
#writeRaster(spruceagg, file="2009_meanSpruceVolOld.tif", format="GTiff")

# # load the aggregated old spruce volumes:
# setwd("/Users/Helen/Documents/work_local/_FUNGI/DATA2/Finland_LUKE/2009")
# spruceagg <- raster("2009_meanSpruceVolOld.tif")


# 2. calculate connectivity for all stands and different alpha values:
# - calculate connectivity as S-i = sum(Sprucevol_j * exp(-alpha * dist_ij));
# - calculate distance from each focal point to all 100m grid cells within buffer distance
# - buffer distance 20 km around focal site
bufferdist=20000
# different alpha values used:
alphas=c(2,1,0.5, 0.2,0.1)  # for 0.5, 1, 2 ,5, 10 km mean dispersal 


# prepare sitedata for output:
sitedata$conn_20km_a2 <- NA
sitedata$conn_20km_a1 <- NA
sitedata$conn_20km_a05 <- NA
sitedata$conn_20km_a02 <- NA
sitedata$conn_20km_a01 <- NA


#start.time <- Sys.time() 

## loop through alphas and sites:
for(k in 1:length(alphas)){
  alpha <- alphas[k]

  for(i in 1:length(sites$StandID)){

    # distance from site (point) to (centroid) of all pixels:
    dist <- distanceFromPoints(spruceagg, sites[i,]) 
    
    # distance m -> km:
    dist <-dist/1000

    # apply distance decay: weight spruceagg separately by distance decay for each point:
    temp <- spruceagg * exp(-alpha*dist)

    # sum values by extracting within buffer from this distance-weighted raster for each point:
    sitedata[i,4+k]  <- extract(temp, sites[i,], method='simple', buffer=bufferdist,fun=sum)
  }
}

## store collected connectivity measures:
write.csv(sitedata,file="Connectivity_20km_100magg_2009.csv", row.names = F)

# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken

rm(list=ls())
gc()



