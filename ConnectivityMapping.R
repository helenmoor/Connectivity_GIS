#####################################################################################
# ConnectivityMapping.R
# Calculation and mapping of connectivity based on species dispersal kernels
# HM, Sept 2017
#####################################################################################
# INPUT data:
# - KNN raster files for forest age and volume of target tree species (spruce, pine,...)
#   - Skogskarta 2010: resolution 25m, coordinate system RT90 (EPSG:3021)
#   - Here, the input raster files have already been clipped to the county border plus a surrounding 
#     buffer area of 20km width (to minimise edge effects).
#   -> AGE_25mRes_RT90_DalWBuffer.tif and Pinevol_25mRes_RT90_DalWBuffer.tif
# - vector file outline of county (here: Dalarna)


# Calculate connectivity for Dalarna using connectivity metric of Mair et al. 2017,
# using a moving window, with function focal().


###  Load required libraries:-------------------------------------
library(raster)
library(gdalUtils)
library(rgdal)


### Required settings:--------------------------------------------


## A) Set threshold radius of search window:
# - threshold radius for moving window (should ideally equal the buffer distance around the county borders);
# - see Table 1 in the protocol for minimum distances, or calculate MinDist below.
# (A greater than needed distance is not a problem except in terms of increased computing time.)

# As an indication of the minimum required distance, if using the negative expontential kernel, calculate the 
# distance from the focal site at which the Probability of dispersal P is <= 1e-06 (depends on Mean dispersal distance):

# - for given mean dispersal distance (= 1/alpha),
MeanDispersalDistance = 5  # unit: km; MeanDispersalDistance = 1/alpha;
# - calculate minimum threshold distance where dispersal probability is below 1e-06 :
MinDist = -log(1e-06)*MeanDispersalDistance
MinDist

# -> Set chosen radius equal or larger MinDist: (unit m, same as map projection)
myradius <- 20000

# Adjust buffer area of input raster accordingly 
# (or work with whole raster if needed, but this will take several hours and might exceed available memory.)


## B) Mean dispersal distances of interest:
# - (dispersal mean distance = 1/alpha)
# - NOTE: choose your threshold radius adjusted to the largest mean dispersal distance used here!
# set alphas:
alphas=c(5, 1, 0.2)  # corresponds to mean distance  = 0.2, 1, 5 km



## C) Load INDATA:
# - 25m resolution, CRS RT90
# - masked to the contour of Dalarna including a  buffer strip around the border covering a distance of myradius or greater
# ( Note: additional indata could be forest continuity, available from http://gpt.vic-metria.nu/data/land/Dalarnas_lan.zip , 
#   see documentation here: http://gpt.vic-metria.nu/data/land/Slutrapport_Kartering_av_kontinuitetsskog_boreal_region_20170117.pdf)

# Navigate to your working directory where raster files are stored:
setwd("DATA/GIS/Sweden/SLU_skogskarta2010_Rt90")

# read in the rasters:
age <- raster('AGE_25mRes_RT90_DalWBuffer.tif')
vol <- raster('PINEVOL_25mRes_RT90_DalwBuffer.tif')

# Coordinate system: 
# - first ensure that your rasters are indeed all in RT90 (or SWEREF99TM)
# - then make sure all rasters use the same definition (may differ in decimal places)

# RT90 2.5 gon V definition: 
rt90 <- crs("+proj=tmerc +lat_0=0 +lon_0=15.808277778 +k=1 +x_0=1500000 +y_0=0 +ellps=bessel +towgs84=414.1,41.3,603.1,-0.855,2.141,-7.023,0
+units=m +no_defs")

# Sweref99TM definition:
# sweref <- crs("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

crs(age) <- rt90
crs(vol) <- rt90



### Connectivity calculation:--------------------------------------

## 1. select old cells and aggregate to 100m:----------------------------------------------------

# select old cells (where age >100 yrs), set remaining cells to 0 (nodata value; don't use NA!):
vol[age<100] <- 0; # this can take a few minutes

# NOTE: Here, old cells could be weighted additionally by their continuity. E.g. multiply volume with weighting value that increases for increased continuity.

# then reduce the resolution by aggregating (sum) to 100m pixels (for computational efficiency):
vol <- aggregate(vol, fact=4, fun=sum, expand=TRUE, na.rm=TRUE)
 
# divide by nr of pixels (4*4=16) for average volume of spruce older than 100 yrs: 
vol <- vol/16


# store this intermediate step (subset to volume of old cells, nodata value = 0) as R format raster, to load again if needed:
#writeRaster(vol, filename = "PINEVolumeOldCells_Res100m.grd", format="raster", overwrite=TRUE)
# load intermediate step raster:
#vol <- raster("PINEVolumeOldCells_Res100m.grd")

# have a look at the raster:
# vol
# plot(vol)



## 2. construct distance weights moving window matrix:------------------------------------------

# create a moving window for distances:
# final unit: kilometer

# current resolution: 100m , crs in meter
myres <- res(vol)[1]
# desired radius: 20km = 20000m  
myradius 

# create an empty raster with desired extent and resolution:
# matrix dimensions: 2 * desired radius (in crs unit)/ resolution (in crs units) + 1 (for center)
# desired extent (in crs units):  -(radius+res/2)   +radius+res/2, setting xmin at half a cell centers the matrix on centroids of cells
win <- raster( xmn=-(myradius + myres/2), xmx= +(myradius + myres/2), ymn=-(myradius + myres/2), ymx=+(myradius + myres/2),  res = myres)
# specify the center point:
p1 <-   SpatialPoints(data.frame(X=0, Y=0))  

# calculate distance matrix:
distances <- distanceFromPoints(win, p1) 
#as.matrix(distances)[195:205, 195:205] # have a look at the center

# convert from m units to km distances:
distances <- distances/1000

# convert to matrix:
distances <- as.matrix(distances)
#distances[201,201] # center of window



## 3. load county border for masking:------------------------------------------------------------------

# load county border outline for cropping/masking:
dal <- readOGR(dsn = "/DATA/GIS/Sweden/.", 
                  layer = "Dalarna_border_RT90",
               p4s="+proj=tmerc +lat_0=0 +lon_0=15.808277778 +k=1 +x_0=1500000 +y_0=0 +ellps=bessel +towgs84=414.1,41.3,603.1,-0.855,2.141,-7.023,0
+units=m +no_defs")
crs(dal) <- crs(conn)


## 4. calculate connectivity and export final masked rasters:------------------------------------------

# set working directory: where do you want to store your output rasters?
setwd("DATA/GIS/Sweden/ConnectivityDalarna/ConnectivityMapsDalarna_v2/PineConn/")

# set threshold radius in km: 
radius = myradius/1000


# Loop over chosen alphas:

for(k in 1:length(alphas)){
  
  alpha<- alphas[k]
  
  ## dispersal kernel: calculate decay weights from distances (depending on alpha):
  weights <- exp(-alpha*distances)
  
  #  radius of 20km - NOTE that the moving window is square, larger at edges -> mask to circle using threshold value
  # Threshold value for weights (depending on threshold radius):
  threshold <- exp(-alpha*radius)
  #  Set weights to zero outside of threshold radius (depends on alpha)
  weights[weights < threshold] <- 0
  
  ## area included:
  myarea <- pi*radius^2
  
  # store in file (not yet standardised or masked to county borders) - to store intermediate step: (can be skipped)
  myfilename = paste("conn_",radius,"km_alpha", alpha, ".tif", sep="")  
  
  # calculate connectivity using weights matrix:
  conn <- focal(vol, w = weights, fun=sum, filename=myfilename, na.rm=TRUE, pad=TRUE, padValue=0) 
  
  # crop and mask by Dalarna real outline
  conn <- crop(conn, extent(dal))
  conn <- mask(conn, dal)

  # standardize by area and store standardised raster (final result):
  conn <- conn/myarea
  
  writeRaster(conn, filename = paste("conn_",radius,"km_alpha", alpha, "_stand.tif", sep=""), format = "GTiff", overwrite=TRUE  )
  
  
}

