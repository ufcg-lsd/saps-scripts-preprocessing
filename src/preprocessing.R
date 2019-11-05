########################################################################################
#                                                                                      #
#                         EU BRAZIL Cloud Connect                                      #
#                                                                                      #
#                                                                                      #
########################################################################################

options(echo=TRUE)
rm(list=ls())

library(gdalUtils)
library(R.utils)
library(raster)
library(rgdal)
library(ncdf4)
library(sp)
library(snow)
library(maptools)

args = commandArgs(trailingOnly=TRUE)
WD <- args[1]
setwd(WD) # Working Directory

# Changing raster tmpdir

rasterOptions(tmpdir=args[2])
# Load the source code in landsat.R to this code
source("landsat.R")

# File that stores the Image Directories (TIFs, MTL, FMask)
dados <- read.csv("dados.csv", sep=";", stringsAsFactors=FALSE)
#################################### Constants ##########################################

k <- 0.41		# Von K?rm?n
g <- 9.81		# Gravity
clusters <- 7		# Number of clusters used in image processing - some raster library methods are naturally coded to run in a clustered way

######################### Reading sensor parameters #####################################

p.s.TM1 <- read.csv("parametros do sensor/parametrosdosensorTM1.csv", sep=";", stringsAsFactors=FALSE)
p.s.TM2 <- read.csv("parametros do sensor/parametrosdosensorTM2.csv", sep=";", stringsAsFactors=FALSE)
p.s.ETM <- read.csv("parametros do sensor/parametrosdosensorETM.csv", sep=";", stringsAsFactors=FALSE)
p.s.LC <- read.csv("parametros do sensor/parametrosdosensorLC.csv", sep=";", stringsAsFactors=FALSE)

# Read relative distance from Sun to Earth
load("d_sun_earth.RData")

# Set projection and spatial resolution
WGS84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84"

######################### Image Information ######################################

fic.dir <- dados$File.images[1]  # Images Directory
MTL <- read.table(dados$MTL[1], skip=0, nrows=140, sep="=", quote = "''", as.is=TRUE) # Reading MTL File

fic <- substr(MTL$V2[MTL$V1 == grep(pattern="LANDSAT_SCENE_ID", MTL$V1, value=T)], 3, 23)

n.sensor <- as.numeric(substr(fic, 3, 3)) # Sensor Number

if (n.sensor==8) MTL <- read.table(dados$MTL[1], skip=0, nrows=-1, sep="=", quote="''", as.is=TRUE, fill=TRUE) # Reading MTL File for Sensor number 8

WRSPR <- substr(fic, 4, 9)						#WRSPR
PATH <- substr(WRSPR, 0, 2)
ROW <- substr(WRSPR, 3, 5)
print (PATH)
print (ROW)
Ano <- as.numeric(substr(fic, 10, 13))			#Images year
Dia.juliano <- as.numeric(substr(fic, 14, 16))	#Julian Day

# Getting the sum elevation at the time of Image Capture
sun_elevation <- as.numeric(MTL$V2[MTL$V1 == grep(pattern="SUN_ELEVATION", MTL$V1, value=TRUE)])
costheta <- sin(sun_elevation*pi/180) # From SUN ELEVATION

# Setting the sensor parameter by the Sattelite sensor type and data
if (n.sensor==8) p.s <- p.s.LC
if (n.sensor==7) p.s <- p.s.ETM
if (Ano < 1992 & n.sensor==5) p.s <- p.s.TM1 
if (Ano > 1992 & n.sensor==5) p.s <- p.s.TM2

# Time image
acquired_date <- as.Date(MTL$V2[MTL$V1==grep(pattern="DATE_ACQUIRED", MTL$V1, value=TRUE)])
daysSince1970 <- as.numeric(acquired_date)
tdim <- ncdim_def("time", "days since 1970-1-1", daysSince1970, unlim=TRUE, create_dimvar=TRUE, "standard", "time")

# Reading image file
# The Images are of the type ".tif" that represents each spectral band captured by the satellite
# Depending on the sattelite the number of spectral bands captured are differents
fichs.imagens <- list.files(path=fic.dir, pattern="*.TIF", full.names=TRUE)

getBandsPath <- function(n.sensor){
  wanted_bands <- NULL
  if(n.sensor == 8) {
    wanted_bands <- c("B2", "B3", "B4", "B5", "B6", "B7", "B10")
  } else if(n.sensor == 7) {
    wanted_bands <- c("B1", "B2", "B3", "B4", "B5", "B6_VCID_1", "B6_VCID_2", "B7")
  } else if(n.sensor == 5) {
    wanted_bands <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7")
  }
  
  bands_path <- list()
  for (i in 1:length(wanted_bands)) {
    for (j in 1:length(fichs.imagens)) {
      if(regexpr(paste(wanted_bands[i], '.TIF', sep=""), fichs.imagens[[j]]) != -1) {
        bands_path[[i]] <- fichs.imagens[[j]]
      }
    }
  }
  
  return(bands_path)
}

# Reading
fic.st <- stack(as.list(getBandsPath(n.sensor)))
wanted_bands.path <- getBandsPath(n.sensor)

print("LABEL - DATA READING")
proc.time()

################################# Fmask ###########################################

# Identifier of clouds and shadows
# The FMask serves to identify if exist any cloud or shadow on the image, the existence of clouds or shadow disturbs the results.

n.fmask <- length(fichs.imagens)
Fmask <- raster(fichs.imagens[[n.fmask]])
fmask <- as.vector(Fmask)

if (n.sensor != 8) mask_filter <- 672 else 
  mask_filter <- 2720
for (i in 1:nlayers(fic.st)) {
  f <- fic.st[[i]][]
  f[fmask != mask_filter] <- NaN
  fic.st[[i]][] <- f 
}

proc.time()

if(n.sensor == 7){
  if (0.99<=(sum(is.na(values(fic.st)))/8)/(fic.st@ncols*fic.st@nrows)) {
    print("Imagem incompativel para o processamento,mais de 99% nuvem e sombra de nuvem")
    quit("no", 1, FALSE)
  }
}else{
  if (0.99<=(sum(is.na(values(fic.st)))/7)/(fic.st@ncols*fic.st@nrows)) { 
    print("Imagem incompativel para o processamento,mais de 99% nuvem e sombra de nuvem")
    quit("no", 1, FALSE)
  }
}

print("LABEL - MASKING CLOUDS")
proc.time()

#WRITING WITHOUT CLOUDS
bands.path <- c()
for(i in 1:length(wanted_bands.path)){
	aux <- sub('\\.TIF', '_noCloud.tif', wanted_bands.path[[i]])
	bands.path <- c(bands.path, aux)
	writeRaster(fic.st[[i]], aux, format='GTiff', NAflag=0, overwrite=TRUE)
}

print("LABEL - WRITING WITHOUT CLOUDS")
proc.time()

# Changing the projection of the images (UTM to GEO)
# This operation can be done in a parallel way by Clusters, projectRaster is implemented to naturally be executed by clusters
# The number of used clusters is given by the 'clusters' constant

s_srs_2 <- paste(crs(fic.st[[1]]))
fic.st <- projectRaster(fic.st, crs=WGS84)

print("LABEL - LANDSAT IMAGES PROJECT")
proc.time()

# Reading Bounding Box
# The Bounding Box area that is important and has less noise in the Image
fic.bounding.boxes <- paste("wrs2_asc_desc/wrs2_asc_desc.shp")
BoundingBoxes <- readOGR(dsn=fic.bounding.boxes)#, proj4string=CRS(WGS84))
BoundingBox <- BoundingBoxes[BoundingBoxes@data$WRSPR == WRSPR, ]

# Reading Elevation
# Read the File that stores the Elevation of the image area, this influence on some calculations
fic.elevation <- paste(fic.dir, "/", WRSPR, ".tif", sep="")
raster.elevation <- raster(fic.elevation)

print("LABEL - BEFORE ELEVATION RESAMPLE")
proc.time()

# Setting the raster elevation resolution as equals to the Fmask raster resolution
tr <- res(fic.st)
s_srs <- paste(crs(raster.elevation))
t_srs <- WGS84

rm(fic.st)
gc()

condition <- paste('PATH=',PATH,' AND ROW=',ROW, sep="")
print (condition)
gdalwarp(fic.elevation, sub('\\.tif', '_RESAMPLED.tif', fic.elevation), s_srs=s_srs, t_srs=t_srs, tr=tr, cutline=fic.bounding.boxes, cwhere=condition, crop_to_cutline=TRUE, overwrite=TRUE, verbose=TRUE)

raster.elevation <- raster(sub('\\.tif', '_RESAMPLED.tif', fic.elevation))

print("LABEL - AFTER ELEVATION RESAMPLE")
proc.time()

#################### Resampling satellite bands images #####################################

# This block of code resample the image based on the Elevation of the terrain captured by the sat
# The Elevation of the terrain needs to be taken into account

fic.st.paths <- c()
for(i in 1:(length(bands.path))){
	fic.st.paths <- c(fic.st.paths, sub('\\_noCloud.tif', '_RESAMPLED.TIF', bands.path[[i]]))
}

tr <- res(raster.elevation)

sapply(bands.path, function(file){
	gdalwarp(file, sub('\\_noCloud.tif', '_RESAMPLED.TIF', file), s_srs=s_srs_2, t_srs=t_srs, tr=tr, cutline=fic.bounding.boxes, srcnodata='0', dstnodata='0', cwhere='PATH=215 AND ROW=65', crop_to_cutline=TRUE, overwrite=TRUE, verbose=TRUE)
})

print("LABEL - BANDS RESAMPLE")
proc.time()

image.rec <- stack()
for(i in 1:length(fic.st.paths)){
	aux <- raster(fic.st.paths[[i]])
	image.rec <- stack(image.rec, aux)
}

print("LABEL - NEW BANDS READING")
proc.time()

############################################################################################

# Reading file Station weather
fic.sw <- dados$File.Station.Weather[1]
table.sw <- (read.csv(fic.sw, sep=";", header=FALSE, stringsAsFactors=FALSE))
hour.image <- (as.numeric(substr(MTL$V2[MTL$V1 == grep(pattern="SCENE_CENTER_TIME", MTL$V1, value=T)], 3, 4))+
                 as.numeric(substr(MTL$V2[MTL$V1 == grep(pattern="SCENE_CENTER_TIME", MTL$V1, value=T)], 6, 7))/60)*100
hour.image.station<-which.min(abs(table.sw$V3[]-hour.image))

# Transmissivity 
tal <- 0.75+2*10^-5*raster.elevation

proc.time()

################## Phase 1: Calculating the image energy balance ##################################

# This block calculate the energy balance of the image

output <- NULL;
outputLandsat <- function() {
  output <- landsat()
  return(output)
}

# timeout before = 2665.151
# timeout now is 7200 (cause: Azure slowness)
res <- NULL;
tryCatch({
  res <- withTimeout({
    output <- outputLandsat();
  }, timeout=7200);
}, TimeoutException=function(ex) {
  cat("Output landsat timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

###########################################################################################

print("LABEL - LANDSAT FUNCTION")
proc.time()

rm(image.rec)
gc()

################## Masking landsat rasters output #########################################

# This block mask the values in the landsat output rasters that has cloud cells and are inside the Bounding Box required
# This block is already Clustered

outputMask <- function() {
  output <- mask(output, BoundingBox)
  return(output)
}

# timeout before = 1716.853
# timeout now is 10800 (cause: Azure slowness)

res <- NULL;
tryCatch({
  res <- withTimeout({
    output <- outputMask();
  }, timeout=10800);
}, TimeoutException=function(ex) {
  cat("Output Fmask timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

##########################################################################################

print("LABEL - MASKING OUTPUT LANDSAT")
proc.time()

################## Write to files landsat output rasters #################################

# This block write landsat outputs rasters to files

output.path<-paste(dados$Path.Output[1], "/", fic, ".nc", sep="")
outputWriteRaster <- function() {
  names(output) <- c("Rn", "TS", "NDVI", "EVI", "LAI", "G", "alb", "SAVI")
  # names(output) <- c("TS", "NDVI","LAI", "alb","SAVI")
  writeRaster(output, output.path, overwrite=TRUE, format="CDF", varname=fic, varunit="daily", longname=fic, xname="lon", yname="lat", bylayer=TRUE, suffix="names")
}

# timeout before = 1708.507
# timeout now is 10800 (cause: Azure slowness)

res <- NULL;
tryCatch({
  res <- withTimeout({
    outputWriteRaster();
  }, timeout=10800);
}, TimeoutException=function(ex) {
  cat("Output write raster timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

####### Saving Albedo ######

# Opening old alb NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_alb.nc", sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

proc.time()

# Getting lat and lon values from old NetCDF
oldLat <- ncvar_get(nc, "lat", start=1, count=raster.elevation@nrows)
oldLon <- ncvar_get(nc, "lon", start=1, count=raster.elevation@ncols)

# Defining latitude and longitude dimensions
dimLatDef <- ncdim_def("lat", "degrees", oldLat, unlim=FALSE, longname="latitude")
dimLonDef <- ncdim_def("lon", "degrees", oldLon, unlim=FALSE, longname="longitude")

proc.time()

# New alb file name
file_output <- paste(dados$Path.Output[1], "/", fic, "_alb.nc", sep="")
oldAlbValues <- ncvar_get(nc, fic)
newAlbValues <- ncvar_def("alb", "daily", list(dimLonDef, dimLatDef, tdim), longname="alb", missval=NaN, prec="double")
nc_close(nc)
newAlbNCDF4 <- nc_create(file_output, newAlbValues)
ncvar_put(newAlbNCDF4, "alb", oldAlbValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newAlbNCDF4)

proc.time()

####### Saving LAI ######

# Opening old LAI NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_LAI.nc", sep="")
nc <- nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# New LAI file name
file_output <- paste(dados$Path.Output[1], "/", fic, "_LAI.nc", sep="")
oldLAIValues <- ncvar_get(nc, fic)
newLAIValues <- ncvar_def("LAI", "daily", list(dimLonDef, dimLatDef, tdim), longname="LAI", missval=NaN, prec="double")
nc_close(nc)
newLAINCDF4 <- nc_create(file_output, newLAIValues)
ncvar_put(newLAINCDF4, "LAI", oldLAIValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newLAINCDF4)

proc.time()

####### Saving NDVI ######

# Opening old NDVI NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_NDVI.nc", sep="")
nc <- nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# New NDVI file name
file_output <- paste(dados$Path.Output[1], "/", fic, "_NDVI.nc", sep="")
oldNDVIValues <- ncvar_get(nc,fic)
newNDVIValues <- ncvar_def("NDVI", "daily", list(dimLonDef, dimLatDef, tdim), longname="NDVI", missval=NaN, prec="double")
nc_close(nc)
newNDVINCDF4 <- nc_create(file_output, newNDVIValues)
ncvar_put(newNDVINCDF4, "NDVI", oldNDVIValues, start=c(1, 1, 1),count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newNDVINCDF4)

proc.time()

###### Saving SAVI ######

#Opening old SAVI NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_SAVI.nc", sep="")
nc <- nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# New SAVI file name
file_output <- paste(dados$Path.Output[1], "/", fic, "_SAVI.nc", sep="")
oldSAVIValues <- ncvar_get(nc, fic)
newSAVIValues <- ncvar_def("SAVI", "daily", list(dimLonDef, dimLatDef, tdim), longname="SAVI", missval=NaN, prec="double")
nc_close(nc)
newSAVINCDF4 <- nc_create(file_output, newSAVIValues)
ncvar_put(newSAVINCDF4, "SAVI", oldSAVIValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newSAVINCDF4)

proc.time()

###### Saving EVI ######

#Opening old EVI NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_EVI.nc", sep="")
nc <- nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# New EVI file name
file_output <- paste(dados$Path.Output[1], "/", fic, "_EVI.nc", sep="")
oldEVIValues <- ncvar_get(nc, fic)
newEVIValues <- ncvar_def("EVI", "daily", list(dimLonDef, dimLatDef, tdim), longname="EVI", missval=NaN, prec="double")
nc_close(nc)
newEVINCDF4 <- nc_create(file_output, newEVIValues)
ncvar_put(newEVINCDF4, "EVI", oldEVIValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newEVINCDF4)

proc.time()

###### Saving TS ######

# Opening old TS NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_TS.nc", sep="")
nc <- nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

#New TS file name
file_output<-paste(dados$Path.Output[1],"/",fic,"_TS.nc",sep="")
oldTSValues<-ncvar_get(nc,fic)
newTSValues<-ncvar_def("TS","daily",list(dimLonDef,dimLatDef,tdim),longname="TS",missval=NaN,prec="double")
nc_close(nc)
newTSNCDF4<-nc_create(file_output,newTSValues)
ncvar_put(newTSNCDF4,"TS",oldTSValues,start=c(1,1,1),count=c(raster.elevation@ncols,raster.elevation@nrows,1))
nc_close(newTSNCDF4)

proc.time()

###### Saving Rn ######

# Opening old Rn NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_Rn.nc", sep="")
nc <- nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

#New Rn file name
file_output<-paste(dados$Path.Output[1],"/",fic,"_Rn.nc",sep="")
oldRnValues<-ncvar_get(nc,fic)
newRnValues<-ncvar_def("Rn","daily",list(dimLonDef,dimLatDef,tdim),longname="Rn",missval=NaN,prec="double")
nc_close(nc)
newRnNCDF4<-nc_create(file_output,newRnValues)
ncvar_put(newRnNCDF4,"Rn",oldRnValues,start=c(1,1,1),count=c(raster.elevation@ncols,raster.elevation@nrows,1))
nc_close(newRnNCDF4)

proc.time()

###### Saving G ######

# Opening old G NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_G.nc", sep="")
nc <- nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

#New G file name
file_output<-paste(dados$Path.Output[1],"/",fic,"_G.nc",sep="")
oldGValues<-ncvar_get(nc,fic)
newGValues<-ncvar_def("G","daily",list(dimLonDef,dimLatDef,tdim),longname="G",missval=NaN,prec="double")
nc_close(nc)
newGNCDF4<-nc_create(file_output,newGValues)
ncvar_put(newGNCDF4,"G",oldGValues,start=c(1,1,1),count=c(raster.elevation@ncols,raster.elevation@nrows,1))
nc_close(newGNCDF4)

proc.time()

###### Saving Elevation raster ######
var_output <- paste(dados$Path.Output[1], "/elevation.tif", sep="")
writeRaster(raster.elevation, var_output, format="GTiff")

print("LABEL - WRITE DATA END")
proc.time()

