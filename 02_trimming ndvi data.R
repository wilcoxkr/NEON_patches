### cleaning and prepping NDVI data
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created July 21, 2021

### Set up workspace
setwd("C:\\Users\\kwilcox4\\Desktop\\NEON_NDVI\\")

library(sf)
library(sp)
library(raster)
library(tidyverse)
library(rgdal)
library(ggmap)
library(plotKML)
library(data.table)
#devtools::install_github("dkahle/ggmap", ref = "tidyup")
register_google(key = "AIzaSyAZWMJxgxkOrqx3XsxqTAAXu3OfBbVoI5o")

## Read in full shape file
sample_shape_full <- readOGR("Field_Sampling_Boundaries_2020\\terrestrialSamplingBoundaries.shp")
flight_shape_full <- readOGR("NEONflightpath\\AOP_flightboxesAllSites-polygon.shp")

summary(flight_shape_full)
## some summary stats on teh shape file
summary(sample_shape_full)
length(sample_shape_full)
head(sample_shape_full)

### HARV ###
## Subset Harvard Forest and turn into data frame
HARV_shape <- fortify(
  subset(sample_shape_full, siteID=="HARV")
)

HARV_location <- c(lon=mean(HARV_shape$long), 
                   lat=mean(HARV_shape$lat))

HARV_map <- get_googlemap(center=HARV_location, maptype="satellite", zoom=9)

ggmap(HARV_map) +
  geom_polygon(data=HARV_outline, aes(x=long, y=lat), fill="grey", alpha=0.2, col="yellow", size=0.2)


### SCBI ###
## Subset site and turn into data frame
SCBI_shape <- fortify(
  subset(sample_shape_full, siteID=="SCBI")
)
SCBI_flight_shape <- fortify(
  subset(flight_shape_full, siteID=="SCBI")
)

SCBI_location <- c(lon=mean(SCBI_flight_shape$long), 
                   lat=mean(SCBI_flight_shape$lat))

SCBI_map <- get_googlemap(center=SCBI_location, maptype="satellite", zoom=12)

ggmap(SCBI_map) +
  geom_polygon(data=subset(SCBI_shape,hole=FALSE), aes(x=long, y=lat, fill=piece), alpha=0.2, col="yellow", size=0.2)

ggmap(SCBI_map) +
  geom_polygon(data=SCBI_flight_shape, aes(x=long, y=lat), fill="grey", alpha=0.2, col="yellow", size=0.2)


### 
### Creating a shape file from flight data
###

scbi_ndvi <- fread("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\NewPatches\\data\\NDVI\\SCBI\\SCBI_ndvi_2016.csv")
str(scbi_ndvi)

scbi_spatial_file <- raster("C:\\Users\\kwilcox4\\Desktop\\NEON_NDVI\\SCBI\\raw\\NEON_indices-veg-spectrometer-mosaic\\NEON.D02.SCBI.DP3.30026.001.2016-07.basic.20210721T213848Z.RELEASE-2021\\NEON_D02_SCBI_DP3_743000_4307000_VegIndices\\NEON_D02_SCBI_DP3_743000_4307000_NDVI.tif")
scbi_crs <- st_crs(scbi_spatial_file) # Grabs coordinate reference system information from existing GEOTIFF at the site

scbi_aerial_locations <- st_as_sf(scbi_ndvi, coords = c("x", "y"), crs = scbi_crs) # converts file to shaper file






