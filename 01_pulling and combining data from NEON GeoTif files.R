### Reading in GeoTIFF files
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### 
### Created July 16, 2021

### set up workspace
library(raster)
library(rgdal)
library(tidyverse)
library(data.table)
#library(feather)
#library(iotools)
#library(httr) # Package for pulling data from NEON website
#library(jsonlite) # Package for pulling data from NEON website

#file.choose()
rm(list=ls())
site_name <- "HARV"
setwd(paste0("C:\\Users\\kwilcox4\\Desktop\\NEON_NDVI\\", 
             site_name, 
             "\\NEON_indices-veg-spectrometer-mosaic\\")) ### This should be the folder one layer outside of the folder all the zipped folder are in
export_ndvi_path <- paste0("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\NewPatches\\data\\NDVI\\",
                           site_name)

# Pull folder names, dates of sampling
mother_folder_vector <- list.files()[grep(site_name, list.files())]
folder_info <- data.frame(full_name=mother_folder_vector) %>%
  mutate(full_name=as.character(full_name)) %>%
  separate(full_name,
           c("data_source", "domain", "site", "dpl", "prod_num", "rev_number", "year", "month", "pkg_type", "access_date", "release", "release_year")) %>%
  mutate(year_month = paste(year, month, sep="_"))
year_month_vector <- as.character(folder_info$year_month)

# Create folder for extracted data
dir.create("Extracted ndvi data\\")

### Loop through sample dates then through folders within each sample date

for(DATE in 1:length(year_month_vector)){
  # Create folder for year and month of data collection
  dir.create(paste0("Extracted ndvi data\\", year_month_vector[DATE], "\\"))
    
  # Create vectors of zipped folders to loop through and pull GeoTIF files from
  filename_vector <- list.files(mother_folder_vector[DATE])
  zipnames_vector <- filename_vector[grep(".zip", filename_vector)]
  
  print(year_month_vector[DATE])  
    ### Loop through zipped folders and pull NDVI data
    for(FOLDER in 1:length(zipnames_vector)){
      zip_path_temp <- paste0(mother_folder_vector[DATE], "\\", zipnames_vector[FOLDER]) # Select zip folder to extract from
      file_name_temp <- gsub("VegIndices.zip", "NDVI.tif", zipnames_vector[FOLDER]) ## Choose NDVI only (saves space)
                        
      extract_folder_path <- paste0("Extracted ndvi data\\", year_month_vector[DATE]) # Choose folder to save extracted file to
      
      unzip(zipfile=zip_path_temp, # Unzip it!
            exdir=extract_folder_path,
            files=file_name_temp)
      
      rm(zip_path_temp, file_name_temp) # clean up between loops
    }
    
    ### Loop to extract NDVI data and error from geotif files
    tif_mother_folder <- paste0("Extracted ndvi data\\", year_month_vector[DATE])
    all_tif_vector <- list.files(tif_mother_folder)
    ndvi_path_vector <- all_tif_vector[grep("_NDVI.tif", all_tif_vector)]
    ndvi_list <- list()
    
    ### prep press bar for FILE cycle
    pb = txtProgressBar(min = 0, max = length(ndvi_path_vector), initial = 0, style=3) 
    
    for(FILE in 1:length(ndvi_path_vector)){
      setTxtProgressBar(pb,FILE) # Update progress bar
      ndvi_path_temp <- paste0(tif_mother_folder, "\\", ndvi_path_vector[FILE])
    
      # NDVI data
      raster_temp <- raster(
        ndvi_path_temp
      )
      
      ndvi_df_temp <- 
        as.data.frame(raster_temp, xy=T) %>%
         rename(NDVI=gsub(".tif", "", ndvi_path_vector[FILE])) #%>%
      
      ndvi_list[[FILE]] <- ndvi_df_temp
      rm(ndvi_path_temp, raster_temp, ndvi_df_temp)
    #  print(paste0(FILE/length(ndvi_path_vector)*100 ,"% complete"))
      } ### End loop through ndvi.tif files

ndvi_df_master <- dplyr::bind_rows(ndvi_list) %>%
  mutate(site=site_name,
         year_month=year_month_vector[DATE]) %>%
  dplyr::select(site, year_month, x, y, NDVI)  

### Write ndvi data frame to file
fwrite(ndvi_df_master, 
       file=paste0(export_ndvi_path, "\\", site_name, "_ndvi_", year_month_vector[DATE], ".csv"),
       row.names=F)

# clean up within outer loop
rm(filename_vector, zipnames_vector, ndvi_list, ndvi_df_master) 

} ### End of loop through sampling dates


master_row_start <- (FILE-1)*1000000+1
master_row_end <- master_row_start+(nrow(ndvi_df_temp)-1)
master_row_nums <- (2-1)*1000000+1    
ndvi_df_master[master_row_start:master_row_end,] <- ndvi_df_temp
ndvi_df_master[51000001:52000000,]
rm(ndvi_path_temp, raster_temp, ndvi_df_temp)
### Checking for duplicate x,y coordinates (THIS WILL TAKE A LONG TIME)
# duplicated(ndvi_df_master[,1:3])
# length(unique(with(ndvi_df_master[,1:4], paste(x, y, sep="_")))) == nrow(ndvi_df_master)
# ### Other code to look at geotif files
#
# unzip("NEON.D06.KONZ.DP3.30026.001.2016-07.basic.20210715T221045Z.RELEASE-2021\\NEON_D06_KONZ_DP3_700000_4328000_VegIndices.zip",
#       exdir="NEON.D06.KONZ.DP3.30026.001.2016-07.basic.20210715T221045Z.RELEASE-2021\\NEON_D06_KONZ_DP3_700000_4328000_VegIndices")
# 
# file <- unz(unz("NEON_indices-veg-spectrometer-mosaic_KNZ.zip", 
#             filename="NEON_indices-veg-spectrometer-mosaic\\NEON.D06.KONZ.DP3.30026.001.2016-07.basic.20210715T221045Z.RELEASE-2021\\NEON_D06_KONZ_DP3_700000_4328000_VegIndices.zip"),
#             filename="NEON_D06_KONZ_DP3_700000_4328000_NDVI.tif")
# 
# file_temp <- unz("NEON.D06.KONZ.DP3.30026.001.2016-07.basic.20210715T221045Z.RELEASE-2021\\NEON_D06_KONZ_DP3_700000_4328000_VegIndices.zip", 
#                  filename="NEON_D06_KONZ_DP3_700000_4328000_NDVI.tif")
# ?unzip
# raster(file_temp)
# KNZ_NDVI_info <- capture.output(
#   GDALinfo("NEON.D06.KONZ.DP3.30026.001.2016-07.basic.20210715T221045Z.RELEASE-2021\\NEON_D06_KONZ_DP3_700000_4328000_VegIndices\\NEON_D06_KONZ_DP3_700000_4328000_NDVI.tif")
#   )
#   
# KNZ_NDVI <- raster("NEON.D06.KONZ.DP3.30026.001.2016-07.basic.20210715T221045Z.RELEASE-2021\\NEON_D06_KONZ_DP3_700000_4328000_VegIndices\\NEON_D06_KONZ_DP3_700000_4328000_NDVI.tif")
# KNZ_NDVI
# summary(KNZ_NDVI)
# summary(KNZ_NDVI, maxsamp=ncell(KNZ_NDVI))
# 
# KNZ_NDVI_df <- as.data.frame(KNZ_NDVI, xy=T) %>%
#   rename(NDVI=NEON_D06_KONZ_DP3_700000_4328000_NDVI)
# str(KNZ_NDVI_df)
ggplot() +
  geom_raster(data=ndvi_df_master, aes(x=x, y=y, fill=NDVI)) +
  scale_fill_viridis_c(na.value='deeppink') +
  coord_quickmap()


# 
# ggplot(data=KNZ_NDVI_df, aes(NDVI)) +
#   geom_histogram()
# crs(KNZ_NDVI)
# nlayers(KNZ_NDVI)



# ### Code to pull data directly from NEON website -- problem currently is that you can't extract data from zip files (at least I can't figure it out)
# 
# data_req <- GET("http://data.neonscience.org/api/v0/products/DP3.30026.001")
# data_req_text <- content(data_req, as="text")
# data_avail <- jsonlite::fromJSON(data_req_text, simplifyDataFrame = T, flatten=T)
# 
# print(data_avail)
# data_avail$data$productDescription
# data_avail$data$siteCodes$availableMonths[[24]]
# KONZ_urls <- unlist(data_avail$data$siteCodes$availableDataUrls[[24]])
# 
# KONZ_input <- GET(KONZ_urls[1])
# 
# KONZ_files <- jsonlite::fromJSON(content(KONZ_input, as="text"))
# 
# path <- KONZ_files$data$file$url[KONZ_files$data$files==KONZ_files$data$files[[134,1]]][1]
# path <- KONZ_files$data$file$url[KONZ_files$data$files==KONZ_files$data$files[[1,1]]][1]
# 
# 
# KONZ_files$data$files[1]
# KONZ_readme <- read.csv(path)
# 
# unzip(path, files="NEON_D06_KONZ_DP3_701000_4336000_NDVI.tif", exdir="C:\\Users\\kwilcox4\\Desktop")
# unzip(KONZ_files$data$files$name[1], files="NEON_D06_KONZ_DP3_701000_4336000_NDVI.tif", exdir="C:\\Users\\kwilcox4\\Desktop\\NEON_indices-veg-spectrometer-mosaic\\Extracted ndvi data\\KONZ")


# head(test)
# colnames(test) <- c("site", "year", "x", "y", "NDVI")
# with(test[1000:5000,], plot(NDVI))
# write_feather(test, )
# test <- read_feather("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\NewPatches\\data\\NDVI\\CPER\\2013\\CPER_ndvi_2013.csv")
# test <- read.csv.raw("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\NewPatches\\data\\NDVI\\CPER\\2013\\CPER_ndvi_2013.csv")
# test_fread <- fread("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\NewPatches\\data\\NDVI\\CPER\\2013\\CPER_ndvi_2013.csv")
# ?fwrite
# filter(test_feather, x==527250. & y==4524500.)
# head(test_feather)
# ?write_feather
# # ?read_feather

KNZ_NDVI_stack <- stack("NEON.D06.KONZ.DP3.30026.001.2016-07.basic.20210715T221045Z.RELEASE-2021\\NEON_D06_KONZ_DP3_700000_4328000_VegIndices\\NEON_D06_KONZ_DP3_700000_4328000_NDVI.tif")
KNZ_NDVI_stack
KNZ_NDVI_stack@layers



?GDALinfo
install.packages("gdalUtils")
library("gdalUtils")
?gd



install.packages("raster")
install.packages("rgdal")

