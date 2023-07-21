{### redoing predictor extraction
  library(lubridate)
  library(sf)
  library(rgdal)
  library(gridExtra)
  library(tidyverse)
  library(raster)
  library(exactextractr)
  library(viridis)
  library(ncdf4)
  #library(plyr)
  
  memory.limit(300000)
  
  # resolve namespace conflicts
  select <- dplyr::select
  map <- purrr::map
  projection <- raster::projection
}

# Climate/Topography - checklists
#standard vars
elev <- raster("gmted/elevation_1KMmd_GMTEDmd.tif")
elev <- round(elev, 1)
data <- read_csv("~/data/surface/habitat_prediction-surface-new.csv")
data <- data[,1:20]
latlon <- data[,c("latitude","longitude")]
data <- st_as_sf(data, coords=c("longitude","latitude"))
elev <- extract(elev, data)
data <- data.frame(cbind(latlon, data, elev))
data$geometry <- NULL
write_csv(data, "~/data/surface/habitat_prediction-surface-new.csv")

bio1 <- raster("chelsa/CHELSA_bio_1_NA30x30_cropped.tif")
bio12 <- raster("chelsa/CHELSA_bio_12_NA30x30_cropped.tif")
bio15 <- raster("chelsa/CHELSA_bio_15_NA30x30_cropped.tif")
cloudsd <- raster("gmted/MODCF_intraannualSD_resampled_masked_NA30x30_cropped.tif")
evisum <- raster("gmted/EVI_Summer_6to8_resampled_NA30x30_cropped.tif")
twi <- raster("gmted/TWI_resampled_NA30x30_cropped.tif")
tri <- raster("gmted/tri_1KMmd_GMTEDmd_resampled_masked_NA30x30_cropped.tif")
charlievars <- stack(bio1, bio12, bio15, cloudsd,
                     evisum, twi, tri)
charlievars <- round(charlievars, 1)
names(charlievars) <- c("bio1","bio12","bio15","cloudsd",
                        "evisum","twi","tri")
data <- read_csv("~/data/surface/habitat_prediction-surface-new.csv")
data <- st_as_sf(data, coords=c("longitude","latitude"))
for (i in 1:7){
  print(i)
  pred <- extract(charlievars[[i]], data)
  if(i==1){allpred <- pred}else{allpred <- cbind(allpred, pred)}
  write_csv(data.frame(cbind(latlon, data, allpred)), 
            "~/data/surface/habitat_prediction-surface-new.csv")
  gc()
}

# reorganize, ditch old lc
data <- read_csv("~/data/surface/habitat_prediction-surface-new.csv")
data <- data[,c(1:2,21,23:29)]
colnames(data)[3:10] <- c("elev","bio1","bio12","bio15","cloudsd",
                          "evisum","twi","tri")
write_csv(data, "~/data/surface/habitat_prediction-surface-new.csv")

# Percent Landcover (lab version)
stack <- list.files("~/data/esa_cci/", full.names=T) %>%
  stack()
stack <- round(stack, 1)

data <- st_as_sf(data, coords=c("longitude","latitude"))
for (i in 1:nlayers(stack)){
  print(i)
  pred <- extract(stack[[i]], data)
  if(i==1){allpred <- pred}else{allpred <- cbind(allpred, pred)}
  write_csv(data.frame(cbind(data, allpred)), 
            "~/data/surface/habitat_prediction-surface-lc.csv")
  gc()
}

# colnames for landcover
data <- data.frame(cbind(latlon, data, allpred))
data$geometry <- NULL
lc_names <- c("pland_10_cropland_rainfed", 
              "pland_100_mosaic_tree_shrub", 
              "pland_11_cropland_rainfed", 
              "pland_110_mosaic_herbacious", 
              "pland_12_cropland_rainfed", 
              "pland_120_shrubland",
              "pland_121_shrubland",
              "pland_122_shrubland",
              "pland_130_grassland", 
              "pland_140_lichens_mosses",
              "pland_150_sparse",
              "pland_152_sparse",
              "pland_153_sparse",
              "pland_160_flooded_freshwater",
              "pland_170_flooded_saltwater",
              "pland_180_flooded_shrub", 
              "pland_190_urban", 
              "pland_20_cropland_irrigated",       
              "pland_200_barren",     
              "pland_201_barren",     
              "pland_202_barren",
              "pland_210_water",
              "pland_220_ice",
              "pland_30_mosaic_cropland",
              "pland_40_mosaic_natural_veg",
              "pland_50_evergreen_broadleaf", 
              "pland_60_deciduous_broadleaf", 
              "pland_61_deciduous_broadleaf", 
              "pland_62_deciduous_broadleaf", 
              "pland_70_evergreen_needleleaf",
              "pland_71_evergreen_needleleaf",
              "pland_72_evergreen_needleleaf", 
              "pland_80_deciduous_needleleaf",
              "pland_81_deciduous_needleleaf",
              "pland_82_deciduous_needleleaf",
              "pland_90_mixed_forest")
colnames(data)[11:ncol(data)] <- lc_names
write_csv(data, "~/data/surface/habitat_prediction-surface-lc.csv")

# New landcover for ebird data
for (season in c(
  # "nonbreeding",
                # "breeding",
                 "prebreeding_migration"
                 # "postbreeding_migration"
                 )){
  print(season)
  ebird <- read_csv(paste0("~/data/processed/ebird_",season,"_habitat.csv"))
  ebird <- ebird[,c("latitude","longitude","checklist_id")]
  latlon <- ebird[,c("latitude","longitude")]
  gc()
  ebird <- st_as_sf(ebird, coords=c("longitude","latitude"))
  pred <- extract(stack, ebird)
  colnames(pred) <- lc_names
  ebird <- read_csv(paste0("~/data/processed/ebird_",season,"_habitat.csv"))
  ebird <- ebird[,c(1:687,706:713)]
  ebird <- cbind(ebird, pred)
  write_csv(ebird, paste0("~/data/processed/ebird_",season,"_new.csv"))
  gc()
}
  
  
  
  
  
  

