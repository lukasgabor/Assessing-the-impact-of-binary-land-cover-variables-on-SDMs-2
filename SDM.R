# 3: Models

{library(maps)
  library(rgdal)
  library(raster)
  library(dismo)
  library(data.table)	
  library(scales)
  library(ranger)
  library(scam)
  library(fields)  
  library(PresenceAbsence) 
  library(pdp)
  library(tidyverse)
  library(gbm)
  library(sf)
  library(rgeos)
  library(plyr)
  library(SDMtune)
  library(ff)
  library(ffbase)
  library(foreach)
  library(parallel)
  library(doParallel)
  
  sessionInfo()
  capabilities() 
}

###############################
#MODELS (RF)
###############################
{ version <- "VS1"
where <- "hpc" # either local or hpc - REMEMBER TO TURN ON FOREACH
set.seed(1)
if(where=="local"){memory.limit(size=400000)}

# paths
if(where=="local"){
  path <- "~/ranges/"
}else{
  path <- "/gpfs/ysm/home/jc3893/30x30/"
}
setwd(path)

# Species names
list <- read_csv("species_list.csv")

# shapefile
countries <- readOGR("shapefiles", "countries")
countries_proj <- spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
aba <- countries[c(9,38),] %>%
  spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>%
  st_as_sf()

{# grid for thinning
  bb <- bbox(countries_proj)
  cs <- c(5000,5000)
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  grd
  sp_grd <- SpatialGridDataFrame(grd,
                                 data=data.frame(id=1:prod(cd)),
                                 proj4string=CRS(proj4string(countries_proj)))
}

# define season - change as needed
season <- "breeding"; months <- 6:8; mindate=152; maxdate=243; scode=2
# season <- "nonbreeding"; months <- c(12,1,2); mindate=335; maxdate=59; scode=3
# season <- "prebreeding_migration"; months <- 3:5; mindate=60; maxdate=151
# season <- "postbreeding_migration"; months <- 9:11; mindate=244; maxdate=334
# season <- "resident"; months <- 1:12; mindate=1; maxdate=365; scode=1


# define scale
scale <- 3 # options are 1, 3, 5, 10, 50, 100
plotres <- 500 # corresponding 500, 300, 200, 100, 50, 30

# ebird data
# ebird data
if (season=="resident"){
for (s in c("breeding","nonbreeding",
            "prebreeding_migration","postbreeding_migration")){
  ebird_ff <- read.csv.ffdf(file=paste0("30x30/ebird_",s,"_CEA.csv"),
                            colClasses=c(rep("factor", 2), rep("numeric",2),
                                         rep("integer",2), "numeric", "integer",
                                         "numeric","integer",rep("integer",677),
                                         rep("numeric",264)))
  assign(paste0("ebird_ff_",s), ebird_ff)
}
}else{
ebird_full <- read.csv.ffdf(file=paste0("ebird_",season, "_CEA.csv"),
                            colClasses=c(rep("factor", 2), rep("numeric",2),
                                         rep("integer",2), "numeric", "integer",
                                         "numeric","integer",rep("integer",677),
                                         rep("numeric",264)))
}

# prediction surface
pred_surface_full <- read.csv.ffdf(file=paste0("habitat_prediction-surface_",scale,".csv"),
                         colClasses=rep("numeric",46))

# standard grid
bbox <- readOGR("shapefiles", "30x30BBox_NAmerica") %>%
  spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>%
  st_as_sf()
grid <- raster(paste0("grids/Global_reference_raster_",scale,"km_CEA.tif")) %>%
  crop(bbox)

# Predictors (with weather data)
all_preds <- c(
  "year","day_of_year","time_observations_started","duration_minutes",
  "effort_distance_km","number_observers","bio1","bio12","bio15","cloudsd",
  "evisum","twi","tri","elev","pland_10_cropland_rainfed", 
  "pland_100_mosaic_tree_shrub","pland_11_cropland_rainfed", 
  "pland_110_mosaic_herbacious","pland_12_cropland_rainfed", 
  "pland_120_shrubland","pland_121_shrubland","pland_122_shrubland",
  "pland_130_grassland","pland_140_lichens_mosses","pland_150_sparse",
  "pland_152_sparse","pland_153_sparse","pland_160_flooded_freshwater",
  "pland_170_flooded_saltwater","pland_180_flooded_shrub", 
  "pland_190_urban","pland_20_cropland_irrigated","pland_200_barren",     
  "pland_201_barren","pland_202_barren","pland_210_water",
  "pland_220_ice","pland_30_mosaic_cropland","pland_40_mosaic_natural_veg",
  "pland_50_evergreen_broadleaf","pland_60_deciduous_broadleaf", 
  "pland_61_deciduous_broadleaf","pland_62_deciduous_broadleaf", 
  "pland_70_evergreen_needleleaf","pland_71_evergreen_needleleaf",
  "pland_72_evergreen_needleleaf","pland_80_deciduous_needleleaf",
  "pland_81_deciduous_needleleaf","pland_82_deciduous_needleleaf",
  "pland_90_mixed_forest")
all_preds_sc <- c(all_preds[1:6], paste0(all_preds[7:50], "_", scale))

# mol range
mol.range <- readOGR("range_polygons", "jetz_ranges_sub")

# directories
if(where=="hpc"){
  dir <- paste0("/gpfs/loomis/pi/jetz/data/results_30x30/sdm_birds/",
                    version,"_30x30/")
}else{
  dir <- paste0("runs/",version,"_30x30/")
}
dir.create(dir)
dir.create(paste0(dir, "plots/"))
dir.create(paste0(dir, "plots/", scale, "/"))
setwd(dir)

if(where=="hpc"){
# parallel computing
#n.cores <- parallel::detectCores()
n.cores <- 4
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
}

}

# ----------------------------------------------------------------------
# RANGER
# ----------------------------------------------------------------------		
# find a species
# which(list[[1]]=="Painted Bunting")
{# loop species
 #series <- 1:10
 series <- 1:nrow(list)
  # series <- 451
  #for (i in series){ 
    foreach (i=series) %dopar% { # TURN ON FOR HPC
    try({
      {library(maps)
        library(rgdal)
        library(raster)
        library(dismo)
        library(data.table)	
        library(scales)
        library(ranger)
        library(scam)
        library(fields)  
        library(PresenceAbsence) 
        library(pdp)
        library(tidyverse)
        library(gbm)
        library(sf)
        library(rgeos)
        library(plyr)
        library(SDMtune)
        library(ff)
        library(ffbase)
        library(foreach)
        library(parallel)
        library(doParallel)
      }
      # cleanup
      rm(sprange)
      gc()
      # Establish names
      common <- as.character(list[i,1])
      sci.name <- as.character(list[i,2])
      sci_name <- gsub(" ","_",sci.name)
      code <- as.character(list[i,5]) 
      sp_status <- ifelse(list[i,13]=="y","waterbird","landbird")
      print(c(i, common, sci_name))
      # check if species has been done
      if(!file.exists(paste0(sci_name,"_complete_",season,"_",scale,".csv"))){
        rm(sprange)
        # Expert range = modeling domain
        setwd(path)
        try({
          total.range <- readOGR(unzip(
            paste0("range_polygons/",code,"-range-2020.gpkg.zip"),
            paste0(code,"-range-mr-2020.gpkg")), "range")
        # Get seasonal range of species
          sprange <- total.range[total.range$season_name==season,]
        })
        if(!exists("sprange")){
          try({
            total.range <- readOGR(unzip(
              paste0("/gpfs/ysm/home/jc3893/30x30/range_polygons/",
                     code,"-range-2021.gpkg.zip"),
              paste0(code,"-range-mr-2021.gpkg")), "range")
            # Get seasonal range of species
            sprange <- total.range[total.range$season==season,]
          })
        }
        if(!exists("sprange")){
          sprange <- mol.range[mol.range$sciname==sci.name & 
                                mol.range$season==scode,]
          sprange$season <- season
        }
        if(nrow(sprange)==0){STOP} # wrong season for this sp or no range exists
        setwd(dir)
          # Proj range, buffer, back to unproj
          sprange <- spTransform(sprange, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
          # if(area(sprange)>200000){
          #   sprange <- buffer(sprange, -5000) # remove small bits
          # }
          sprange <- buffer(sprange, 200000)
          box <- bbox(sprange)
          if(where=="local"){plot(sprange, col="bisque")}
          # Back to SP dataframe
          range.df <- data.frame(ID=1, row.names="buffer")
          sprange <- SpatialPolygonsDataFrame(sprange, range.df, match.ID=F)
          # ST range
          rangest <- st_as_sf(sprange)
          # check that range is at least partly in the ABA
          sf::sf_use_s2(FALSE)
          intersect <- st_intersects(aba, rangest, sparse=T)
          if(sum(lengths(intersect))>0){
          print(paste("START", i, common, season))
          # directories
          dir.create(paste0(sci_name,"/"))
          dir.create(paste0(sci_name,"/predictions/"))
          dir.create(paste0(sci_name,"/predictions/ROR/"))
          dir.create(paste0(sci_name,"/predictions/PA/"))
          dir.create(paste0(sci_name,"/points/"))
          dir.create(paste0(sci_name,"/accuracy/"))
          dir.create(paste0(sci_name,"/rangemap/"))
          # Save range
          writeOGR(sprange, paste0(sci_name,"/rangemap"),
                   paste0(sci_name,"_",season), driver="ESRI Shapefile",
                   overwrite_layer=T)
          
          # ebird - get cropped rows, columns needed
          # different path for residents vs migrants
          if (season=="resident"){
            for (sson in c("breeding","nonbreeding",
                        "prebreeding_migration","postbreeding_migration")){
              ebird_full <- get(paste0("ebird_ff_",sson))
              ebird_index <- ffwhich(ebird_full, 
                                     y<box[2,2] &
                                       y>box[2,1] &
                                       x<box[1,2] &
                                       x>box[1,1])
              ebird_full <- data.frame(ebird_full[ebird_index,
                                                  c(gsub(" ",".",sci.name),all_preds,
                                                    "y","x","checklist_id")])
              # NA removal (only a few rows)
              ebird_full <- ebird_full[complete.cases(ebird_full),]
              
              # Cut ebird to range
              print("processing ebird data")
              ebird_full <- st_as_sf(ebird_full, coords = c("x", "y"), 
                                     crs = crs(rangest))
              ebird.int <- st_intersects(ebird_full, rangest, sparse=F)
              ebird_full <- ebird_full[as.vector(ebird.int),]
              #plot(ebird[sample.int(nrow(ebird), nrow(ebird)/20),], max.plot=1)
              ebird_full <- data.frame(sf:::as_Spatial(ebird_full))
              ebird_full$optional <- NULL
              names(ebird_full)[(ncol(ebird_full)-1):ncol(ebird_full)] <- 
                c("x", "y")
              rm(ebird.int)
              
              # Impose data size limit
              if(nrow(ebird_full) > 1000000){
                ebird_full <- ebird_full[sample(nrow(ebird_full),1000000),]
              }
              
              # add to other seasons
              if(sson=="breeding"){ebird <- ebird_full}else{
                ebird <- rbind(ebird, ebird_full)
              }
            }
            
            # randomize rows
            ebird <- ebird[sample(nrow(ebird)),]
            
          }else{ # migrants
          # shorebirds- exclude August for breeding models
          if (season=="breeding" & list[i,10]=="Charadriiformes"){
            ebird_index <- ffwhich(ebird_full, 
                                   y<box[2,2] &
                                     y>box[2,1] &
                                     x<box[1,2] &
                                     x>box[1,1] &
                                     day_of_year<213)
          }else{
            ebird_index <- ffwhich(ebird_full, 
                                   y<box[2,2] &
                                     y>box[2,1] &
                                     x<box[1,2] &
                                     x>box[1,1])
          }
          ebird <- data.frame(ebird_full[ebird_index,
                              c(gsub(" ",".",sci.name),all_preds_sc,
                             "y","x","checklist_id")])
          
          # Remove scale from names
          colnames(ebird) <- gsub(paste0("_",scale,"$"),"",colnames(ebird))
          
          # NA removal (only a few rows)
          ebird <- ebird[complete.cases(ebird),]
          
          # Cut ebird to range
          print("processing ebird data")
          ebird <- st_as_sf(ebird, coords = c("x", "y"), 
                            crs = crs(rangest))
          ebird.int <- st_intersects(ebird, rangest, sparse=F)
          ebird <- ebird[as.vector(ebird.int),]
          if(where=="local"){
            plot(ebird[sample.int(nrow(ebird), nrow(ebird)/20),], max.plot=1)}
          ebird <- data.frame(sf:::as_Spatial(ebird))
          ebird$optional <- NULL
          names(ebird)[(ncol(ebird)-1):ncol(ebird)] <- c("x", "y")
          rm(ebird.int)
          }
          
          # Check if there's enough data
          if(nrow(ebird)<50){
            writeLines("stop", paste0(sci_name,"_",common,"_failed_too_few_points_",scale,".txt"))
            STOP
          }
          
          # Column for presence/absence
          ebird$species_observed <- ebird[,gsub(" ",".",sci.name)]
          
          # Cut pred surface to range
          print("processing prediction surface")
          surface_index <- ffwhich(pred_surface_full, 
                                   y<box[2,2] &
                                     y>box[2,1] &
                                     x<box[1,2] &
                                     x>box[1,1])
          pred_surface <- data.frame(pred_surface_full[surface_index,])
          pred_surface <- st_as_sf(pred_surface, coords = c("x", "y"), 
                                   crs = crs(sprange))
          surface.int <- st_intersects(pred_surface, rangest, sparse=F)
          pred_surface <- pred_surface[as.vector(surface.int),]
          pred_surface <- data.frame(sf:::as_Spatial(pred_surface))
          pred_surface$optional <- NULL
          names(pred_surface)[(ncol(pred_surface)-1):ncol(pred_surface)] <- c("x", "y")
          pred_surface <- pred_surface[complete.cases(pred_surface),]
          rm(surface.int)
          gc()
          
          # Check if there's enough surface
          if(nrow(pred_surface)<10){
            writeLines("stop", paste0(sci_name,"_",common,"_failed_too_few_cells_",scale,".txt"))
            STOP
          }
          
          # Populate with values
          pred_surface$duration_minutes <- 60
          pred_surface$effort_distance_km <- 1
          pred_surface$number_observers <- 1
          pred_surface$year <- 2020
          # Find max hr of day for spatial predictions
          hrs <- ebird[,c("time_observations_started","species_observed")] %>%
            group_by(round(time_observations_started)) %>%
            summarize_at(vars(species_observed), mean)
          pred_surface$time_observations_started <- as.numeric(
            hrs[which.max(hrs$species_observed),1])
          # set day for spatial predictions
          if(season=="nonbreeding"){
            pred_surface$day_of_year <- sample(
              c(1:maxdate,mindate:365), nrow(pred_surface), replace=T)
          }else{
          pred_surface$day_of_year <- sample(
            mindate:maxdate, nrow(pred_surface), replace=T)
          }
          
          # Splits
          print("splitting and filtering data")
          {#training and testing samples
            index <- runif(nrow(ebird))
            train <- ebird[which(index>.4),]
            test <- ebird[which(index<=.4 & index>.2),]
            train.oob <- ebird[which(index<=.2),]
            rm(index)
          }
          
          {# Spatiotemporal thinning
            train$id <- NULL
            coordinates(train) <- c("x", "y")
            crs(train) <- crs(sp_grd)
            over <- over(train, sp_grd)
            train <- cbind(data.frame(train), over)
            # get cell id and week number for each checklist
            train$optional <- NULL
            checklist_cell <- train %>% 
              mutate(cell = id,
                     week = ceiling(day_of_year/7))
            # sample one checklist per grid cell per week
            # sample detection/non-detection independently 
            train <- checklist_cell %>% 
              group_by(species_observed, year, week, cell) %>% 
              sample_n(size = 1) %>% 
              ungroup()
          }
          {# Assemble TrainOOB
            train.oob$id <- NULL
            coordinates(train.oob) <- c("x","y")
            crs(train.oob) <- crs(sp_grd)
            over <- over(train.oob, sp_grd)
            train.oob <- cbind(data.frame(train.oob), over)
            train$optional <- NULL
            checklist_cell_oob <- train.oob %>% 
              mutate(cell = id,
                     week = ceiling(day_of_year/7))
            train.oob <- checklist_cell_oob %>% 
              group_by(species_observed, year, week, cell) %>% 
              sample_n(size = 1) %>% 
              ungroup()
          }
          
          {# remove oversampling for occurrence model, limit size
          train <- train[!duplicated(train$checklist_id), ]
          # Binary Occurrence Response, coded as numeric
          train$pres_abs <- as.numeric(train[[which(names(train) == "species_observed")]] > 0)
          # balanced sampling
          pos_fraction <- mean(as.numeric(train$pres_abs))
          # Check that positives ARE the minority class
          if (pos_fraction > 0.5) pos_fraction <- 1 - pos_fraction
          # Ranger binary response model requires code response as factor
          train$pres_abs <- as.factor(as.numeric(train$pres_abs))
          }
          
          # Save points alone (only once)
          if(scale==1){
          write_csv(train[,c("y","x","pres_abs")], 
                    paste0(sci_name,"/points/pts_binary_occurrences_",
                           sci_name,"_",season,"_00.csv"))
          }
          
          # For water birds, loop continuous or binary landcover
          # Multiple thresholds of defining binary values
          if(sp_status=="landbird"){scenarios="cont"}else{
          scenarios=c("cont","binary.01","binary.10","binary.20","binary.50")}
          for (water in scenarios){
            try({
            if(water=="binary.01"){
            train_pland <- train$pland_210_water
            test_pland <- test$pland_210_water
            oob_pland <- train.oob$pland_210_water
            pred_pland <- pred_surface$pland_210_water}
              if(grepl("binary", water, fixed=T)){
                thresh <- as.numeric(str_sub(water, 7, 9))
                train$pland_210_water <- ifelse(train_pland>thresh,1,0)
                test$pland_210_water <- ifelse(test_pland>thresh,1,0)
                train.oob$pland_210_water <- ifelse(oob_pland>thresh,1,0)
                pred_surface$pland_210_water <- ifelse(pred_pland>thresh,1,0)
             }
            print(paste("model", water))
            
            {#MODEL
              # Model Formula, threads
              m_formula <- paste("pres_abs ~", paste(all_preds, collapse = "+"))
              n_threads <- 7
              # Balanced Ranger Occurrence Model
              rf_occ <- ranger::ranger(
                formula =  m_formula, 
                num.trees = 100, 
                #max.depth = 0, 
                importance = "impurity",
                num.threads = n_threads,
                respect.unordered.factors = "order",
                always.split.variables = NULL,
                probability = TRUE,
                replace = TRUE, 
                sample.fraction = c(pos_fraction, pos_fraction),
                data = train)
            }
            
            # Test Set PPMS
            test$occ_pred <- predict(
              rf_occ, data = test, 
              type = "response",
              num.threads = n_threads)$predictions[,2]
            
            {#Select Optimal Threshold on train OOB data
              trainOOB.pred <- predict(rf_occ, data=train.oob, type="response")
              pa.df <- data.frame("space.holder.tag",
                                  obs = train.oob$species_observed, 
                                  ppp = trainOOB.pred$predictions[,2] )
              pa.metrics <- presence.absence.accuracy(
                pa.df,
                threshold = quantile(
                  pa.df$ppp,
                  probs = seq(from=0, to=1, length=1000), na.rm =T),
                na.rm = T, st.dev = F)
              pa.metrics$PCC <- pa.metrics$Kappa <- NULL
              pa.metrics$TSS <- (pa.metrics$sensitivity + pa.metrics$specificity) - 1
              
              # Save confusion matrix
              write_csv(pa.metrics, 
                        paste0(sci_name,"/accuracy/accuracy_matrix_",
                               sci_name,"_",season,"_",scale,"_",water,"_00.csv"))
              
              # optimal_thresh_position <- which.max(pa.metrics$Kappa) 
              optimal_thresh_position <- which.max(pa.metrics$sensitivity +
                                                     pa.metrics$specificity)
            }
            {# Compute Test set PPMs
              test.pred <- predict(rf_occ, data=test, type="response")
              pa.df <- data.frame("space.holder.tag",
                                  obs = test$species_observed, 
                                  ppp = test.pred$predictions[,2] )
              pa.metrics <- presence.absence.accuracy(
                pa.df,
                threshold = pa.metrics$threshold[ optimal_thresh_position ] ,
                na.rm = T, st.dev = F)
              thresh <- pa.metrics$threshold
              pa.metrics["PCC"] <- pa.metrics["Kappa"] <- NULL
              pa.metrics["TSS"] <- (pa.metrics$sensitivity + pa.metrics$specificity) - 1
            }
            
            # Spatial predictions
            pred_surface$occ_pred <- predict(
              rf_occ, data = pred_surface, type = "response",
              num.threads = n_threads)$predictions[,2]
            pred_surface$occ_thresh <- ifelse(
              pred_surface$occ_pred>pa.metrics$threshold,1,0)
            
            # Save rasterized predictions
            coordinates(pred_surface) <- c("x","y")
            ras <- rasterize(pred_surface, grid, pred_surface$occ_pred, fun='first')
            writeRaster(ras, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                    sci_name,"_",season,"_",scale,"_",water,"_00.tif"),
                        overwrite=T)
            pred_ones <- pred_surface[pred_surface$occ_thresh==1,]
            ras <- rasterize(pred_ones, grid, pred_ones$occ_thresh, fun='first')
            writeRaster(ras, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                    sci_name,"_",season,"_",scale,"_",water,"_00.tif"),
                        overwrite=T)
            pred_surface <- data.frame(pred_surface)
            pred_ones <- data.frame(pred_ones)
            
            occ.probs <- pred_surface[,c("x","y","occ_pred")]
            write_csv(occ.probs, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                        sci_name,"_",season,"_",scale,"_",water,"_00.csv"))
            occ.probs <- pred_ones[,c("x","y","occ_thresh")]
            write_csv(occ.probs, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                        sci_name,"_",season,"_",scale,"_",water,"_00.csv"))
            rm(pred_ones)
            
            # Predictions projected
              coordinates(pred_surface) <- c("x","y")
              pred_surface <- data.frame(pred_surface)
            pred_surface$occ_pred <- predict(
              rf_occ, data = pred_surface, type = "response",
              num.threads = n_threads)$predictions[,2]
            pred_surface$occ_thresh <- ifelse(
              pred_surface$occ_pred>pa.metrics$threshold,1,0)
            
            {# Plot maps
              jpeg(paste0("plots/",scale,"/PO_spatial_predictions_",
                          sci_name,"_",season,"_",scale,"_",water,"_00.jpeg"),
                   width=500, height=500, units="px")
              # Plot prob of occurrence
              quilt.plot(
                x = pred_surface$x,
                y = pred_surface$y,
                z = pred_surface$occ_pred,
                nrow = plotres,
                ncol = plotres,
                na.rm = T,
                main = paste(common,"-",sci.name,"\n",season,"-",scale))
              plot(countries_proj, col="transparent", border="black",
                   usePolypath=F, add=T)
              dev.off()
              # Plot thresholded occurrence
              jpeg(paste0("plots/",scale,"/Thresholded_spatial_predictions_",
                          sci_name,"_",season,"_",scale,"_",water,"_00.jpeg"),
                   width=500, height=500, units="px")
              quilt.plot(
                x = pred_surface$x,
                y = pred_surface$y,
                z = pred_surface$occ_thresh,
                nrow = plotres,
                ncol = plotres,
                na.rm = T,
                add.legend = F,
                col = c("transparent", "darkgreen"),
                main = paste(common,"-",sci.name,"\n",season,"-",scale))
              plot(countries_proj, col="transparent", border="black", 
                   usePolypath=F, add=T)
              dev.off()
              }
            
            # Species paths
            modpath <- paste0("/mnt/Work/Yale/MOL/results_30x30/sdm_birds/",
                           version,"_30x30/",sci_name)
            modName <- paste0(sci_name,"_",season,"_",scale,"_",water)
            modURL <- modPDF <- bad_ok_good <- comments <- 
              rangeOffset <- elevOffset <- deltaAUC <- NA
            modPathROR <- paste0(modpath,"/predictions/ROR/PO_spatial_predictions_",
                                 sci_name,"_",season,"_",scale,"_",water,"_00.jpeg")
            modPathPA <- paste0(modpath,"/predictions/PA/Thresholded_spatial_predictions_",
                                sci_name,"_",season,"_",scale,"_",water,"_00.jpeg")
            rangePath <- paste0(modpath,"/rangemap/",sci_name,"_",season,".shp")
            ptsPath <- paste0(modpath,"/points/pts_occurrences_",season,".csv")
            ptsBgPath <- NA
            confMatPath <- paste0(modpath,"/accuracy/accuracy_matrix_",
                                  sci_name,"_",season,"_",scale,"_",water,"_00.csv")
            POpredsPath <- paste0(modpath,"/predictions/ROR/PO_predictions_",
                                  sci_name,"_",season,"_",scale,"_",water,"_00.csv")
            ThreshpredsPath <- paste0(modpath,"/predictions/PA/Thresholded_predictions_",
                                      sci_name,"_",season,"_",scale,"_",water,"_00.csv")
            POpredsrasPath <- paste0(modpath,"/predictions/ROR/PO_predictions_",
                                  sci_name,"_",season,"_",scale,"_",water,"_00.tif")
            ThreshpredsrasPath <- paste0(modpath,"/predictions/PA/Thresholded_predictions_",
                                      sci_name,"_",season,"_",scale,"_",water,"_00.tif")
            envVars <- RDataPath <- rorOrigPath <- paOrigPath <- rangeOrigPath <-
              ptsBgOrigPath <- statsOrigPath <- NA
            
            # Species values
            imp.scores <- rf_occ$variable.importance
            model.output <- unlist(c(common, sci_name, sp_status, season, scale, water, version, 
                                     modName, modURL, modPDF, bad_ok_good, 
                                     comments, rangeOffset, elevOffset,
                                     nrow(train), nrow(occ.probs)*(scale^2),deltaAUC, 
                                     pa.metrics[5],pa.metrics[c(2,6,3,4)],
                                     modPathROR,modPathPA,rangePath,ptsPath,
                                     ptsBgPath,confMatPath,POpredsPath,ThreshpredsPath,
                                     POpredsrasPath,ThreshpredsrasPath,
                                     envVars,imp.scores,RDataPath,rorOrigPath,paOrigPath,
                                     rangeOrigPath,ptsBgOrigPath,statsOrigPath))
            names(model.output)<- c("common_name","sciname", "status",
                                    "season","scale","water_cover","version","modName",
                                    "modURL","modPDF","bad_ok_good","comments",
                                    "rangeOffset","elevOffset","noPts","range_area",
                                    "deltaAUC","AUC","SPSthreshold","TSS",
                                    "Sensitivity","Specificity",
                                    "modPathROR","modPathPA",
                                    "rangePath","ptsPath","ptsBgPath","confMatPath",
                                    "POPredsPath","threshPredsPath","POPredsRasPath",
                                    "threshPredsRasPath","envVars",
                                    paste0(all_preds,"_imp"),
                                    "RDataPath","rorOrigPath","paOrigPath",	
                                    "rangeOrigPath","ptsBgOrigPath","statsOrigPath")            
            if(sp_status=="landbird" | water=="cont"){
              sp.output <- model.output}
            if(sp_status=="waterbird" & water!="cont"){
              sp.output <- rbind(sp.output, model.output)}
            })} #end cont/binary water cover loop
          # get delta AUC
          if(is.matrix(sp.output)){
            sp.output <- data.frame(sp.output)
          AUCs <- as.numeric(as.character(sp.output$AUC))
          deltaAUC <- round(max(AUCs) - AUCs, 4)
          sp.output$deltaAUC <- deltaAUC}
          if(ncol(data.frame(sp.output))==1){
            sp.output <- t(data.frame(sp.output))
            rownames(sp.output) <- NULL
          }
          write_csv(data.frame(sp.output), 
                    paste0(sci_name,"_complete_",season,"_",scale,".csv"))
        } # end range exists in bbox check
      } #end species completion check
  }) # end try
} #end sp
} #end models
  