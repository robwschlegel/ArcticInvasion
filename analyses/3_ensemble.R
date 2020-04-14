# 3_ensemble.R
# The purose of this script is to take all of the results from Bbiomod2
# and MaxEnt and merge them together into one final binary map for each species.
# 1: Setup the environment
# 2: Functions
# 3: Calculate ensemble models

# 1: Setup ----------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))
library(tidyverse)
library(biomod2)
library(sp)
library(dtplyr)
library(doParallel); registerDoParallel(cores = 50)

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# The species occurrence data
sps_names <- str_remove(dir("data/occurrence", full.names = F), pattern = "_near.csv")

# Choose species
# sps <- sps_names[13]


# 2: Functions ------------------------------------------------------------

# convenience function for loading MaxEnt data
load_maxent <- function(maxent_file){
  res <- as.data.frame(raster(maxent_file), xy = TRUE) %>%
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(model = "MaxEnt") %>% 
    dplyr::select(model, x, y, z)
}

# Function that converts a raster layer to a binary dataframe
raster_to_df_binary <- function(layer_num, proj_x, res_table){
  
  # Convert to Dataframe
  proj_sub_df <- as.data.frame(proj_x@layers[[layer_num]], xy = T)
  
  # Find TSS cutoff
  sub_cutoff <- strsplit(colnames(proj_sub_df[3]), "_")[[1]]
  res_table_sub <- res_table %>% 
    filter(run == sub_cutoff[3], model == sub_cutoff[4])
  
  # Stop if specificity is 0
  if(res_table_sub$Specificity == 0) return()
  
  # Stop run if the model TSS is not >=0.7
  if(res_table_sub$score[1] < 0.7) return()
  
  # Convert to binary
  proj_sub_df_binary <- proj_sub_df %>% 
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(z = if_else(z >= res_table_sub$Cutoff[1], 1, 0),
           run = sub_cutoff[3],
           model = sub_cutoff[4]) %>% 
    dplyr::select(model, run, x, y, z)
}

# Function that creates the final mean binary value for a projection
proj_binary_mean <- function(proj_biomod, proj_maxent, res_table, proj_name){
  
  # Extract all of the models as a long data.frame
  df_biomod <- plyr::ldply(1:length(proj_biomod@layers), raster_to_df_binary, 
                           .parallel = F, proj_biomod, res_table)
  
  # Create rounded mean binary values per model
  # system.time(
  df_biomod_mean <- lazy_dt(df_biomod) %>% 
    dplyr::select(-run) %>% 
    group_by(model, x, y) %>% 
    summarise(z = round(mean(z))) %>% 
    ungroup() %>% 
    data.frame()
  # ) # 340 seconds
  saveRDS(df_biomod_mean, paste0("data/biomod/",res_table$sps[1],"_df_",proj_name,".Rds"))
  
  # Add MaxEnt binary and calculate final binary result
  df_res <- rbind(df_biomod_mean, proj_maxent)
  # system.time(
  df_res <- lazy_dt(df_res) %>% 
    group_by(x, y) %>% 
    summarise(z = round(mean(z))) %>% 
    ungroup() %>% 
    data.frame()
  # ) # 63 seconds
  rm(df_biomod_mean); gc()
  return(df_res)
}

# Function for converting dataframe to raster
df_to_raster <- function(df_proj){
  
  # coerce to SpatialPointsDataFrame
  spg <- df_proj
  coordinates(spg) <- ~ x + y
  
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  
  # coerce to raster
  raster_df <- raster(spg)
}

# Function that crawls through species files, finds cutoffs, and creates binary output per model
sps_binary <- function(sps){
  
  print(paste0("Began run on ",sps," at ", Sys.time()))
  
  # The cutoff table
  res_table <- plyr::adply(get_evaluations(loadRData(paste0(sps,"/",sps,".",sps,".models.out"))), c(1,3,4,5)) %>% 
    dplyr::rename(test = X1, model = X2, run = X3, PA = X4, score = Testing.data) %>% 
    mutate(sps = sps) %>% 
    dplyr::select(sps, everything()) %>% 
    filter(test == "TSS")
  
  # Present projections
  print(paste0("Present projections for ",sps," at ",Sys.time()))
  proj_present <- proj_binary_mean(get_predictions(loadRData(paste0(sps,"/proj_present/",sps,".present.projection.out"))),
                                   load_maxent(paste0("data/maxent/",sps,"_avg_binary.tif")),
                                   res_table, "present")
  writeRaster(df_to_raster(proj_present), paste0("data/biomod/",sps,"_binary_present.asc"), overwrite = TRUE)
  rm(proj_present); gc()
  
  # 2050 projections
  print(paste0("2050 projections for ",sps," at ",Sys.time()))
  proj_2050 <- proj_binary_mean(get_predictions(loadRData(paste0(sps,"/proj_2050/",sps,".2050.projection.out"))),
                                load_maxent(paste0("data/maxent/",sps,"_2050_45_avg_binary.tif")),
                                res_table, "2050")
  writeRaster(df_to_raster(proj_2050), paste0("data/biomod/",sps,"_binary_2050.asc"), overwrite = TRUE)
  rm(proj_2050); gc()
  
  # 2100 projections
  print(paste0("2100 projections for ",sps," at ",Sys.time()))
  proj_2100 <- proj_binary_mean(get_predictions(loadRData(paste0(sps,"/proj_2100/",sps,".2100.projection.out"))), 
                                load_maxent(paste0("data/maxent/",sps,"_2100_45_avg_binary.tif")), 
                                res_table, "2100")
  writeRaster(df_to_raster(proj_2100), paste0("data/biomod/",sps,"_binary_2100.asc"), overwrite = TRUE)
  rm(proj_2100); gc()
  
  # Give the machine a breather
  Sys.sleep(30)
}
  

# 3: Calculate ensemble models --------------------------------------------

# Run one
# sps_binary(sps_names[1])

# Run all
  # NB: Uses to much RAM when running more than a few at a time
registerDoParallel(cores = 4)
plyr::l_ply(sps_names, sps_binary, .parallel = TRUE)

