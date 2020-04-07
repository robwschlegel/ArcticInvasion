# 3_ensemble.R
# The purose of this script is to take all of the results from Bbiomod2
# and MaxEnt and merge them together into one final binary map for each species.
# 1: Setup the environment
# 2: Load biomod2 and MaxEnt results
# 3: Find cutoff values
# 4: Average all results
# 5: Convert the binary results to raster format

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
# sps <- sps_names[5]


# convenience function for loading MaxEnt data
load_maxent <- function(maxent_file){
  res <- as.data.frame(raster(maxent_file), xy = TRUE) %>%
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(model = "MaxEnt") %>% 
    dplyr::select(model, x, y, z)
}

# Function that converts a raster layer to a binary dataframe
raster_to_df_binary <- function(layer_num, proj_x){
  
  # Convert to Dataframe
  proj_sub_df <- as.data.frame(proj_x@layers[[layer_num]], xy = T)
  
  # Find TSS cutoff
  sub_cutoff <- strsplit(colnames(proj_sub_df[3]), "_")[[1]]
  res_table_sub <- res_table %>% 
    filter(run == sub_cutoff[3], model == sub_cutoff[4])
  
  # Convert to binary
  proj_sub_df_binary <- proj_sub_df %>% 
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(z = if_else(z >= res_table_sub$Cutoff[1], 1, 0),
           run = sub_cutoff[3],
           model = sub_cutoff[4]) %>% 
    dplyr::select(model, run, x, y, z)
}

# Function that creates the final mean binary value for a porjection
proj_binary_mean <- function(proj_biomod, proj_maxent){
  
  # Extract all of the models as a long data.frame
  df_biomod <- plyr::ldply(1:length(proj_biomod@layers), raster_to_df_binary,
                           .parallel = F, proj_x = proj_biomod)
  
  # Create rounded mean binary values per model
  # df_biomod_mean <- plyr::ddply(df_present, c("model", "x", "y"), mean, .parallel = T)
  # system.time(
  df_biomod_mean <- lazy_dt(df_biomod) %>% 
    dplyr::select(-run) %>% 
    group_by(model, x, y) %>% 
    summarise(z = round(mean(z))) %>% 
    ungroup() %>% 
    data.frame()
  # ) # 340 seconds
  
  # Add MaxEnt binary and calculate final binary result
  df_res <- rbind(df_biomod_mean, proj_maxent)
  # system.time(
  df_res <- lazy_dt(df_res) %>% 
    group_by(x, y) %>% 
    summarise(z = round(mean(z))) %>% 
    ungroup() %>% 
    data.frame()
  # ) # 63 seconds
  return(df_res)
}

# Function that crawls through species files, finds cutoffs, and creates binary output per model
sps_binary <- function(sps){
  
}


# 2: Load -----------------------------------------------------------------

# The MaxEnt binary results
maxent_present <- load_maxent(paste0("data/maxent/",sps,"_avg_binary.tif"))
maxent_2050 <- load_maxent(paste0("data/maxent/",sps,"_2050_45_avg_binary.tif"))
maxent_2100 <- load_maxent(paste0("data/maxent/",sps,"_2100_45_avg_binary.tif"))

# The biomod non-binary results
biomod_present <- get_predictions(loadRData(paste0(sps,"/proj_present/",sps,".present.projection.out")))
biomod_2050 <- get_predictions(loadRData(paste0(sps,"/proj_2050/",sps,".2050.projection.out")))
biomod_2100 <- get_predictions(loadRData(paste0(sps,"/proj_2100/",sps,".2100.projection.out")))

# The biomod model results
biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))


# 3: Cutoffs --------------------------------------------------------------

# The cutoff table
res_table <- plyr::adply(get_evaluations(biomod_model), c(1,3,4,5)) %>% 
  dplyr::rename(test = X1, model = X2, run = X3, PA = X4, score = Testing.data) %>% 
  mutate(sps = sps) %>% 
  dplyr::select(sps, everything()) %>% 
  filter(test == "TSS")


# 4: Mean binary results --------------------------------------------------

proj_present <- proj_binary_mean(biomod_present, maxent_present)
proj_2050 <- proj_binary_mean(biomod_2050, maxent_2050)
proj_2100 <- proj_binary_mean(biomod_2100, maxent_2100)


# 5: Convert to raster ----------------------------------------------------

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


raster_present <- df_to_raster(proj_present)
writeRaster(raster_present, paste0("data/biomod/",sps,"_binary.asc"))

# test
test_raster <- raster(paste0("data/biomod/",sps,"_binary.asc"))
plot(test_raster)
