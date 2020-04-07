# 3_ensemble.R
# The purose of this script is to take all of the results from Bbiomod2
# and MaxEnt and merge them together into one final binary map for each species.
# 1: Setup the environment
# 2: Load biomod2 and MaxEnt results
# 3: Find TSS cutoff points
# 4: Average all results

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


# 2: Load -----------------------------------------------------------------

# convenience function for loading MaxEnt data
load_MaxEnt <- function(MaxEnt_file){
  res <- as.data.frame(raster(MaxEnt_file), xy = TRUE) %>%
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(model = "MaxEnt") %>% 
    dplyr::select(model, x, y, z)
}

# Function that crawls through species files, finds cutoffs, and creates binary output per model

# The MaxEnt results
maxent_present <- load_MaxEnt(paste0("data/maxent/",sps,"_avg_binary.tif"))
maxent_2050 <- load_MaxEnt(paste0("data/maxent/",sps,"_2050_45_avg_binary.tif"))
maxent_2100 <- load_MaxEnt(paste0("data/maxent/",sps,"_2100_45_avg_binary.tif"))

# The biomod model results
biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))

# The cutoff table
res_table <- plyr::adply(get_evaluations(biomod_model), c(1,3,4,5)) %>% 
  dplyr::rename(test = X1, model = X2, run = X3, PA = X4, score = Testing.data) %>% 
  mutate(sps = sps) %>% 
  dplyr::select(sps, everything()) %>% 
  filter(test == "TSS")

# raster_present <- raster(paste0(sps,"/proj_present/proj_present_",sps,".gri"))
# plot(res_present)

# The raster results
res_present <- loadRData(paste0(sps,"/proj_present/",sps,".present.projection.out"))
proj_present <- get_predictions(res_present)

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

# Extract all of the models as a long data.frame
df_present <- plyr::ldply(1:length(proj_present@layers), raster_to_df_binary, 
                          .parallel = F, proj_x = proj_present)

# Create rounded mean binary values per model
df_present_mean <- lazy_dt(df_present) %>% 
  group_by(model, x, y) %>% 
  summarise(z = round(mean(z))) %>% 
  ungroup() %>% 
  data.frame()

# Add MaxEnt binary and calculate final binary result
df_present_all

