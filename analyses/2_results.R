# 2_results.R
# The purpose of this script is to take the outputs of 1_biomod.R
# and process them into a format that may be compared to 
# Jesi's Maxent results
# 1: Setup the environment
# 2: Create table of variables used
# 3: Create table of model results
# 4: Create visuals


# 1: Setup ----------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))
library(tidyverse)
library(biomod2)
library(sp)
library(doParallel); registerDoParallel(cores = 50)

# Function for re-loading .RData files as necessary
#loads an RData file, and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# The species occurrence data
sps_names <- str_remove(dir("data/occurrence", full.names = F), pattern = "_near.csv")

# Choose species
# sps <- sps_names[5]


# 2: Create table of variables used ---------------------------------------

# The wrapper functions
biomod_var_table <- function(sps){
  
  # The full model results
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Get list of variables used
  biomod_var <- as.data.frame(biomod_model@expl.var.names) %>% 
    `colnames<-`(c("var")) %>%
    mutate(sps = sps,
           var_count = 1:n()) %>% 
    pivot_wider(names_from = var_count, values_from = var, names_prefix = "Var")
  return(biomod_var)
}

# Run it all
all_var_table <- plyr::ldply(sps_names, biomod_var_table, .parallel = T)

# Save
write_csv(all_var_table, "metadata/all_var_table.csv")


# 3: Create table of model results ----------------------------------------

# The wrapper functions
biomod_res_table <- function(sps){
  
  # The full model results
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Get the TSS scores, cutoffs for binary presence/absence, and specificity/sensitivity
  biomod_cutoff <- plyr::adply(get_evaluations(biomod_model), c(1,3,4,5)) %>% 
    dplyr::rename(test = X1, model = X2, run = X3, PA = X4, score = Testing.data) %>% 
    mutate(sps = sps) %>% 
    dplyr::select(sps, everything()) %>% 
    group_by(sps, test, model) %>% 
    summarise_if(is.numeric, mean) %>% 
    ungroup() %>% 
    mutate_if(is.numeric, round, 2) %>% 
    mutate(Cutoff = round(Cutoff))
  return(biomod_cutoff)
}

# Run it all
all_res_table <- plyr::ldply(sps_names, biomod_res_table, .parallel = T)

# Save
write_csv(all_res_table, "metadata/all_res_table.csv")


# 4: Visuals --------------------------------------------------------------

# Visualise quality of different models
# models_scores_graph(biomod_model)

# Load rasters
# NB: Not all species have their MaxEnt data uploaded yet so this needs to be confirmed first
if(file.exists(paste0("data/maxent/",sps_name,"_binary.tif"))){
  ensemble_raster <- raster(paste0(sps_name,"/proj_present/proj_present_",sps_name,"_TSSbin.gri"))
  maxent_raster <- raster(paste0("data/maxent/",sps_name,"_binary.tif"))
  ensemble_2050_raster <- raster(paste0(sps_name,"/proj_2050/proj_2050_",sps_name,"_TSSbin.gri"))
  maxent_2050_raster <- raster(paste0("data/maxent/",sps_name,"_2050_45_avg_binary.tif"))
  ensemble_2100_raster <- raster(paste0(sps_name,"/proj_2100/proj_2100_",sps_name,"_TSSbin.gri"))
  maxent_2100_raster <- raster(paste0("data/maxent/",sps_name,"_2100_45_avg_binary.tif"))
}

# Convert raster to data frame
ensemble_df <- as.data.frame(ensemble_raster, xy = TRUE) %>% 
  mutate(x = plyr::round_any(x, 0.25), 
         y = plyr::round_any(y, 0.25)) %>% 
  group_by(x, y) %>% 
  summarise(ensemble = round(mean(layer.1, na.rm = TRUE))) %>% 
  ungroup() %>% 
  na.omit()
maxent_df <- as.data.frame(maxent_raster, xy = TRUE) %>% 
  mutate(x = plyr::round_any(x, 0.25), 
         y = plyr::round_any(y, 0.25)) %>% 
  group_by(x, y) %>% 
  summarise(maxent = round(mean(Aebu_binary, na.rm = TRUE))) %>% 
  ungroup() %>% 
  na.omit()
both_df <- left_join(ensemble_df, maxent_df, by = c("x", "y")) %>% 
  mutate(both = ifelse(ensemble == 1 & maxent == 1, 1, 0)) %>% 
  pivot_longer(cols = ensemble:both, names_to = "Model", values_to = "Binary")

comparison_plot <- ggplot(data = both_df, aes(x = x, y = y)) +
  borders(fill = "grey20", colour = "black") +
  geom_tile(aes(fill = as.factor(Binary))) +
  labs(x = NULL, y = NULL, fill = "Presence") +
  scale_fill_manual(values = c("grey80", "forestgreen")) +
  facet_wrap(~Model) +
  coord_quickmap(expand = F)
ggsave(plot = comparison_plot, filename = "graph/Aebu_comparison.png", width = 20, height = 4)


