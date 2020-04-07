# 2_results.R
# The purpose of this script is to take the outputs of 1_biomod.R
# and process them into a format that may be compared to Jesi's Maxent results
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
  
  # Visualise quality of different models
  # models_scores_graph(biomod_model)
  
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

# Convenience function to process a raster into a ggplot friendly dataframe
raster_to_df <- function(raster_file, model_name){
  res <- as.data.frame(raster(raster_file), xy = TRUE) %>%
    na.omit() %>% 
    mutate(x = plyr::round_any(x, 0.25), 
           y = plyr::round_any(y, 0.25)) %>% 
    `colnames<-`(c("x", "y", "val")) %>% 
    group_by(x, y) %>% 
    summarise(val = round(mean(val, na.rm = TRUE))) %>% 
    ungroup() %>% 
    `colnames<-`(c("x", "y", model_name))
  return(res)
}

# Convenience function for plotting
comparison_plot <- function(df){
  comparison_fig <- ggplot(data = df, aes(x = x, y = y)) +
    borders(fill = "grey20", colour = "black") +
    geom_tile(aes(fill = as.factor(Binary))) +
    labs(x = NULL, y = NULL, fill = "Presence") +
    scale_fill_manual(values = c("grey80", "forestgreen")) +
    facet_wrap(~Model) +
    coord_quickmap(expand = F)
  return(comparison_fig)
}

# Wraper to run visuals for all species
biomod_visuals <- function(sps){
  # Load rasters
  # NB: Not all species have their MaxEnt data uploaded yet so this needs to be confirmed first
  if(file.exists(paste0("data/maxent/",sps,"_avg_binary.tif"))){
    
    # Load raster files as formatted dataframes
    ensemble_df <- raster_to_df(paste0(sps,"/proj_present/proj_present_",sps,"_TSSbin.gri"), "ensemble")
    maxent_df <- raster_to_df(paste0("data/maxent/",sps,"_avg_binary.tif"), "maxent")
    ensemble_2050_df <- raster_to_df(paste0(sps,"/proj_2050/proj_2050_",sps,"_TSSbin.gri"), "ensemble_2050")
    maxent_2050_df <- raster_to_df(paste0("data/maxent/",sps,"_2050_45_avg_binary.tif"), "maxent_2050")
    ensemble_2100_df <- raster_to_df(paste0(sps,"/proj_2100/proj_2100_",sps,"_TSSbin.gri"), "ensemble_2100")
    maxent_2100_df <- raster_to_df(paste0("data/maxent/",sps,"_2100_45_avg_binary.tif"), "maxent_2100")
    
    # Combine files for plotting
    both_df <- left_join(ensemble_df, maxent_df, by = c("x", "y")) %>% 
      mutate(both = ifelse(ensemble == 1 & maxent == 1, 1, 0)) %>% 
      pivot_longer(cols = ensemble:both, names_to = "Model", values_to = "Binary")
    both_2050_df <- left_join(ensemble_2050_df, maxent_2050_df, by = c("x", "y")) %>% 
      mutate(both_2050 = ifelse(ensemble_2050 == 1 & maxent_2050 == 1, 1, 0)) %>% 
      pivot_longer(cols = ensemble_2050:both_2050, names_to = "Model", values_to = "Binary")
    both_2100_df <- left_join(ensemble_2100_df, maxent_2100_df, by = c("x", "y")) %>% 
      mutate(both_2100 = ifelse(ensemble_2100 == 1 & maxent_2100 == 1, 1, 0)) %>% 
      pivot_longer(cols = ensemble_2100:both_2100, names_to = "Model", values_to = "Binary")
    
    # Create figures
    comparison_present <- comparison_plot(both_df)
    ggsave(plot = comparison_present, filename = paste0("graph/comparison/",sps,"_comparison.png"), width = 20, height = 4)
    comparison_2050 <- comparison_plot(both_2050_df)
    ggsave(plot = comparison_2050, filename = paste0("graph/comparison/",sps,"_comparison_2050.png"), width = 20, height = 4)
    comparison_2100 <- comparison_plot(both_2100_df)
    ggsave(plot = comparison_2100, filename = paste0("graph/comparison/",sps,"_comparison_2100.png"), width = 20, height = 4)
  }
}

# Run them all
plyr::l_ply(sps_names, biomod_visuals, .parallel = T)

