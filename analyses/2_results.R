# 2_results.R
# The purpose of this script is to take the outputs of 1_biomod.R
# and process them into a format that may be compared to Jesi's Maxent results
# 1: Setup the environment
# 2: Create table of variables used
# 3: Create table of model results
# 4: Create basic comparisons
# 5: Create multi-model comparisons
# 6: Create non-binary comparisons


# 1: Setup ----------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))
library(tidyverse)
library(ggpubr)
library(biomod2)
library(sp)
library(dtplyr)
library(doParallel); registerDoParallel(cores = 50)

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# The base map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>%
  dplyr::rename(lon = long) %>%
  mutate(group = ifelse(lon > 180, group+9999, group),
         lon = ifelse(lon > 180, lon-360, lon))

# The species occurrence data
sps_files <- dir("data/occurrence", full.names = T)

# The species abbreviations
sps_names <- str_remove(dir("data/occurrence", full.names = F), pattern = "_near.csv")

# The species depth masks and full names
sps_depths <- read_csv("metadata/sps_depth.csv")

# Choose species
# sps <- sps_names[4]


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
# all_var_table <- plyr::ldply(sps_names, biomod_var_table, .parallel = T)

# Save
# write_csv(all_var_table, "metadata/all_var_table.csv")


# 3: Create table of model results ----------------------------------------

# The wrapper functions
biomod_res_table <- function(sps){
  
  # The full model results
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  model_scores <- models_scores_graph(biomod_model, plot = F)
  # ggsave(filename = paste0("graph/model_scores/",sps,"_scores.png"), plot = model_scores)
  
  # Get the TSS scores, cutoffs for binary presence/absence, and specificity/sensitivity
  biomod_cutoff <- plyr::adply(get_evaluations(biomod_model), c(1,3,4,5)) %>% 
    dplyr::rename(test = X1, model = X2, run = X3, PA = X4, score = Testing.data) %>% 
    mutate(sps = sps) %>% 
    dplyr::select(sps, everything())
  
  return(biomod_cutoff)
}

# Run it all
# all_res_table <- plyr::ldply(sps_names, biomod_res_table, .parallel = T)
# write_csv(all_res_table, "metadata/all_res_table.csv")

# Create a table that shows the mean results
# mean_res_table <- all_res_table %>%
#   group_by(sps, test, model) %>%
#   summarise_if(is.numeric, mean) %>%
#   ungroup() %>%
#   mutate_if(is.numeric, round, 2) %>%
#   mutate(Cutoff = round(Cutoff))
# write_csv(mean_res_table, "metadata/mean_res_table.csv")


# 3: Create basic comparisons ---------------------------------------------

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
    geom_tile(aes(fill = as.factor(z))) +
    labs(x = NULL, y = NULL, fill = "Presence") +
    scale_fill_manual(values = c("grey80", "forestgreen")) +
    facet_wrap(~model) +
    coord_quickmap(expand = F) +
    theme(legend.position = "bottom")
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
      pivot_longer(cols = ensemble:both, names_to = "model", values_to = "z")
    both_2050_df <- left_join(ensemble_2050_df, maxent_2050_df, by = c("x", "y")) %>% 
      mutate(both_2050 = ifelse(ensemble_2050 == 1 & maxent_2050 == 1, 1, 0)) %>% 
      pivot_longer(cols = ensemble_2050:both_2050, names_to = "model", values_to = "z")
    both_2100_df <- left_join(ensemble_2100_df, maxent_2100_df, by = c("x", "y")) %>% 
      mutate(both_2100 = ifelse(ensemble_2100 == 1 & maxent_2100 == 1, 1, 0)) %>% 
      pivot_longer(cols = ensemble_2100:both_2100, names_to = "model", values_to = "z")
    
    # Create figures
    comparison_present <- comparison_plot(both_df)
    ggsave(plot = comparison_present, filename = paste0("graph/comparison/",sps,"_comparison.png"), width = 20, height = 4)
    comparison_2050 <- comparison_plot(both_2050_df)
    ggsave(plot = comparison_2050, filename = paste0("graph/comparison/",sps,"_comparison_2050.png"), width = 20, height = 4)
    comparison_2100 <- comparison_plot(both_2100_df)
    ggsave(plot = comparison_2100, filename = paste0("graph/comparison/",sps,"_comparison_2100.png"), width = 20, height = 4)
  }
}

# Run one
# biomod_visuals(sps_names[1])

# Run them all
# plyr::l_ply(sps_names, biomod_visuals, .parallel = T)


# 5: Create multi-model comparisons ---------------------------------------
# NB: This section relies on data created in 'analyses/3_ensemble.R'

# Function that loads an .Rds file and rounds it to the nearest 0.25 degree resolution
readRDS_0.25 <- function(file_name, projection_name){
  df <- readRDS(file_name) %>% 
    na.omit() %>% 
    mutate(x = plyr::round_any(x, 0.25), 
           y = plyr::round_any(y, 0.25),
           projection = projection_name) %>% 
    group_by(model, projection, x, y) %>% 
    summarise(z = round(mean(z, na.rm = T))) %>% 
    ungroup()
}

# Convenience function to process a raster into a long dataframe
raster_to_long <- function(raster_file, model_name, projection_name){
  res <- as.data.frame(raster(raster_file), xy = TRUE) %>%
    na.omit() %>% 
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(x = plyr::round_any(x, 0.25), 
           y = plyr::round_any(y, 0.25),
           model = model_name,
           projection = projection_name) %>% 
    group_by(model, projection, x, y) %>% 
    summarise(z = round(mean(z, na.rm = TRUE))) %>% 
    ungroup()
  return(res)
}

# Load depth for result screening
depth_long <- raster_to_long("data/present/depthclip.asc", "depth", "present") %>% 
  dplyr::rename(depth = z) %>% 
  dplyr::select(x, y, depth)

# Convenience function for multi-plotting
single_plot <- function(df){
  comp_multi_fig <- ggplot(data = df, aes(x = x, y = y)) +
    geom_polygon(data = map_base, fill = "grey70",
                 aes(x = lon, y = lat, group = group)) +
    geom_tile(aes(fill = as.factor(z))) +
    # geom_tile(fill = "red4") +
    labs(x = NULL, y = NULL, 
         fill = "Predicted habitat suitability") +
    scale_fill_manual(values = c("red4")) +
    facet_wrap(~model) +
    coord_quickmap(expand = F) +
    # theme_void() + # This interferes with ggarrange()
    theme(legend.position = "bottom", 
          title = element_text(size = 16),
          strip.text = element_text(size = 12),
          legend.text = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank(), 
          panel.grid = element_blank())
  return(comp_multi_fig)
}

comp_multi_plot <- function(df, sps){
  
  # Find species depth limit
  sps_depth <- filter(sps_depths, Sps == sps)
  
  # Load species presence data
  sps_data <- read_csv(sps_files[which(sps == sps_names)]) %>% 
    `colnames<-`(c("Sps", "lon", "lat")) %>% 
    mutate(lon =  plyr::round_any(lon, 0.25),
           lat =  plyr::round_any(lat, 0.25))
  
  # Prep data for plotting
  df_sub <- df %>% 
    left_join(depth_long, by = c("x", "y")) %>% 
    mutate(z = ifelse(depth > sps_depth$Depth, 0, z),
           model = factor(model, levels = unique(model))) %>% 
    filter(z != 0)
  
  # Create each panel
  plot_ANN <- single_plot(filter(df_sub, model == "ANN")) +
    geom_point(data = sps_data, aes(x = lon, y = lat), 
               fill = "white", shape = 21, size = 0.5, stroke = 0.2) #+
    # labs(title = sps_depth$Species)
  plot_GAM <- single_plot(filter(df_sub, model == "GAM"))
  plot_GLM <- single_plot(filter(df_sub, model == "GLM"))
  plot_RF <- single_plot(filter(df_sub, model == "RF"))
  plot_MaxEnt <- single_plot(filter(df_sub, model == "MaxEnt"))
  plot_Ensemble <- single_plot(filter(df_sub, model == "Ensemble")) +
    theme(panel.background = element_rect(colour = "black", size = 2))
  
  # Steek'em
  plot_multi <- ggarrange(plot_ANN, plot_GAM, plot_GLM, plot_RF, plot_MaxEnt, plot_Ensemble, 
                          ncol = 3, nrow = 2, common.legend = T, legend = "bottom", align = "hv")
  plot_multi_title <- annotate_figure(p = plot_multi, fig.lab = df_sub$projection[1], fig.lab.face = "bold", fig.lab.size = 14,
                                      top = text_grob(sps_depth$Species, color = "black", face = "italic", size = 16))
  
  # Save
  ggsave(plot = plot_multi_title, width = 21, height = 8,
         filename = paste0("graph/comparison_multi/",sps,"_",df$projection[1],".png"))
}

# Wraper to run visuals for all species
biomod_multi_visuals <- function(sps){
  
  # Create present figure
  df_present <- rbind(readRDS_0.25(paste0("data/biomod/",sps,"_df_present.Rds"), "present"),
                      raster_to_long(paste0("data/maxent/",sps,"_avg_binary.tif"), "MaxEnt", "present"),
                      raster_to_long(paste0(sps,"/proj_present/proj_present_",sps,"_TSSbin.gri"), "Ensemble", "present"))
  comp_multi_present <- comp_multi_plot(df_present, sps)
  rm(df_present, comp_multi_present); gc()
  
  # Create 2050 figure
  df_2050 <- rbind(readRDS_0.25(paste0("data/biomod/",sps,"_df_2050.Rds"), "2050"),
                   raster_to_long(paste0("data/maxent/",sps,"_2050_45_avg_binary.tif"), "MaxEnt", "2050"),
                   raster_to_long(paste0(sps,"/proj_2050/proj_2050_",sps,"_TSSbin.gri"), "Ensemble", "2050"))
  comp_multi_2050 <- comp_multi_plot(df_2050, sps)
  rm(df_2050, comp_multi_2050); gc()
  
  # Create 2100 figure
  df_2100 <- rbind(readRDS_0.25(paste0("data/biomod/",sps,"_df_2100.Rds"), "2100"),
                        raster_to_long(paste0("data/maxent/",sps,"_2100_45_avg_binary.tif"), "MaxEnt", "2100"),
                        raster_to_long(paste0(sps,"/proj_2100/proj_2100_",sps,"_TSSbin.gri"), "Ensemble", "2100"))
  comp_multi_2100 <- comp_multi_plot(df_2100, sps)
  rm(df_2100, comp_multi_2100); gc()
}

# Run one
# system.time(biomod_multi_visuals(sps_names[1])) # 223 seconds

# Run them all
plyr::l_ply(sps_names, biomod_multi_visuals, .parallel = T)


# 6: Create non-binary comparisons ----------------------------------------

# Function that converts a raster layer to a dataframe
raster_to_df_0.25 <- function(layer_num, proj_x, proj_name){
  
  # Base data
  proj_base_df <- as.data.frame(proj_x@layers[[layer_num]], xy = T)
  
  # Info
  df_info <- strsplit(colnames(proj_base_df[3]), "_")[[1]]
  
  # Prep for datatable
  proj_prep_df <- proj_base_df %>% 
    na.omit() %>% 
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(x = plyr::round_any(x, 0.25), 
           y = plyr::round_any(y, 0.25),
           run = df_info[3],
           model = paste0(df_info[4],"_",proj_name))
  
  # Result
  proj_df <- lazy_dt(proj_prep_df) %>% 
    group_by(model, run, x, y) %>% 
    summarise(z = round(mean(z, na.rm = TRUE))) %>% 
    ungroup() %>% 
    data.frame()
  return(proj_df)
}

# Function that combines every non-binary model proj for each species
proj_combine_0.25 <- function(proj_biomod, proj_maxent, proj_name, sps){
  
  # Extract all of the models as a long data.frame
  # system.time(
  df_biomod <- plyr::ldply(1:length(proj_biomod@layers), raster_to_df_0.25, 
                           .parallel = F, proj_biomod, proj_name)
  # ) # 186 seconds
  
  # Create rounded mean non-binary values per model
  # system.time(
  df_biomod_mean <- lazy_dt(df_biomod) %>% 
    dplyr::select(-run) %>% 
    group_by(model, x, y) %>% 
    summarise(z = round(mean(z))) %>% 
    ungroup() %>% 
    data.frame()
  # ) # 20 seconds
  saveRDS(df_biomod_mean, paste0("data/biomod/",sps,"_df_non_",proj_name,".Rds"))
  # df_biomod_mean <- readRDS(paste0("data/biomod/",sps,"_df_non_",proj_name,".Rds"))
  rm(df_biomod); gc()
  
  # Prep maxent df for row binding
  df_maxent <- proj_maxent %>% 
    `colnames<-`(c("x", "y", "z")) %>% 
    mutate(model = paste0("MaxEnt_",proj_name),
           z = z*1000) %>% 
    dplyr::select(model, x, y, z)
  
  # Add MaxEnt binary and calculate final binary result
  df_res <- rbind(df_biomod_mean, df_maxent)
  rm(df_biomod_mean, df_maxent); gc()
  return(df_res)
}

# Convenience function for plotting
non_binary_plot <- function(df){
  non_binary_fig <- ggplot(data = df, aes(x = x, y = y)) +
    borders(fill = "grey20", colour = "black") +
    geom_tile(aes(fill = z)) +
    labs(x = NULL, y = NULL, fill = "Suitability") +
    scale_fill_gradient(low = "grey80", high = "forestgreen") +
    facet_wrap(~model) +
    coord_quickmap(expand = F) +
    theme(legend.position = "bottom")
  return(non_binary_fig)
}

# Function for creating all of the non-binary visuals
biomod_non_binary_visuals <- function(sps){
  
  # Load raster files and create figures
  ## Present
  both_present_df <- proj_combine_0.25(get_predictions(loadRData(paste0(sps,"/proj_present/",sps,".present.projection.out"))),
                                       raster_to_df(paste0("data/maxent/",sps,"_avg_binary.tif"), "maxent"), "present", sps)
  non_binary_present <- non_binary_plot(both_present_df)
  ggsave(plot = non_binary_present, filename = paste0("graph/comparison_non_binary/",sps,"_present.png"), width = 20, height = 8)
  rm(both_present_df, non_binary_present); gc()
  ## 2050
  both_2050_df <- proj_combine_0.25(get_predictions(loadRData(paste0(sps,"/proj_2050/",sps,".2050.projection.out"))),
                                    raster_to_df(paste0("data/maxent/",sps,"_2050_45_avg_binary.tif"), "maxent"), "2050", sps)
  non_binary_2050 <- non_binary_plot(both_2050_df)
  ggsave(plot = non_binary_2050, filename = paste0("graph/comparison_non_binary/",sps,"_2050.png"), width = 20, height = 8)
  rm(both_2050_df, non_binary_2050); gc()
  ## 2100
  both_2100_df <- proj_combine_0.25(get_predictions(loadRData(paste0(sps,"/proj_2100/",sps,".2100.projection.out"))),
                                    raster_to_df(paste0("data/maxent/",sps,"_2100_45_avg_binary.tif"), "maxent"), "2100", sps)
  non_binary_2100 <- non_binary_plot(both_2100_df)
  ggsave(plot = non_binary_2100, filename = paste0("graph/comparison_non_binary/",sps,"_2100.png"), width = 20, height = 8)
  rm(both_2100_df, non_binary_2100); gc()
}

# Run it all
# plyr::l_ply(sps_names, biomod_non_binary_visuals, .parallel = T)

