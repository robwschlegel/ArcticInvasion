# 2_results.R
# The purpose of this script is to take the outputs of 1_biomod.R
# and process them into a format that may be compared to 
# Jesi's Maxent results
# 1: Setup the environment
# 2: Load biomod results
# 3: Look at results
# 4: Create visuals
# 5: Process results into usable outputs
# 6: Save processed outputs


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
sps_files_short <- dir("data/occurrence", full.names = F)

# Choose species
sps_choice <- sps_files_short[5]

# Extract name abreviation
sps_name <- str_remove(sps_choice, pattern = "_near.csv")


# 2: Load biomod ----------------------------------------------------------


# biomod_ensemble_projection <- loadRData("Aebu/proj_present/Aebu.present.ensemble.projection.out")
# plot(biomod_ensemble_projection)

# biomod_data <- readRDS("Bsch/Bsch.base.Rds")
# file_name <- "Bsch/Bsch.Bsch.models.out"
biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
# biomod_projection <- loadRData("Bsch/proj_present/Bsch.present.projection.out")


# 3: Results --------------------------------------------------------------


# Get list of variables used
biomod_var <- biomod_model@expl.var.names

# Get the TSS scores, cutoffs for binary presence/absence, and specificity/sensitivity
biomod_cutoff <- get_evaluations(biomod_model)

# Visualise quality of different models
# models_scores_graph(biomod_model)


# 4: Visuals --------------------------------------------------------------

# Load rasters
ensemble_raster <- raster("Aebu/proj_present/proj_present_Aebu_TSSbin.gri")
maxent_raster <- raster("Aebu/Aebu_binary.tif")

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


# 5: Process --------------------------------------------------------------


# 6: Saved ----------------------------------------------------------------


