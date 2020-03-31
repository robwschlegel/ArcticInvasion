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
# library(FNN)
library(doParallel); registerDoParallel(cores = 50)


# 2: load biomod ----------------------------------------------------------

# Function for re-loading .RData files as necessary
#loads an RData file, and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

biomod_ensemble_projection <- loadRData("Aebu/proj_present/Aebu.present.ensemble.projection.out")
plot(biomod_ensemble_projection)

test_raster <- raster("Aebu/proj_present/proj_present_Aebu_TSSbin.gri")
plot(test_raster)

writeRaster(test_raster, "Aebu/proj_present/proj_present_Aebu_TSSbin.asc")


# 3: Results --------------------------------------------------------------

# check data format
# biomod_data

# Check plot of data
# plot(biomod_data)

# Have a look at the outputs
# biomod_model

# Relative importance of exploratory variables
# variable_importances <- get_variables_importance(biomod_model)
# variable_importances

# Get all models evaluation
# evaluate()
# biomod_eval <- get_evaluations(biomod_model)
# dimnames(biomod_eval)

# biomod_eval[,,,"RUN1","PA5"]

# Visualise quality of different models
# models_scores_graph(biomod_model)

# biomod_model@models.computed

# Plot projections
# plot(biomod_projection)

# Get predictions
# NB: There are many biomod2::get_ functions for looking at results more closely
# present_predictions <- get_predictions(biomod_projection)
# present_predictions

# Look at particular aspects of predictions
# biomod2::free()

# 4: Visuals --------------------------------------------------------------

# Convert raster to data frame
test_df <- as.data.frame(test_raster, xy = TRUE) %>% 
  mutate(x = plyr::round_any(x, 0.25), 
         y = plyr::round_any(y, 0.25)) %>% 
  group_by(x, y) %>% 
  summarise(mean_layer = round(mean(layer.1, na.rm = TRUE))) %>% 
  ungroup() %>% 
  na.omit()

ggplot(data = test_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = as.factor(mean_layer)))


# 5: Process --------------------------------------------------------------


# 6: Saved ----------------------------------------------------------------


