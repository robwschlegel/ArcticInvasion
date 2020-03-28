# 1_biomod.R
# The purpose of this script is to run the BIOMOD model ensemble
# on the many potential alien invasive species in the Arctic.
# The data have already been prepared by Jesi for use in her Maxent model.
# The order of operations is:
# 1: Setup the environment for biomod_pipeline()
# 2: Load data
# 3: Prep data for modelling
# 4: Run model ensemble
# 5: Present projections
# 6: Future projections
# 7: Run the full pipeline


# 1: Setup ----------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))
library(tidyverse)
library(biomod2)
library(sp)
library(FNN)
library(doParallel); registerDoParallel(cores = 50)

# NB: The data are housed on dropbox on Jesi Goldsmit's professional account
# They have been downloaded locally to the dropbox folder on my machine
# I then created a symbolic link from there to this project folder

# The species occurrence data
sps_files <- dir("data/occurrence", full.names = T)

# The environmental file pathways
var_files <- dir("data/present", full.names = T)
var_2050_files <- dir("data/2050", full.names = T)
var_2100_files <- dir("data/2100", full.names = T)

# Add the nutrient files to the future files
var_2050_files <- c(var_2050_files,
                    var_files[which(!sapply(str_split(var_files, "/"), "[[", 3) %in% sapply(str_split(var_2050_files, "/"), "[[", 3))])
var_2100_files <- c(var_2100_files,
                    var_files[which(!sapply(str_split(var_files, "/"), "[[", 3) %in% sapply(str_split(var_2100_files, "/"), "[[", 3))])

# Global coords from Jesi's data
global_coords <- as.data.frame(sp::read.asciigrid(var_files[1]), xy = T)
global_coords$env_index <- 1:nrow(global_coords)

# The best variables per species
top_var <- read_csv("metadata/top_var.csv") %>% 
  dplyr::select(Code:var6) %>% 
  pivot_longer(cols = var1:var6) %>% 
  dplyr::select(-name) %>% 
  na.omit()

# The full biomod pipeline in one function
biomod_pipeline <- function(sps_choice){
  
  # 2: Load data ------------------------------------------------------------
  
  # Load the species
  sps <- read_csv(sps_choice) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("s1", "s2")]),
                                            as.matrix(.[,c("lon", "lat")]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(sps, s1, s2) %>%
    dplyr::rename(lon = s1, lat = s2)
  
  # Filter out the top variables
  top_var_sub <- top_var %>% 
    filter(str_remove(Code, pattern = "_near") == sps$sps[1]) %>% 
    mutate(value = paste0(value,".asc"))
  
  # Load the top variables for the species
  expl <- raster::stack(var_files[which(sapply(str_split(var_files, "/"), "[[", 3) %in% top_var_sub$value)])
  expl_2050 <- raster::stack(var_2050_files[which(sapply(str_split(var_2050_files, "/"), "[[", 3) %in% top_var_sub$value)])
  expl_2100 <- raster::stack(var_2100_files[which(sapply(str_split(var_2100_files, "/"), "[[", 3) %in% top_var_sub$value)])
  
  # 3: Prep data ------------------------------------------------------------
  
  # Prep data for modelling
  biomod_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(sps)),
    resp.xy = as.matrix(sps[,2:3]),
    resp.name = sps$sps[1],
    # eval.resp.var = rep(1, nrow(sps_test)), # Doesn't work with presence only...
    # eval.resp.xy = as.matrix(sps_test[,2:3]),
    expl.var = expl, 
    PA.nb.rep = 5,
    PA.strategy = "sre",
    PA.sre.quant = 0.1)
  
  # check data format
  # biomod_data
  
  # Check plot of data
  # plot(biomod_data)
  
  # Model options
  biomod_option <- BIOMOD_ModelingOptions()
  
  # Cross validation
  # biomod_cv <- BIOMOD_cv(biomod_data)
  
  # 4: Model ----------------------------------------------------------------
  
  # Run the model
  biomod_model <- BIOMOD_Modeling(
    biomod_data,
    models = c('GLM', 'GAM', 'ANN', 'SRE', 'RF'),
    models.options = biomod_option,
    NbRunEval = 3,
    DataSplit = 70,
    VarImport = 1,
    models.eval.meth = c('TSS','ROC', 'KAPPA', 'ACCURACY', 'BIAS'),
    rescal.all.models = TRUE,
    do.full.models = T,
    modeling.id = sps$sps[1])
  
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
  
  # 5. Present projections --------------------------------------------------
  
  # biomod_model@models.computed
  
  # For ensemble forecast from all models
  # BIOMOD_EnsembleForecasting()
  
  # Create projections
  biomod_projection <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl,
    proj.name = 'present',
    # selected.models = biomod_model@models.computed,#[13:16],
    binary.meth = 'TSS',
    Bin.trans = TRUE,
    slot = biomod_model@models.computed,
    compress = FALSE,
    build.clamping.mask = FALSE)
  
  # Plot projections
  # plot(biomod_projection)
  
  # Get predictions
  # NB: There are many biomod2::get_ functions for looking at results more closely
  # present_predictions <- get_predictions(biomod_projection)
  # present_predictions
  
  # Look at particular aspects of predictions
  # biomod2::free()
  
  # 6: Future projections ---------------------------------------------------
  
  # Run 2050 projections
  biomod_projection_2050 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl_2050,
    proj.name = '2050',
    # selected.models = biomod_model@models.computed,#[13:16],
    binary.meth = 'TSS',
    Bin.trans = TRUE,
    slot = biomod_model@models.computed,
    compress = FALSE,
    build.clamping.mask = FALSE)
  
  # Plot projections
  # plot(biomod_projection_2050)
  
  # Get projections
  # projection_2050 <- get_predictions(biomod_projection_2050)
  # projection_2050
  
  # Somehow use this to compare the projections over time
  # BIOMOD_RangeSize()
  
  # Run 2100 projections
  biomod_projection_2100 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl_2100,
    proj.name = '2100',
    # selected.models = biomod_model@models.computed,#[13:16],
    binary.meth = 'TSS',
    Bin.trans = TRUE,
    slot = biomod_model@models.computed,
    compress = FALSE,
    build.clamping.mask = FALSE)
}

# 7: Run the full pipeline ------------------------------------------------

plyr::l_ply(sps_files, biomod_pipeline, .parallel = T)

# Testing the loading of files saved automatically to disk
# load("Aebu/models/1585165610/Aebu_PA1_RUN1_ANN")
# get_formal_model(Aebu_PA1_RUN1_ANN) 
# summary(get_formal_model(Aebu_PA1_RUN1_ANN))
# load("Aebu/Aebu.1585165851.models.out")
# test_model <- BIOMOD_LoadModels(Aebu.1585165851.models.out, models = 'RF')

