# 1_biomod.R
# The purpose of this script is to run the BIOMOD model ensemble
# on the many potential alien invasive species in the Arctic.
# The data have already been prepared by Jesi for use in her Maxent model.
# The order of operations is:
# 1: Setup the environment
# 2: Load data
# 3: Prep data for modelling
# 4: Run model ensemble
# 5: Present projections
# 6: Future projections
# 7: Run the pipeline


# 1: Setup ----------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))
library(tidyverse)
library(biomod2)
library(sp)
library(FNN)
library(doParallel); registerDoParallel(cores = 50)

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

# Function for re-loading .RData files as necessary
  #loads an RData file, and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Choose a species
# sps_choice <- sps_files[18]
# NB: Focus on zooplankton group first
# c(2, 4, 18)

# The full pipeline wrapped into a function
biomod_pipeline <- function(sps_choice){
  
  print(paste0("Began run on ",sps_choice))
  
  
  # 2: Load data ------------------------------------------------------------
  
  # Load the species
  sps <- read_csv(sps_choice) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("s1", "s2")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sps, s1, s2) %>%
    dplyr::rename(lon = s1, lat = s2)
  
  # The species abbreviation
  sps_name <- sps$Sps[1]
  
  # Determine number of pseudo-absences to use
  if(nrow(sps) <= 1000){
    PA_absence_count <- 1000
  } else{
    PA_absence_count <- 10000
  }
  
  # Filter out the top variables
  top_var_sub <- top_var %>% 
    filter(str_remove(Code, pattern = "_near") == sps_name) %>% 
    mutate(value = paste0(value,".asc"))
  
  # Load the top variables for the species
  expl <- raster::stack(var_files[which(sapply(str_split(var_files, "/"), "[[", 3) %in% top_var_sub$value)])
  expl_2050 <- raster::stack(var_2050_files[which(sapply(str_split(var_2050_files, "/"), "[[", 3) %in% top_var_sub$value)])
  expl_2100 <- raster::stack(var_2100_files[which(sapply(str_split(var_2100_files, "/"), "[[", 3) %in% top_var_sub$value)])
  
  # Set temp folder save locations
  # http://www.r-forge.r-project.org/forum/forum.php?thread_id=30946&forum_id=995&group_id=302
  dir.create (file.path(sps_name), showWarnings = FALSE)
  dir.create (file.path(sps_name,"/Temp"), showWarnings = FALSE)
  rasterOptions(tmpdir = paste0(sps_name,"/Temp"))
  
  
  # 3: Prep data ------------------------------------------------------------
  
  # Prep data for modelling
  biomod_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(sps)),
    resp.xy = as.matrix(sps[,2:3]),
    resp.name = sps_name,
    expl.var = expl, 
    PA.nb.rep = 5,
    PA.nb.absences = PA_absence_count)
  # biomod_data <- readRDS(paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # Model options
  biomod_option <- BIOMOD_ModelingOptions()
  
  
  # 4: Model ----------------------------------------------------------------
  
  # Run the model
  biomod_model <- BIOMOD_Modeling(
    biomod_data,
    models = c('GLM', 'ANN', 'SRE', 'RF'),
    models.options = biomod_option,
    NbRunEval = 3,
    DataSplit = 70,
    VarImport = 0,
    models.eval.meth = c('KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS'),
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = sps_name)
  # biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
  
  # Build the ensemble models
    # It looks like the ensembles aren't going to be needed
    # This is because the projections produce the needed binary output files
  # biomod_ensemble <- BIOMOD_EnsembleModeling(
  #   modeling.output = biomod_model,
  #   eval.metric = 'TSS',
  #   eval.metric.quality.threshold = 0.7,
  #   models.eval.meth = 'TSS'
  # )
  # biomod_ensemble <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"ensemble.models.out"))
  
  # Save the pre-model data for possible later use
  saveRDS(biomod_data, file = paste0(sps_name,"/",sps_name,".base.Rds"))
  
  
  # 5. Present projections --------------------------------------------------
  
  # Create projections
  biomod_projection <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl,
    proj.name = 'present',
    binary.meth = 'TSS',
    compress = "xz",
    build.clamping.mask = FALSE)
  
  # Create ensemble projections
  # biomod_ensemble_projection <- BIOMOD_EnsembleForecasting(
  #   EM.output = biomod_ensemble,
  #   projection.output = biomod_projection)
  
  # Clean out some space
  rm(biomod_projection); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  # dir(tempdir())
  
  
  # 6: Future projections ---------------------------------------------------
  
  # Run 2050 projections
  biomod_projection_2050 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl_2050,
    proj.name = '2050',
    binary.meth = 'TSS',
    compress = 'xz',
    build.clamping.mask = FALSE)
  
  # Create 2050 ensemble projections
  # biomod_ensemble_projection_2050 <- BIOMOD_EnsembleForecasting(
  #   EM.output = biomod_ensemble,
  #   projection.output = biomod_projection_2050)
  
  # Clean out 2050
  rm(biomod_projection_2050); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Run 2100 projections
  biomod_projection_2100 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl_2100,
    proj.name = '2100',
    binary.meth = 'TSS',
    compress = 'xz',
    build.clamping.mask = FALSE)
  
  # Create 2100 ensemble projections
  # biomod_ensemble_projection_2100 <- BIOMOD_EnsembleForecasting(
  #   EM.output = biomod_ensemble,
  #   projection.output = biomod_projection_2100)
  
  # Clean out 2100
  rm(biomod_projection_2100); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Unlink Temp folder
  unlink(paste0(sps_name,"/Temp"), recursive = TRUE)
}


# 7: Run the pipeline -----------------------------------------------------

plyr::l_ply(sps_files[7:8], biomod_pipeline, .parallel = TRUE)

