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
sps_choice <- sps_files[2]
print(paste0("Began run on ",sps_choice))

# NB: Focus on zooplankton group first
# c(2, 4, 18)
  

# 2: Load data ------------------------------------------------------------

# Load the species
sps <- read_csv(sps_choice) %>% 
  mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("s1", "s2")]),
                                          as.matrix(.[,c("lon", "lat")]), k = 1))) %>%
  left_join(global_coords, by = "env_index") %>% 
  dplyr::select(Sps, s1, s2) %>%
  dplyr::rename(lon = s1, lat = s2)

# Filter out the top variables
top_var_sub <- top_var %>% 
  filter(str_remove(Code, pattern = "_near") == sps$Sps[1]) %>% 
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
  resp.name = sps$Sps[1],
  # eval.resp.var = rep(1, nrow(sps_test)), # Doesn't work with presence only...
  # eval.resp.xy = as.matrix(sps_test[,2:3]),
  expl.var = expl, 
  PA.nb.rep = 5,
  PA.strategy = "sre",
  PA.sre.quant = 0.1)
# biomod_data <- readRDS(paste0(sps$Sps[1],"/",sps$Sps[1],".base.Rds"))

# Model options
biomod_option <- BIOMOD_ModelingOptions()


# 4: Model ----------------------------------------------------------------

# Run the model
biomod_model <- BIOMOD_Modeling(
  biomod_data,
  models = c('GLM', 'ANN', 'SRE', 'RF'),#'GAM', ,
  models.options = biomod_option,
  NbRunEval = 3,
  DataSplit = 70,
  VarImport = 0,
  models.eval.meth = c('KAPPA', 'TSS', 'ROC', 'ACCURACY', 'BIAS'),
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = sps$Sps[1])
biomod_model <- loadRData(paste0(sps$Sps[1],"/",sps$Sps[1],".",sps$Sps[1],".models.out"))

# Build the ensemble models
biomod_ensemble <- BIOMOD_EnsembleModeling(
  modeling.output = biomod_model,
  eval.metric = 'TSS',
  eval.metric.quality.threshold = 0.7,
  models.eval.meth = 'TSS'#,
  # prob.ci = TRUE
)
# biomod_ensemble <- loadRData(paste0(sps$Sps[1],"/",sps$Sps[1],".",sps$Sps[1],"ensemble.models.out"))

# Save the pre-model data for possible later use
saveRDS(biomod_data, file = paste0(sps$Sps[1],"/",sps$Sps[1],".base.Rds"))


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
biomod_ensemble_projection <- BIOMOD_EnsembleForecasting(
  EM.output = biomod_ensemble,
  projection.output = biomod_projection)

# Clean out some space
rm(biomod_projection, biomod_ensemble_projection); gc()


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
biomod_ensemble_projection_2050 <- BIOMOD_EnsembleForecasting(
  EM.output = biomod_ensemble,
  projection.output = biomod_projection_2050)

# Clean out 2050
rm(biomod_projection_2050, biomod_ensemble_projection_2050); gc()

# Run 2100 projections
biomod_projection_2100 <- BIOMOD_Projection(
  modeling.output = biomod_model,
  new.env = expl_2100,
  proj.name = '2100',
  binary.meth = 'TSS',
  compress = 'xz',
  build.clamping.mask = FALSE)

# Create 2100 ensemble projections
biomod_ensemble_projection_2100 <- BIOMOD_EnsembleForecasting(
  EM.output = biomod_ensemble,
  projection.output = biomod_projection_2100)

# Clean out 2100
rm(biomod_projection, biomod_ensemble_projection); gc()

