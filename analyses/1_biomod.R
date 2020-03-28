# 1_biomod.R
# The purpose of this script is to run the BIOMOD model ensemble
# on the many potential alien invasive species in the Arctic.
# The data have already been prepared by Jesi for use in her Maxent model.
# The order of operations is:
# 1: Load libraries
# 2: Load data
# 3: Prep data for modelling
# 4: Run model ensemble
# 5: Present projections
# 6: Future projections


# 1: Libraries ------------------------------------------------------------

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
var_2050_files <- dir("data/future/2050", full.names = T)
var_2100_files <- dir("data/future/2100", full.names = T)


# 2: Load data ------------------------------------------------------------

# Arctic coords from BioOracle
# load("~/ArcticKelp/data/Arctic_BO.RData")
# ggplot(data = Arctic_BO, aes(x = lon, y = lat)) +
#   geom_raster(aes(fill = BO2_templtmin_bdmax))

# Global coords from Jesi's data
global_coords <- as.data.frame(sp::read.asciigrid(var_files[1]), xy = T)
global_coords$env_index <- 1:nrow(global_coords)

# The best variables per species
## Need to create this from supp table 3

# Load a test species
sps <- read_csv(sps_files[1]) %>% 
  mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("s1", "s2")]),
                                          as.matrix(.[,c("lon", "lat")]), k = 1))) %>%
  left_join(global_coords, by = "env_index") %>% 
  dplyr::select(sps, s1, s2) #%>%
  # dplyr::rename(lon = s1, lat = s2)

# Convert it to spatial points dataframe
sps_spatial <- left_join(global_coords, sps, by = c("s1", "s2")) %>% 
  dplyr::select(s1, s2, sps)
coordinates(sps_spatial) <- ~ s1 + s2
gridded(sps_spatial) <- TRUE

# Visualise
# plot(sps_spatial)
# ggplot(data = sps, aes(x = lon, y = lat)) +
#   borders() + geom_point(colour = "red") +
#   coord_quickmap(expand = F) + theme_void()

# Split into training and testing
# train <- sample(1:nrow(sps), 0.7*nrow(sps), replace = FALSE)
# sps_train <- sps[train,]
# sps_test <- sps[-train,]

# Wrapper function for loading an ASCII file as a data.frame
# load_asc_to_df <- function(file_name){
#   res <- as.data.frame(read.asciigrid(file_name), xy = T)
# }

# Load the present variables
## Need to screen the files loaded based on the species being modelled
# expl <- load_asc_to_df(var_files[1]) %>% 
#   dplyr::rename(lon = s1, lat = s2)
# expl <- map_dfc(var_files[1:2], load_asc_to_df)
# expl <- raster::stack(var_files)
expl <- raster::stack("data/present/Bottom.Temp.max.asc",
                      "data/present/SST.min.asc",
                      "data/present/landdist_clip.asc")
# plot(expl)


# 3: Prep data ------------------------------------------------------------

# Prep data for modelling
biomod_data <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(sps)),
  resp.xy = as.matrix(sps[,2:3]),
  # resp.var = sps_spatial,
  resp.name = sps$sps[1],
  # eval.resp.var = rep(1, nrow(sps_test)), # Doesn't work with presence only...
  # eval.resp.xy = as.matrix(sps_test[,2:3]),
  expl.var = expl, 
  PA.nb.rep = 5,
  PA.strategy = "sre",
  PA.sre.quant = 0.1)

# check data format
biomod_data

# Check plot of data
plot(biomod_data)

# Model options
biomod_option <- BIOMOD_ModelingOptions()

# Cross validation
biomod_cv <- BIOMOD_cv(biomod_data)

# Need to look into this
# biomod_tuning <- BIOMOD_tuning(biomod_data, # Doesn't run... 
#                                env.ME = expl,
#                                n.bg.ME = ncell(expl))

# 4: Model ----------------------------------------------------------------

# Run the model
biomod_model <- BIOMOD_Modeling(
  biomod_data,
  models = c('GLM', 'GAM', 'ANN', 'SRE', 'RF'),
  models.options = biomod_option,
  NbRunEval = 3,
  DataSplit = 70,
  # Yweights = NULL,
  VarImport = 1,
  models.eval.meth = c('TSS','ROC', 'KAPPA', 'ACCURACY', 'BIAS'),
  # SaveObj = FALSE,
  rescal.all.models = TRUE,
  do.full.models = F)
# save(model_out, file = paste0("data/biomod_",sps$sps[1],".Rdata"))
# load(paste0("data/biomod_",sps$sps[1],".Rdata"))

# Have a look at the outputs
biomod_model

# Relative importance of exploratory variables
variable_importances <- get_variables_importance(biomod_model)
variable_importances

# Get all models evaluation
# evaluate()
biomod_eval <- get_evaluations(biomod_model)
dimnames(biomod_eval)

biomod_eval[,,,"RUN1","PA5"]

# Visualise quality of different models
models_scores_graph(biomod_model)


# 5. Present projections --------------------------------------------------

# prints out whichever other variable you need (See dimensions to pick which evaluations) 
# best: model_eval["ROC",,,"Full","PA1"]

biomod_model@models.computed

# For ensemble forecast from all models
# BIOMOD_EnsembleForecasting()

# For presence only
# BIOMOD_presenceonly()

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
  build.clamping.mask = FALSE)#,
  # SaveObj = TRUE)
# save(projection, file = "data/projection_Aebu.Rdata")
# load("data/projection_Aebu.Rdata")

# Plot projections
plot(biomod_projection)

# Get predictions
  # NB: There are many biomod2::get_ functions for looking at results more closely
present_predictions <- get_predictions(biomod_projection)
present_predictions

# Look at particular aspects of predictions
biomod2::free()

# present_GAM <- raster(present_projection, layer = "ECKMAX_PA1_Full_GAM") # doesn't work
# present_MAXENT <- raster(present_projection, layer = "Aebu_PA1_Full_MAXENT.Phillips")
present_GLM <- raster(present_projection, layer = "Aebu_PA1_Full_GLM")
present_GBM <- raster(present_projection, layer = "Aebu_PA1_Full_GBM")
present_MARS <- raster(present_projection, layer = "Aebu_PA1_Full_MARS")
present_CTA <- raster(present_projection, layer = "Aebu_PA1_Full_CTA")
present_RF <- raster(present_projection, layer = "Aebu_PA1_Full_RF")

# Save individual projections
writeRaster(present_GLM, filename = 'data/ECKMAX/proj_GGM/GLM.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(present_GBM, filename = 'data/ECKMAX/proj_GGM/GBM.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(present_MARS, filename = 'data/ECKMAX/proj_GGM/MARS.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(present_CTA, filename = 'data/ECKMAX/proj_GGM/CTA.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(present_RF, filename = 'data/ECKMAX/proj_GGM/RF.asc',
            format = "ascii", overwrite = TRUE )


# 6: Future projections ---------------------------------------------------

# Currently only looking at 50 year projections
expl_50 <- raster::stack(
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/future/2050/Bottom.Temp.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/future/2050/SST.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/future/2050/Bottom.Salinity.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/future/2050/SSS.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/future/2050/Ice.thick.Mean.asc"
)
plot(expl_50)

# Run 50 in situ projections
projection_50 <- BIOMOD_Projection(
  modeling.output = model_out,
  new.env = expl_50,
  proj.name = 'GGM50',
  xy.new.env = as.matrix(Arctic_BO[,1:2]),
  selected.models = model_out@models.computed[13:16],
  Bin.trans = TRUE,
  slot = model_out@models.computed,
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  SaveObj = TRUE)
save(projection_50, file = "data/projection_Aebu_50.Rdata")
load("data/projection_Aebu_50.Rdata")

# Plot projections
plot(projection_50)

# Get projections
current_projection_50 <- get_predictions(projection_50)
current_projection_50

# Somehow use this to compare the projections over time
BIOMOD_RangeSize()

# current_GAM <- raster(current_projection_50, layer = "ECKMAX_PA1_Full_GAM") # doesn't work
# current_MAXENT <- raster(current_projection_50, layer = "ECKMAX_PA1_Full_MAXENT.Tsuruoka")
current_GLM <- raster(current_projection_50, layer = "Aebu_PA1_Full_GLM")
current_GBM <- raster(current_projection_50, layer = "Aebu_PA1_Full_GBM")
current_MARS <- raster(current_projection_50, layer = "Aebu_PA1_Full_MARS")
current_CTA <- raster(current_projection_50, layer = "Aebu_PA1_Full_CTA")
current_RF <- raster(current_projection_50, layer = "Aebu_PA1_Full_RF")

# Save individual projections
writeRaster(current_GLM, filename = 'data/ECKMAX/proj_GGM_50/GLM.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_GBM, filename = 'data/ECKMAX/proj_GGM_50/GBM.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_MARS, filename = 'data/ECKMAX/proj_GGM_50/MARS.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_CTA, filename = 'data/ECKMAX/proj_GGM_50/CTA.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_RF, filename = 'data/ECKMAX/proj_GGM_50/RF.asc',
            format = "ascii", overwrite = TRUE )

