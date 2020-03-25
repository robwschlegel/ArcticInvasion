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

library(tidyverse)
library(biomod2)
library(FNN)

# NB: The data are housed on dropbox on Jesi Goldsmit's professional account
# They have been downloaded locally to the dropbox folder on my machine
# I then created a symbolic link from there to this project folder

# The species occurrence data
sps_files <- dir("data/occurrence", full.names = T)

# The environmental file pathways
var_files <- dir("data/present", full.names = T)


# 2: Load data ------------------------------------------------------------

# Arctic coords from BioOracle
load("~/ArcticKelp/data/Arctic_BO.RData")
ggplot(data = Arctic_BO, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = BO2_templtmin_bdmax))

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
  dplyr::select(sps, s1, s2) %>%
  dplyr::rename(lon = s1, lat = s2)

# Visualise
ggplot(data = sps, aes(x = lon, y = lat)) +
  borders() + geom_point(colour = "red") +
  coord_quickmap(expand = F) + theme_void()

# Split into training and testing
train <- sample(1:nrow(sps), 0.7*nrow(sps), replace = FALSE)
sps_train <- sps[train,]
sps_test <- sps[-train,]

# Wrapper function for loading an ASCII file as a data.frame
load_asc_to_df <- function(file_name){
  res <- as.data.frame(sp::read.asciigrid(file_name), xy = T)
}



# Load the present variables
## Need to screen the files loaded based on the species being modelled
expl <- load_asc_to_df(var_files[1]) %>% 
  dplyr::rename(lon = s1, lat = s2)
# expl <- map_dfc(var_files[1:2], load_asc_to_df)
expl <- raster::stack(var_files)
# plot(expl)

test_join <- left_join(sps, expl)


# 3: Prep data ------------------------------------------------------------

# Prep data for modelling
biomod_data <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(sps_train)),
  resp.xy = as.matrix(sps_train[,2:3]),
  resp.name = sps$sps[1],
  expl.var = expl)

# check data format
biomod_data

# Check plot of data
plot(biomod_data)

# Cross validation
BIOMOD_cv()


# 4: Model ----------------------------------------------------------------

# Model options
biomod_option <- BIOMOD_ModelingOptions(GLM = list(type = 'polynomial', interaction.level = 1))

# Need to lk into this
BIOMOD_tuning()

# Run the model
model_out <- BIOMOD_Modeling( # No need to run a second time
  biomod_data,
  models = c('GLM', 'GBM','CTA','RF'),
  models.options = biomod_option,
  NbRunEval = 3,
  DataSplit = 70,
  Yweights = NULL,
  VarImport = 3,
  models.eval.meth = c('TSS','ROC', 'KAPPA'),
  SaveObj = TRUE,
  rescal.all.models = TRUE)
save(model_out, file = "data/model_out_Aebu.Rdata")
load("data/model_out_Aebu.Rdata")

# Have a look at the outputs
model_out

# Relative importance of exploratory variables
variable_importances <- get_variables_importance(model_out)
variable_importances

# Get all models evaluation
evaluate()
model_eval <- get_evaluations(model_out)
dimnames(model_eval)

model_eval[,,,"Full",]

# Visualise quality of different models
biomod2::models_scores_graph()


# 5. Present projections --------------------------------------------------

# prints out whichever other variable you need (See dimensions to pick which evaluations) 
# best: model_eval["ROC",,,"Full","PA1"]

model_out@models.computed

# For ensemble forecast from all models
BIOMOD_EnsembleForecasting()

# For presence only
BIOMOD_presenceonly()

# Create projections
projection <- BIOMOD_Projection(
  modeling.output = model_out,
  new.env = expl,
  proj.name = 'GGM',
  xy.new.env = as.matrix(Arctic_BO[,1:2]),
  selected.models = model_out@models.computed[13:16],
  Bin.trans = TRUE,
  slot = model_out@models.computed,
  binary.meth ='TSS',
  compress = 'xz',
  clamping.mask = F,
  SaveObj = TRUE)
save(projection, file = "data/projection_Aebu.Rdata")
load("data/projection_Aebu.Rdata")

# Plot projections
plot(projection)

# Get projections
  # NB: There are many biomod2::get_ functions for looking at results more closely
current_projection <- get_predictions(projection)
current_projection

# Look at particular aspects of projections
biomod2::free()

# current_GAM <- raster(current_projection, layer = "ECKMAX_PA1_Full_GAM") # doesn't work
# current_MAXENT <- raster(current_projection, layer = "Aebu_PA1_Full_MAXENT.Phillips")
current_GLM <- raster(current_projection, layer = "Aebu_PA1_Full_GLM")
current_GBM <- raster(current_projection, layer = "Aebu_PA1_Full_GBM")
current_MARS <- raster(current_projection, layer = "Aebu_PA1_Full_MARS")
current_CTA <- raster(current_projection, layer = "Aebu_PA1_Full_CTA")
current_RF <- raster(current_projection, layer = "Aebu_PA1_Full_RF")

# Save individual projections
writeRaster(current_GLM, filename = 'data/ECKMAX/proj_GGM/GLM.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_GBM, filename = 'data/ECKMAX/proj_GGM/GBM.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_MARS, filename = 'data/ECKMAX/proj_GGM/MARS.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_CTA, filename = 'data/ECKMAX/proj_GGM/CTA.asc',
            format = "ascii", overwrite = TRUE )

writeRaster(current_RF, filename = 'data/ECKMAX/proj_GGM/RF.asc',
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

