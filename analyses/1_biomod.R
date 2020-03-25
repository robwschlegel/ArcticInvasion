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


# 2: Load data ------------------------------------------------------------

# Arctic coords
load("~/ArcticKelp/data/Arctic_BO.RData")
ggplot(data = Arctic_BO, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = ))

# NB: The data are housed on dropbox on Jesi Goldsmit's professional account
# They have been downloaded locally to the dropbox folder on my machine

# Load a test species
Aebu <- read_csv("~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/occurrence/Aebu_near.csv")

# Load a few test variables
expl <- raster::stack(
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/present/Bottom.Temp.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/present/SST.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/present/Bottom.Salinity.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/present/SSS.Mean.asc",
  "~/modeling team! Dropbox/hotspot modelling/modelling_data/environmental/present/Ice.thick.Mean.asc"
)
plot(expl)


# 3: Prep data ------------------------------------------------------------

# Create pseudoabsence points
  # NB: Need to verify these argument choices
biomod_data <- BIOMOD_FormatingData(
  resp.var = as.matrix(rep(1, nrow(Aebu))),
  expl.var = expl,
  resp.xy = as.matrix(Aebu[,2:3]),
  resp.name = Aebu$sps[1],
  PA.nb.rep = 1,
  PA.nb.absences = 2100,
  PA.strategy = 'sre',
  PA.dist.min = 1,
  PA.dist.max = NULL,
  PA.sre.quant = 0.01)

# check data format
biomod_data

# Check plot of data
plot(biomod_data)


# 4: Model ----------------------------------------------------------------

# Model options
biomod_option <- BIOMOD_ModelingOptions(GLM = list(type = 'polynomial', interaction.level = 1))

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
model_eval <- get_evaluations(model_out)
dimnames(model_eval)

model_eval[,,,"Full",]


# 5. Present projections --------------------------------------------------

# prints out whichever other variable you need (See dimensions to pick which evaluations) 
# best: model_eval["ROC",,,"Full","PA1"]

model_out@models.computed

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
current_projection <- get_predictions(projection)
current_projection

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

