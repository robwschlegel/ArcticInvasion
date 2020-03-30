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


# 2: load biomod ----------------------------------------------------------


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


# 5: Process --------------------------------------------------------------


# 6: Saved ----------------------------------------------------------------


