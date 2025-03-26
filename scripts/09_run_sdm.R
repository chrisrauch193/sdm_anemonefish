# scripts/09_run_sdm.R
# This script defines the core SDM functions using ENMeval.  It DOES NOT run the models.
# That is handled in the anemone and anemonefish specific scripts

library(ENMeval)
library(raster)
library(dplyr)
library(sf)

# Source the helper functions
source("helpers/sdm_helpers.R")

# Load the environmental stack (if not already in memory)
# env_stack <- readRDS("data/env/selected_env_stack.rds") #Uncomment if you saved

if (!exists("env_stack")) {
  stop("env_stack not found. Please run 05_1_env_variable_selection.R first.")
}

# Function to run ENMeval
run_enmeval <- function(occurrences_thinned, env_stack, background_points, species_name) {
  # Prepare occurrence and background data for ENMeval
  occs <- st_coordinates(occurrences_thinned)
  colnames(occs) <- c("x", "y") # Necessary for ENMeval
  
  # Run ENMeval
  enmeval_results <- ENMevaluate(
    occs = occs,
    envs = env_stack,
    bg = background_points,
    algorithm = "maxnet",
    partitions = "randomkfold",
    tune.args = list(fc = c("L", "LQ", "H", "LQH"), rm = 1:5),  # Example tuning parameters
    parallel = FALSE # Set to TRUE if you have a parallel backend configured
  )
  
  return(enmeval_results)
}

# Function to generate a prediction raster
generate_prediction_raster <- function(enmeval_results, env_stack) {
  # Select the best model based on AICc
  best_model_index <- which.min(enmeval_results@results$AICc)
  best_model <- enmeval_results@models[[best_model_index]]
  best_tune.args <- enmeval_results@results$tune.args[best_model_index]
  
  # Generate the prediction raster
  prediction_raster <- raster::predict(
    object = env_stack,
    model = best_model,
    type = "logistic"  # Or "cloglog" if appropriate for Maxent
  )
  
  return(prediction_raster)
}

# Putting it all together as a single function
run_sdm <- function(species_name, output_folder, env_stack) {
  tryCatch({  # Use tryCatch to handle errors and continue to the next species
    
    # 1. Get the thinned occurrences and background points
    occurrences_thinned <- get(paste0(gsub(" ", "_", species_name), "_thinned"))
    background_points <- get(paste0(species_name,"_background_points")) #background_points_list[[gsub(" ", "_", species_name)]]
    
    # 2. Run ENMeval
    enmeval_results <- run_enmeval(occurrences_thinned, env_stack, background_points, species_name)
    
    # 3. Generate prediction raster
    prediction_raster <- generate_prediction_raster(enmeval_results, env_stack)
    
    # Save the prediction raster
    output_filename <- file.path(output_folder, paste0(gsub(" ", "_", species_name), "_prediction.tif"))
    raster::writeRaster(prediction_raster, filename = output_filename, overwrite = TRUE)
    
    # Plot and save the plot
    plot_filename <- file.path(output_folder, paste0(gsub(" ", "_", species_name), "_prediction.png"))
    png(filename = plot_filename, width = 800, height = 600)
    plot(prediction_raster, main = paste("SDM Prediction for", species_name))
    points(st_coordinates(occurrences_thinned), pch = 16, col = "red")
    dev.off()
    
    print(paste("SDM completed and saved for species:", species_name))
    
    return(prediction_raster)
    
  }, error = function(e) {
    # Print error message and species name, and then return NULL
    warning(paste("Error running SDM for species:", species_name, ":", e$message))
    return(NULL)  # Return NULL so the loop can continue
  })
}