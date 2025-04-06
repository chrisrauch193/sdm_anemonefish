# scripts/sdm_runs/anemone_sdm.R

library(ENMeval)
library(raster)
library(dplyr)
library(sf)
library(tidyr)
library(terra)
library(readxl)

# 1. Data Loading and Setup

# Assuming you have already run scripts 02, 03, and 04

# Load species lists
anemone_list <- read.csv("final_anemone_species_list.csv")

# Define paths
env_folder <- "data/env/current"
output_folder <- "data/sdm_output"  # Create this folder if it doesn't exist

# Create output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Load environmental data (current conditions) - Load this ONCE, not in the loop
# List all raster files in the env_folder
raster_files <- list.files(env_folder, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# Stack the raster files into a single RasterStack
env_stack <- raster::stack(raster_files)

# Print the names of the layers in the RasterStack
print(names(env_stack))

crs <- st_crs(env_stack)  # Get CRS *before* functions, use in functions


# 2. Data Cleaning Function

clean_occurrence_data <- function(occurrences) {
  # Remove records with missing coordinates
  occurrences <- occurrences %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude))
  
  # Convert to sf object
  occurrences_sf <- st_as_sf(occurrences, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) # Assumes WGS84, change if needed
  occurrences_sf <- st_transform(occurrences_sf, crs = st_crs(env_stack))
  return(occurrences_sf)
}

# 3. Spatial Thinning Function

thin_occurrence_data <- function(occurrences_sf, env_stack) {
  # Extract cell numbers for each occurrence
  occs_cells <- raster::extract(env_stack[[1]], st_coordinates(occurrences_sf), cellnumbers = TRUE)
  
  # Identify duplicated cell numbers
  occs_cellDups <- duplicated(occs_cells[, 1])
  
  # Remove duplicated occurrences
  occurrences_thinned <- occurrences_sf[!occs_cellDups, ]
  
  return(occurrences_thinned)
}

# 4. Background Point Generation Function
generate_background_points <- function(occurrences_thinned, env.rast, n = 10000) {
  
  # Generate random background points
  background_points <- dismo::randomPoints(env.rast, n = n) %>%
    as.data.frame()
  colnames(background_points) <- c("x", "y")
  
  
  return(background_points)
}


# 5. SDM with ENMeval Function

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

# 6. Prediction Function

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

# 7. Putting it all together as a single function

run_sdm <- function(species_name, species_aphia_id, species_folder, env_stack, env_folder, output_folder) {
  tryCatch({  # Use tryCatch to handle errors and continue to the next species
    
    # 1. Load occurrence data
    occurrence_file <- list.files(species_folder, pattern = paste0(species_aphia_id, ".csv"), full.names = TRUE)
    if (length(occurrence_file) == 0) {
      warning(paste("No occurrence file found for species:", species_name, "Aphia ID:", species_aphia_id, "in folder:", species_folder))
      return(NULL)  # Skip to the next species
    }
    occurrences <- read.csv(occurrence_file)
    
    # 2. Clean occurrence data
    occurrences_cleaned <- clean_occurrence_data(occurrences)
    
    #check to see if we have any valid occurances
    if(nrow(occurrences_cleaned) == 0) {
      warning(paste("No valid occurances after cleaning for species:", species_name, "Aphia ID:", species_aphia_id))
      return(NULL)  # Skip to the next species
    }
    
    # 3. Thin occurrence data
    occurrences_thinned <- thin_occurrence_data(occurrences_cleaned, env_stack)
    
    #check to see if we have any valid occurances
    if(nrow(occurrences_thinned) == 0) {
      warning(paste("No valid occurances after thinning for species:", species_name, "Aphia ID:", species_aphia_id))
      return(NULL)  # Skip to the next species
    }
    
    # 4. Generate background points
    background_points <- generate_background_points(occurrences_thinned, env_stack)
    
    # 5. Run ENMeval
    enmeval_results <- run_enmeval(occurrences_thinned, env_stack, background_points, species_name)
    
    # 6. Generate prediction raster
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
    # Print error message, Aphia ID, and species name, and then return NULL
    warning(paste("Error running SDM for species:", species_name, "Aphia ID:", species_aphia_id, ":", e$message))
    return(NULL)  # Return NULL so the loop can continue
  })
}




# 8. Run SDM for all anemone species
species_folder <- "data/occurrence/anemone"
print("Starting SDMs for anemones...")
for (i in 1:nrow(anemone_list)) {
  species_name <- anemone_list$scientificName[i]
  species_aphia_id <- anemone_list$AphiaID[i]
  
  print(paste("Starting SDM for anemone:", species_name, " (", i, "/", nrow(anemone_list), ")"))
  
  final_prediction_raster <- run_sdm(species_name, species_aphia_id, species_folder, env_stack, env_folder, output_folder)
  
  if (is.null(final_prediction_raster)) {
    print(paste("Skipping to the next anemone species."))
  }
}
print("SDM processing completed for all anemone species.")