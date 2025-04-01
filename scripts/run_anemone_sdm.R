# scripts/run_anemone_sdm.R
# This script runs the SDM for all anemone species.

library(raster)
library(sf)

# Define paths
output_folder <- "data/sdm_output"  # Create this folder if it doesn't exist

# Load the environmental stack (if not already in memory)
# env_stack <- readRDS("data/env/selected_env_stack.rds") #Uncomment if you saved

if (!exists("env_stack")) {
  stop("env_stack not found. Please run 05_env_variable_selection.R first.")
}

# Load species lists
anemone_list <- read.csv("final_anemone_species_list.csv")

#Run the main SDM script
source("scripts/09_run_sdm.R")

# Run SDM for all anemone species
print("Starting SDMs for anemones...")
for (i in 1:nrow(anemone_list)) {
  species_name <- anemone_list$scientificName[i]
  species_aphia_id <- anemone_list$AphiaID[i]
  
  #Retrieve the thinned background occurances.
  background_points <- background_points_list[[gsub(" ", "_", species_name)]]
  
  # Assign thinned background occurances to global object for SDM steps.
  assign(paste0(species_name,"_background_points"), background_points, envir = .GlobalEnv)
  
  print(paste("Starting SDM for anemone:", species_name, " (", i, "/", nrow(anemone_list), ")"))
  
  final_prediction_raster <- run_sdm(species_name, output_folder, env_stack)
  
  if (is.null(final_prediction_raster)) {
    print(paste("Skipping to the next anemone species."))
  }
}
print("SDM processing completed for all anemone species.")