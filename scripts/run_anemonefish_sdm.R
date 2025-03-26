# scripts/run_anemonefish_sdm.R
# This script runs the SDM for all anemonefish species.

library(raster)
library(sf)

# Define paths
output_folder <- "data/sdm_output"  # Create this folder if it doesn't exist

# Load the environmental stack (if not already in memory)
# env_stack <- readRDS("data/env/selected_env_stack.rds") #Uncomment if you saved

if (!exists("env_stack")) {
  stop("env_stack not found. Please run 05_1_env_variable_selection.R first.")
}

# Load species lists
anemonefish_list <- read.csv("final_anemonefish_species_list.csv")

#Run the main SDM script
source("scripts/05_5_run_sdm.R")

# Run SDM for all anemonefish species
print("Starting SDMs for anemonefish...")
for (i in 1:nrow(anemonefish_list)) {
  species_name <- anemonefish_list$scientificName[i]
  species_aphia_id <- anemonefish_list$AphiaID[i]
  
  #Retrieve the thinned background occurances.
  background_points <- background_points_list[[gsub(" ", "_", species_name)]]
  
  # Assign thinned background occurances to global object for SDM steps.
  assign(paste0(species_name,"_background_points"), background_points, envir = .GlobalEnv)
  
  print(paste("Starting SDM for anemonefish:", species_name, " (", i, "/", nrow(anemonefish_list), ")"))
  
  final_prediction_raster <- run_sdm(species_name, output_folder, env_stack)
  
  if (is.null(final_prediction_raster)) {
    print(paste("Skipping to the next anemonefish species."))
  }
}
print("SDM processing completed for all anemonefish species.")