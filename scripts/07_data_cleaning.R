# scripts/07_data_cleaning.R
# This script cleans the occurrence data (removes NAs, converts to sf)
# and performs spatial thinning.

library(dplyr)
library(sf)
library(raster)

# Source the helper functions
source("helpers/sdm_helpers.R")

# Load the environmental stack (if not already in memory)
# env_stack <- readRDS("data/env/selected_env_stack.rds") #Uncomment if you saved

if (!exists("env_stack")) {
  stop("env_stack not found. Please run 05_1_env_variable_selection.R first.")
}

# Function to process species data
process_species_data <- function(species_list_file, occurrence_folder, output_folder, env_stack) {
  # Read the species list
  species_list <- read.csv(species_list_file)
  
  # Iterate through each species
  for (i in 1:nrow(species_list)) {
    species_name <- species_list$scientificName[i]
    species_aphia_id <- species_list$AphiaID[i]
    cat("\n","-----------------------------","\n")
    cat("-----------------------------","\n")
    cat("\n","Processing species:", species_name, "(", i, "/", nrow(species_list), ")","\n")
    cat("-----------------------------","\n")
    cat("-----------------------------","\n")
    
    # Construct the occurrence file path
    occurrence_file <- file.path(occurrence_folder, paste0(species_aphia_id, ".csv"))
    
    # Check if the occurrence file exists
    if (!file.exists(occurrence_file)) {
      warning(paste("No occurrence file found for species:", species_name, "Aphia ID:", species_aphia_id))
      next # Skip to the next species
    }
    
    # Load occurrence data
    occurrences <- read.csv(occurrence_file)
    
    # Clean occurrence data
    occurrences_cleaned <- clean_occurrence_data(occurrences)
    
    # Check for valid occurrences after cleaning
    if (nrow(occurrences_cleaned) == 0) {
      warning(paste("No valid occurrences after cleaning for species:", species_name, "Aphia ID:", species_aphia_id))
      next # Skip to the next species
    }
    
    # Spatial thinning
    occurrences_thinned <- thin_occurrence_data(occurrences_cleaned, env_stack)
    
    # Check for valid occurrences after thinning
    if (nrow(occurrences_thinned) == 0) {
      warning(paste("No valid occurrences after thinning for species:", species_name, "Aphia ID:", species_aphia_id))
      next # Skip to the next species
    }
    
    # Optionally save the cleaned and thinned occurrences
    # output_file <- file.path(output_folder, paste0(gsub(" ", "_", species_name), "_thinned.rds"))
    # saveRDS(occurrences_thinned, file = output_file)
    # cat("\n","Saved thinned occurrences to:", output_file,"\n")
    #-------------------------------------------------------------------------------
    # Assign thinned occurances to global object for the remaining scripts
    #-------------------------------------------------------------------------------
    
    assign(paste0(gsub(" ", "_", species_name), "_thinned"), occurrences_thinned, envir = .GlobalEnv)
    #-------------------------------------------------------------------------------
  }
}

# Define species lists and folder paths
species_lists <- list(
  anemonefish = list(species_list_file = "final_anemonefish_species_list.csv",
                     occurrence_folder = "data/occurrence/anemonefish"),
  anemone = list(species_list_file = "final_anemone_species_list.csv",
                 occurrence_folder = "data/occurrence/anemone")
)
# Perform data cleaning and thinning for both anemonefish and anemone species
for (group in names(species_lists)) {
  process_species_data(species_lists[[group]]$species_list_file,
                       species_lists[[group]]$occurrence_folder,
                       "data/sdm_output", # Output folder for saving results
                       env_stack)
}