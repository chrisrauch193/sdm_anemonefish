# scripts/08_background_generation.R
# This script generates background points for each species.

library(raster)
library(sf)

# Source the helper functions
source("helpers/sdm_helpers.R")

# Load the environmental stack (if not already in memory)
# env_stack <- readRDS("data/env/selected_env_stack.rds") #Uncomment if you saved

if (!exists("env_stack")) {
  stop("env_stack not found. Please run 05_1_env_variable_selection.R first.")
}

# List to store background points for each species
background_points_list <- list()

# Function to generate background points for a single species
generate_species_background <- function(species_name, env_stack, n = 10000) {
  
  # Get the thinned occurrences from global assignment in script 05_3_data_cleaning.R
  occurrences_thinned <- get(paste0(gsub(" ", "_", species_name), "_thinned"))
  
  # Generate background points
  background_points <- generate_background_points(occurrences_thinned, env_stack, n = n)
  cat("\n","Generated background points for species:", species_name,"\n")
  return(background_points)
}

# Load species lists
anemonefish_list <- read.csv("final_anemonefish_species_list.csv")
anemone_list <- read.csv("final_anemone_species_list.csv")

# Generate background points for anemonefish species
for (i in 1:nrow(anemonefish_list)) {
  species_name <- anemonefish_list$scientificName[i]
  background_points <- generate_species_background(species_name, env_stack)
  
  # Store background points in the list
  background_points_list[[gsub(" ", "_", species_name)]] <- background_points
  cat("\n","Stored background points for species:", species_name,"\n")
  
}

# Generate background points for anemone species
for (i in 1:nrow(anemone_list)) {
  species_name <- anemone_list$scientificName[i]
  background_points <- generate_species_background(species_name, env_stack)
  
  # Store background points in the list
  background_points_list[[gsub(" ", "_", species_name)]] <- background_points
  cat("\n","Stored background points for species:", species_name,"\n")
}

# Optionally, save the background points list
saveRDS(background_points_list, file = "data/sdm_output/background_points_list.rds")