# helpers/sdm_helpers.R
# Helper functions for SDM workflow

library(sf)
library(raster)

# Data Cleaning Function
clean_occurrence_data <- function(occurrences) {
  # Remove records with missing coordinates
  occurrences <- occurrences %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude))
  
  # Convert to sf object
  occurrences_sf <- st_as_sf(occurrences, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) # Assumes WGS84, change if needed
  #occurrences_sf <- st_transform(occurrences_sf, crs = st_crs(env_stack)) #moved to later script to avoid missing reference stack
  return(occurrences_sf)
}

# Spatial Thinning Function
thin_occurrence_data <- function(occurrences_sf, env_stack) {
  # Extract cell numbers for each occurrence
  occs_cells <- raster::extract(env_stack[[1]], st_coordinates(occurrences_sf), cellnumbers = TRUE)
  
  # Identify duplicated cell numbers
  occs_cellDups <- duplicated(occs_cells[, 1])
  
  # Remove duplicated occurrences
  occurrences_thinned <- occurrences_sf[!occs_cellDups, ]
  
  return(occurrences_thinned)
}

# Background Point Generation Function
generate_background_points <- function(occurrences_thinned, env.rast, n = 10000) {
  
  # Generate random background points
  background_points <- dismo::randomPoints(env.rast, n = n) %>%
    as.data.frame()
  colnames(background_points) <- c("x", "y")
  
  
  return(background_points)
}