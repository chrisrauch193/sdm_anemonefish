# helpers/extract_env.R

# Helper functions for data extraction and environmental data loading.

library(raster)
library(terra)
library(dplyr)
library(tools)  # file_path_sans_ext

#' @title Load and prepare environmental data
#' @description Loads raster files from a specified directory and performs
#'   preprocessing steps, including cropping and masking to a defined region.
#' @param env_folder Path to the directory containing environmental raster files.
#' @param coral_shapefile Path to the shapefile defining the analysis region.
#' @return A terra::SpatRaster object.
load_and_prepare_env_data <- function(env_folder, coral_shapefile) {
  tryCatch({
    # Load raster files
    raster_files <- list.files(env_folder, pattern = "\\.tif$",
                               full.names = TRUE, recursive = TRUE)
    env_stack    <- raster::stack(raster_files)  # raster::stack
    env_rast     <- terra::rast(env_stack)        # terra::rast
    
    cat("Available environmental variables:\n")
    print(names(env_rast))
    
    # Load coral reef shapefile and apply spatial restrictions
    coral_areas <- terra::vect(coral_shapefile)
    env_rast    <- terra::crop(env_rast, coral_areas)
    env_rast    <- terra::mask(env_rast, coral_areas)
    
    return(env_rast)
  }, error = function(e) {
    cat("\nError loading/preparing environmental data:", e$message, "\n")
    return(NULL)
  })
}

#' @title Load occurrence data and extract environmental values
#' @description Loads occurrence data from CSV files, extracts environmental
#'   values at occurrence locations, and combines the results.
#' @param occurrence_folder Path to the directory containing occurrence data files.
#' @param env_rast SpatRaster object containing the environmental data.
#' @param occurrence_crs Coordinate reference system of the occurrence data.
#' @return A data frame containing the extracted environmental values.
load_occurrence_and_extract_env <- function(occurrence_folder, env_rast,
                                            occurrence_crs) {
  # List occurrence files
  occurrence_files <- list.files(occurrence_folder, pattern = "\\.csv$",
                                 full.names = TRUE)
  
  # Create an empty list to store the extracted environmental values
  all_env_values <- list()
  
  # Loop through each occurrence file
  for (file in occurrence_files) {
    tryCatch({
      # Read the occurrence data
      occurrences <- read.csv(file)
      
      # Check for longitude and latitude columns
      if (!("decimalLongitude" %in% colnames(occurrences) &&
            "decimalLatitude" %in% colnames(occurrences))) {
        stop(paste("Longitude/latitude columns missing in", file))
      }
      
      # Filter out rows with NA in longitude or latitude
      occurrences <- occurrences %>%
        filter(!is.na(decimalLongitude) & !is.na(decimalLatitude))
      
      # Create a SpatVector from the occurrence data
      occurrences_spatvector <- terra::vect(occurrences,
                                            geom = c("decimalLongitude",
                                                     "decimalLatitude"),
                                            crs  = occurrence_crs)
      
      # Reproject if necessary
      raster_crs <- terra::crs(env_rast)
      if (terra::crs(occurrences_spatvector) != raster_crs) {
        occurrences_spatvector <- terra::project(occurrences_spatvector,
                                                 raster_crs)
      }
      
      # Extract environmental values
      env_values <- terra::extract(env_rast, occurrences_spatvector)
      env_values <- env_values[, -1, drop = FALSE] # drop = FALSE to keep it a dataframe even if only 1 column
      
      # Remove NA rows
      env_values <- na.omit(env_values)
      
      #Add the file name as the name
      filename <- tools::file_path_sans_ext(basename(file))
      all_env_values[[filename]] <- env_values
    }, error = function(e) {
      cat("\nError processing", file, ":", e$message, "\n")
    })
  }
  
  # Combine all environmental values into one dataframe
  env_values_df <- do.call(rbind, all_env_values) # Combine dataframes in the list
  
  return(env_values_df)
}