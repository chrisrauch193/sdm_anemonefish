# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for running individual species SDMs using ENMeval
# Adapts logic from user's original scripts (05_data_standardisation, 07, 08, 09)
#-------------------------------------------------------------------------------
pacman::p_load(dplyr, sf, terra, dismo, ENMeval, readr, tools) # Ensure needed packages loaded

#' Load and Clean Individual Species Occurrence Data
#'
#' Loads a single species CSV, cleans coordinates, converts to sf.
#'
#' @param species_aphia_id The AphiaID of the species.
#' @param occurrence_folder Path to the directory containing the species CSV file.
#' @param config Configuration list (for occurrence_crs).
#'
#' @return An sf object with cleaned occurrence points, or NULL if file not found or no valid points.
#' @export
load_clean_individual_occ <- function(species_aphia_id, occurrence_folder, config) {
  
  occurrence_file <- file.path(occurrence_folder, paste0(species_aphia_id, ".csv"))
  
  if (!file.exists(occurrence_file)) {
    warning("Occurrence file not found for AphiaID ", species_aphia_id, " in ", occurrence_folder, call. = FALSE)
    return(NULL)
  }
  
  tryCatch({
    occurrences <- readr::read_csv(occurrence_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  }, error = function(e){
    warning("Error reading occurrence file ", occurrence_file, ": ", e$message, call.=FALSE)
    return(NULL)
  })
  
  if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occurrences))) {
    warning("Required columns 'decimalLongitude'/'decimalLatitude' not found in ", occurrence_file, call.=FALSE)
    return(NULL)
  }
  
  # Clean NAs and invalid coords
  occ_clean <- occurrences %>%
    mutate(
      decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
      decimalLatitude = suppressWarnings(as.numeric(decimalLatitude))
    ) %>%
    filter(
      !is.na(decimalLongitude), !is.na(decimalLatitude),
      decimalLongitude >= -180, decimalLongitude <= 180,
      decimalLatitude >= -90, decimalLatitude <= 90
    )
  
  if (nrow(occ_clean) == 0) {
    # cat("  No valid coordinates after cleaning for AphiaID:", species_aphia_id, "\n")
    return(NULL)
  }
  
  # Convert to sf
  tryCatch({
    occ_sf <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs)
    return(occ_sf)
  }, error = function(e){
    warning("Error converting occurrences to sf for AphiaID ", species_aphia_id, ": ", e$message, call.=FALSE)
    return(NULL)
  })
}


#' Thin Individual Species Occurrence Data (Cell-based)
#'
#' Performs spatial thinning using the cell-based method from the original script.
#' Requires the environmental raster stack (predictors) for the target scenario.
#'
#' @param occurrences_sf sf object of cleaned points for a single species.
#' @param predictor_stack SpatRaster object (e.g., PCA layers) for the target scenario.
#' @param config Configuration list.
#'
#' @return An sf object with thinned occurrence points, or NULL if input invalid or thinning results in no points.
#' @export
thin_individual_occ <- function(occurrences_sf, predictor_stack, config) {
  
  if (is.null(occurrences_sf) || nrow(occurrences_sf) == 0 || is.null(predictor_stack)) {
    return(NULL)
  }
  
  # Ensure CRS match before thinning/extracting
  target_crs_obj <- terra::crs(predictor_stack)
  if (sf::st_crs(occurrences_sf) != target_crs_obj) {
    # cat("    Transforming occurrences for thinning...\n") # Can be verbose
    tryCatch({
      occurrences_sf <- sf::st_transform(occurrences_sf, crs = target_crs_obj)
    }, error = function(e){
      warning("Thinning Error: Failed to transform CRS: ", e$message, call.=FALSE)
      return(NULL) # Cannot proceed if transform fails
    })
  }
  
  # Perform cell-based thinning (using the logic from the group thinning function)
  tryCatch({
    occs_spatvector <- terra::vect(occurrences_sf)
    # Extract cell numbers from the *first layer* of the predictor stack
    occs_cells <- terra::extract(predictor_stack[[1]], occs_spatvector, cells = TRUE)
    
    valid_cells_indices <- which(!is.na(occs_cells$cell))
    if(length(valid_cells_indices) == 0) {
      # warning("No occurrences fall within the predictor raster extent.", call.=FALSE)
      return(NULL)
    }
    
    # Find indices of unique cells among the valid points
    unique_cell_indices_in_valid <- which(!duplicated(occs_cells$cell[valid_cells_indices]))
    # Map these back to the original row indices of the sf object
    original_indices_to_keep <- valid_cells_indices[unique_cell_indices_in_valid]
    
    occs_thinned <- occurrences_sf[original_indices_to_keep, ]
    
    if(nrow(occs_thinned) == 0) return(NULL) # Should not happen if logic correct
    
    return(occs_thinned)
    
  }, error = function(e){
    warning("Error during individual cell-based thinning: ", e$message, call.=FALSE)
    return(NULL) # Return NULL if thinning fails
  })
}


#' Generate Background Points
#'
#' Generates random background points within the extent of the predictor rasters.
#' Uses dismo::randomPoints.
#'
#' @param predictor_stack SpatRaster object (e.g., PCA layers) for the target scenario.
#' @param n Number of background points to generate (from config).
#' @param config Configuration list.
#'
#' @return A data frame with columns 'x' and 'y' for background point coordinates, or NULL on error.
#' @export
generate_sdm_background <- function(predictor_stack, n, config) {
  cat("  Generating", n, "background points...\n")
  if (is.null(predictor_stack)) return(NULL)
  
  tryCatch({
    # randomPoints needs a RasterLayer/Stack/Brick, use first layer of SpatRaster
    # Need to convert SpatRaster to RasterLayer for dismo compatibility
    # Ensure raster package is loaded if needed
    if (!requireNamespace("raster", quietly = TRUE)) pacman::p_load(raster)
    
    # predictor_raster <- raster::raster(predictor_stack[[1]]) # Convert first layer
    # Alternatively, dismo might work directly with SpatRaster now? Test this.
    # bg_points <- dismo::randomPoints(predictor_raster, n = n)
    
    # Safer alternative using terra directly:
    # 1. Create a mask (cells with non-NA values in the first layer)
    mask_rast <- !is.na(predictor_stack[[1]])
    # 2. Sample random points from the valid cells
    # points=TRUE gives coordinates directly
    bg_points_terra <- terra::spatSample(mask_rast, size = n, method = "random", na.rm = TRUE, xy = TRUE, warn=FALSE)
    
    # Check if enough points were generated
    if(nrow(bg_points_terra) < n) {
      warning("Only generated ", nrow(bg_points_terra), " background points (less than requested ", n, ").", call.=FALSE)
    }
    if(nrow(bg_points_terra) == 0) {
      warning("Failed to generate any background points.", call.=FALSE)
      return(NULL)
    }
    
    # Return as data frame with x, y columns as expected by ENMeval
    bg_df <- as.data.frame(bg_points_terra[, c("x", "y")])
    return(bg_df)
    
  }, error = function(e) {
    warning("Error generating background points: ", e$message, call. = FALSE)
    return(NULL)
  })
}


#' Run SDM using ENMeval
#'
#' Executes ENMeval with specified settings. Requires cleaned+thinned occurrences (sf)
#' and background points (data.frame). Assumes predictor_stack uses terra.
#'
#' @param occs_thinned_sf sf object of thinned occurrence points.
#' @param predictor_stack SpatRaster object of predictors (e.g., PCA layers).
#' @param background_df Data frame with 'x', 'y' columns for background points.
#' @param config Configuration list (for SDM settings).
#'
#' @return An ENMevaluation object, or NULL on error.
#' @export
run_sdm_enmeval <- function(occs_thinned_sf, predictor_stack, background_df, config) {
  cat("  Running ENMeval...\n")
  
  if (is.null(occs_thinned_sf) || nrow(occs_thinned_sf) == 0 || is.null(predictor_stack) || is.null(background_df)) {
    warning("Invalid inputs for running ENMeval.", call. = FALSE)
    return(NULL)
  }
  
  # Prepare occurrence coordinates for ENMeval (needs data frame or matrix)
  occs_coords <- sf::st_coordinates(occs_thinned_sf)
  colnames(occs_coords) <- c("x", "y")
  
  # Ensure predictor_stack is a RasterStack/Brick if ENMeval requires it, or if terra works directly
  # Convert SpatRaster to RasterStack for ENMeval if needed
  # Check ENMeval version requirements - newer versions might support terra
  if (!requireNamespace("raster", quietly = TRUE)) pacman::p_load(raster)
  tryCatch({
    predictor_rasterstack <- raster::stack(predictor_stack)
  }, error = function(e){
    warning("Failed to convert predictor SpatRaster to RasterStack for ENMeval. Error: ", e$message, "\n  Check ENMeval compatibility or raster package.", call.=FALSE)
    return(NULL)
  })
  
  
  # Setup parallel processing if enabled
  if (config$use_parallel && config$num_cores > 1) {
    doParallel::registerDoParallel(cores = config$num_cores)
    cat("    Using", config$num_cores, "cores for parallel execution.\n")
    parallel_enabled <- TRUE
  } else {
    parallel_enabled <- FALSE
  }
  
  # Run ENMevaluate
  enmeval_results <- NULL
  tryCatch({
    enmeval_results <- ENMeval::ENMevaluate(
      occs = occs_coords,
      envs = predictor_rasterstack, # Use RasterStack here
      bg = background_df,
      algorithm = config$sdm_algorithm,
      partitions = config$sdm_partitions,
      partition.settings = list(kfolds = config$sdm_n_folds), # Pass kfolds if using kfold
      tune.args = config$sdm_tune_args,
      parallel = parallel_enabled,
      numCores = if (parallel_enabled) config$num_cores else NULL # Specify cores if parallel
      # Other arguments like clamp, etc., can be added if needed
    )
  }, error = function(e) {
    warning("ENMevaluate failed: ", e$message, call. = FALSE)
    enmeval_results <- NULL
  }, finally = {
    # Stop parallel backend if it was started
    if(parallel_enabled) {
      # Check if doParallel is registered before stopping
      # This avoids errors if ENMevaluate failed early
      if(foreach::getDoParRegistered()) {
        doParallel::stopImplicitCluster()
        # foreach::registerDoSEQ() # Optionally register sequential backend
      }
    }
  })
  
  return(enmeval_results)
}


#' Select Best Model and Predict SDM
#'
#' Selects the best model from ENMeval results based on the specified metric
#' and generates a prediction raster.
#'
#' @param enmeval_results The output object from ENMeval::ENMevaluate.
#' @param predictor_stack SpatRaster object used for prediction.
#' @param config Configuration list (for evaluation metric).
#'
#' @return A SpatRaster object with the prediction, or NULL on error or if no best model found.
#' @export
predict_sdm_best <- function(enmeval_results, predictor_stack, config) {
  cat("  Selecting best model and predicting...\n")
  if (is.null(enmeval_results) || !inherits(enmeval_results, "ENMevaluation") || nrow(enmeval_results@results) == 0) {
    warning("Invalid ENMeval results provided.", call. = FALSE)
    return(NULL)
  }
  if (is.null(predictor_stack)) {
    warning("Predictor stack is required for prediction.", call. = FALSE)
    return(NULL)
  }
  
  # Select best model based on metric
  eval_results_df <- enmeval_results@results
  metric <- config$sdm_evaluation_metric
  
  if(!metric %in% names(eval_results_df)) {
    warning("Evaluation metric '", metric, "' not found in ENMeval results. Available: ", paste(names(eval_results_df), collapse=", "), call.=FALSE)
    # Fallback to AICc if available
    if("AICc" %in% names(eval_results_df)) {
      metric <- "AICc"
      cat("    Falling back to AICc for model selection.\n")
    } else {
      warning("AICc also not found. Cannot select best model.", call.=FALSE)
      return(NULL)
    }
  }
  
  # Lower AICc is better, higher AUC/TSS is better
  select_best <- if (metric == "AICc") which.min else which.max
  # Handle potential NA values in the metric column
  valid_rows <- !is.na(eval_results_df[[metric]])
  if(sum(valid_rows) == 0){
    warning("All values for metric '", metric, "' are NA. Cannot select best model.", call.=FALSE)
    return(NULL)
  }
  best_model_index <- select_best(eval_results_df[[metric]][valid_rows])
  # Adjust index if NAs were present
  original_indices <- which(valid_rows)
  best_model_row_index <- original_indices[best_model_index]
  
  
  best_model_params <- eval_results_df[best_model_row_index, ]
  best_model_obj <- enmeval_results@models[[best_model_row_index]] # Get the actual model object
  
  cat("    Best model selected based on", metric, ":", best_model_params$tune.args, "\n")
  # print(best_model_params) # Print details of best model
  
  # Predict using the best model object
  # Ensure raster package loaded for predict generic if needed
  if (!requireNamespace("raster", quietly = TRUE)) pacman::p_load(raster)
  prediction_raster <- NULL
  tryCatch({
    # Use terra::predict directly with the model object
    # Need to ensure the model object is compatible
    # Maxnet models from ENMeval might work? Check documentation/examples
    # If not, might need raster::predict with RasterStack version
    prediction_raster_terra <- terra::predict(
      predictor_stack,
      model = best_model_obj,
      # type = "logistic", # Maxnet usually outputs cloglog or exponential, ENMeval often standardizes to logistic? Check ENMeval docs.
      # Assuming ENMeval provides a model compatible with predict() generic
      # May need specific arguments depending on model type
      na.rm = TRUE
    )
    # Sometimes predict returns multiple layers, assume first is suitability
    prediction_raster <- prediction_raster_terra[[1]]
    names(prediction_raster) <- "suitability"
    
  }, error = function(e){
    warning("Prediction using terra::predict failed: ", e$message, call.=FALSE)
    # Fallback using raster::predict if terra failed
    cat("    Attempting fallback prediction using raster::predict...\n")
    tryCatch({
      predictor_rasterstack <- raster::stack(predictor_stack)
      prediction_raster_raster <- raster::predict(
        object = predictor_rasterstack,
        model = best_model_obj,
        # type = "cloglog", # Adjust type based on algorithm output if needed
        na.rm = TRUE
      )
      # Convert back to SpatRaster
      prediction_raster <- terra::rast(prediction_raster_raster)
      names(prediction_raster) <- "suitability"
    }, error = function(e2){
      warning("Fallback prediction using raster::predict also failed: ", e2$message, call.=FALSE)
      prediction_raster <- NULL
    })
  })
  
  if(!is.null(prediction_raster)) {
    cat("    Prediction raster generated.\n")
  }
  
  return(prediction_raster)
}

#-------------------------------------------------------------------------------