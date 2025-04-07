# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for running individual species SDMs using sdmtune
#-------------------------------------------------------------------------------
pacman::p_load(dplyr, sf, terra, dismo, SDMtune, readr, tools, stringr)

#' Load and Clean Individual Species Occurrence Data
#' (Unchanged from previous version)
load_clean_individual_occ <- function(species_aphia_id, occurrence_folder, config) {
  occurrence_file <- file.path(occurrence_folder, paste0(species_aphia_id, ".csv"))
  if (!file.exists(occurrence_file)) {
    warning("Occurrence file not found for AphiaID ", species_aphia_id, " in ", occurrence_folder, call. = FALSE); return(NULL)
  }
  tryCatch({ occurrences <- readr::read_csv(occurrence_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  }, error = function(e){ warning("Error reading occurrence file ", occurrence_file, ": ", e$message, call.=FALSE); return(NULL) })
  if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occurrences))) {
    warning("Required columns 'decimalLongitude'/'decimalLatitude' not found in ", occurrence_file, call.=FALSE); return(NULL)
  }
  occ_clean <- occurrences %>% mutate(decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)), decimalLatitude = suppressWarnings(as.numeric(decimalLatitude))) %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude), decimalLongitude >= -180, decimalLongitude <= 180, decimalLatitude >= -90, decimalLatitude <= 90)
  if (nrow(occ_clean) == 0) return(NULL)
  tryCatch({ occ_sf <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs); return(occ_sf)
  }, error = function(e){ warning("Error converting occurrences to sf for AphiaID ", species_aphia_id, ": ", e$message, call.=FALSE); return(NULL) })
}

#' Thin Individual Species Occurrence Data (Cell-based)
#' (Unchanged from previous version)
thin_individual_occ <- function(occurrences_sf, predictor_stack, config) {
  if (is.null(occurrences_sf) || nrow(occurrences_sf) == 0 || is.null(predictor_stack)) return(NULL)
  target_crs_obj <- terra::crs(predictor_stack)
  if (sf::st_crs(occurrences_sf) != target_crs_obj) {
    tryCatch({ occurrences_sf <- sf::st_transform(occurrences_sf, crs = target_crs_obj)
    }, error = function(e){ warning("Thinning Error: Failed to transform CRS: ", e$message, call.=FALSE); return(NULL) })
  }
  tryCatch({
    occs_spatvector <- terra::vect(occurrences_sf); occs_cells <- terra::extract(predictor_stack[[1]], occs_spatvector, cells = TRUE)
    valid_cells_indices <- which(!is.na(occs_cells$cell)); if(length(valid_cells_indices) == 0) return(NULL)
    unique_cell_indices_in_valid <- which(!duplicated(occs_cells$cell[valid_cells_indices]))
    original_indices_to_keep <- valid_cells_indices[unique_cell_indices_in_valid]
    occs_thinned <- occurrences_sf[original_indices_to_keep, ]; if(nrow(occs_thinned) == 0) return(NULL)
    return(occs_thinned)
  }, error = function(e){ warning("Error during individual cell-based thinning: ", e$message, call.=FALSE); return(NULL) })
}

#' Generate Background Points
#' (Unchanged from previous version)
generate_sdm_background <- function(predictor_stack, n, config) {
  cat("  Generating", n, "background points...\n")
  if (is.null(predictor_stack)) return(NULL)
  tryCatch({
    mask_rast <- !is.na(predictor_stack[[1]]); bg_points_terra <- terra::spatSample(mask_rast, size = n, method = "random", na.rm = TRUE, xy = TRUE, warn=FALSE)
    if(nrow(bg_points_terra) < n) warning("Only generated ", nrow(bg_points_terra), " background points (less than requested ", n, ").", call.=FALSE)
    if(nrow(bg_points_terra) == 0) { warning("Failed to generate any background points.", call.=FALSE); return(NULL) }
    bg_df <- as.data.frame(bg_points_terra[, c("x", "y")]); return(bg_df)
  }, error = function(e) { warning("Error generating background points: ", e$message, call. = FALSE); return(NULL) })
}

#' Run SDM Tuning using sdmtune::gridSearch
#'
#' Executes sdmtune grid search with specified settings.
#'
#' @param occs_thinned_sf sf object of thinned occurrence points.
#' @param predictor_stack SpatRaster object of predictors.
#' @param background_df Data frame with 'x', 'y' columns for background points.
#' @param config Configuration list (for SDM settings).
#'
#' @return An SDMtune object containing tuning results, or NULL on error.
#' @export
run_sdm_sdmtune_grid <- function(occs_thinned_sf, predictor_stack, background_df, config) {
  cat("  Running sdmtune::gridSearch...\n")
  
  if (is.null(occs_thinned_sf) || nrow(occs_thinned_sf) == 0 || is.null(predictor_stack) || is.null(background_df)) {
    warning("Invalid inputs for running sdmtune.", call. = FALSE); return(NULL)
  }
  
  occs_coords <- sf::st_coordinates(occs_thinned_sf)
  
  # Prepare SWD object (Species With Data)
  swd_data <- tryCatch({
    sdmtune::prepareSWD(
      species = "species", # sdmtune expects a species name column
      p = occs_coords,
      a = background_df,
      env = predictor_stack
    )
  }, error = function(e) {
    warning("Failed to prepare SWD object: ", e$message, call. = FALSE); return(NULL)
  })
  
  if (is.null(swd_data)) return(NULL)
  
  # Define the model type
  model <- sdmtune::Maxnet() # Using Maxnet as specified in config
  
  # Define the hyperparameter grid (using updated names from config)
  hyper_grid <- config$sdm_tune_grid
  
  # Run grid search with cross-validation
  tuned_model <- NULL
  tryCatch({
    # sdmtune uses 'kfold' for random k-fold
    tuned_model <- sdmtune::gridSearch(
      swd_data,
      hypers = hyper_grid,
      metric = config$sdm_evaluation_metric, # e.g., "AUC"
      method = config$sdm_partitions, # e.g., "kfold"
      k = config$sdm_n_folds, # Number of folds
      # parallel = config$use_parallel, # sdmtune handles parallel differently, often internal
      # numCores = config$num_cores # Not a direct arg here, check sdmtune docs if needed
      save_models = TRUE # Keep model objects for prediction
    )
    cat("  sdmtune::gridSearch completed.\n")
    # print(sdmtune::results(tuned_model)) # Print tuning results
  }, error = function(e) {
    warning("sdmtune::gridSearch failed: ", e$message, call. = FALSE)
    tuned_model <- NULL
  })
  
  return(tuned_model)
}


#' Select Best Model and Predict SDM using sdmtune
#'
#' Selects the best model from sdmtune results and generates a prediction raster.
#'
#' @param sdmtune_results The output object from sdmtune::gridSearch or train.
#' @param predictor_stack SpatRaster object used for prediction.
#' @param config Configuration list.
#'
#' @return A SpatRaster object with the prediction, or NULL on error.
#' @export
predict_sdm_sdmtune <- function(sdmtune_results, predictor_stack, config) {
  cat("  Selecting best model and predicting using sdmtune...\n")
  
  if (is.null(sdmtune_results) || !inherits(sdmtune_results, "SDMtune")) {
    warning("Invalid sdmtune results object provided.", call. = FALSE); return(NULL)
  }
  if (is.null(predictor_stack)) {
    warning("Predictor stack is required for prediction.", call. = FALSE); return(NULL)
  }
  
  # sdmtune results already contain the best model info based on the metric used
  # Access the best model directly if gridSearch stored it
  best_model_obj <- tryCatch({
    # sdmtune::results(sdmtune_results) # View results table
    # sdmtune_results@models # Access list of models if saved
    # The gridSearch object itself often represents the best configuration found
    # OR use optimizeModel if you ran that instead of gridSearch
    # Assuming gridSearch object holds the best config implicitly or use the first model if multiple are equal best
    if(length(sdmtune_results@models) > 0) {
      # Find the row with the best metric score
      res_df <- sdmtune::results(sdmtune_results)
      metric <- config$sdm_evaluation_metric
      # sdmtune uses metric names like AUC_cv, TSS_cv etc. Need to check results() output
      metric_cv_name <- paste0(metric, "_cv") # Likely name format
      if (!metric_cv_name %in% colnames(res_df)) {
        # Fallback or alternative metrics if primary not found
        if("AUC_cv" %in% colnames(res_df)) metric_cv_name <- "AUC_cv"
        else if ("TSS_cv" %in% colnames(res_df)) metric_cv_name <- "TSS_cv"
        else {
          warning("Cannot find metric '", metric_cv_name, "' or fallbacks in sdmtune results. Using first model.", call.=FALSE)
          return(sdmtune_results@models[[1]])
        }
        warning("Metric '", config$sdm_evaluation_metric, "_cv' not found, using '", metric_cv_name, "' instead.", call.=FALSE)
      }
      
      best_row_index <- which.max(res_df[[metric_cv_name]])
      if(length(best_row_index) == 0 || best_row_index > length(sdmtune_results@models)) {
        warning("Could not determine best model index. Using first model.", call.=FALSE)
        best_row_index <- 1
      }
      cat("    Best model hyperparameters (from gridSearch):", paste(names(sdmtune_results@hypers), sdmtune_results@hypers[best_row_index,], collapse=", "), "\n")
      sdmtune_results@models[[best_row_index]]
    } else {
      warning("No models found within the sdmtune results object.", call.=FALSE)
      NULL
    }
  }, error = function(e) {
    warning("Error accessing best model from sdmtune results: ", e$message, call. = FALSE); NULL
  })
  
  if (is.null(best_model_obj)) {
    warning("Could not retrieve the best model object.", call.=FALSE); return(NULL)
  }
  
  # Predict using sdmtune::predict
  prediction_raster <- NULL
  tryCatch({
    prediction_raster <- sdmtune::predict(
      object = best_model_obj, # The model object itself
      data = predictor_stack,
      type = "cloglog" # Maxnet default often cloglog, check model output if needed
    )
    # Rename the layer for consistency
    names(prediction_raster) <- "suitability"
    cat("    Prediction raster generated.\n")
  }, error = function(e) {
    warning("sdmtune::predict failed: ", e$message, call. = FALSE)
    prediction_raster <- NULL
  })
  
  return(prediction_raster)
}

#-------------------------------------------------------------------------------