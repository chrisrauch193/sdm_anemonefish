# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow)
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stats) # Ensure all needed packages are loaded

#' Load, Clean, and Get Coordinates for Individual Species Occurrences
#' (Simplified version focusing on coordinates needed for SWD)
#'
#' @param species_aphia_id AphiaID of the target species.
#' @param occurrence_dir Directory containing individual species CSV files.
#' @param config Project configuration list.
#' @return A matrix of cleaned coordinates (x, y) or NULL on error.
load_clean_individual_occ_coords <- function(species_aphia_id, occurrence_dir, config) {
  # cat("    Loading occurrences for AphiaID:", species_aphia_id, "\n")
  occ_file <- file.path(occurrence_dir, paste0(species_aphia_id, ".csv"))
  if (!file.exists(occ_file)) {
    warning("Occurrence file not found for AphiaID:", species_aphia_id, call. = FALSE)
    return(NULL)
  }
  
  tryCatch({
    occ_df <- readr::read_csv(occ_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
    
    # Basic cleaning for coordinates
    if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occ_df))) {
      warning("Missing coordinate columns in file: ", basename(occ_file), call. = FALSE)
      return(NULL)
    }
    
    occ_clean <- occ_df %>%
      dplyr::mutate(
        decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
        decimalLatitude = suppressWarnings(as.numeric(decimalLatitude))
      ) %>%
      dplyr::filter(
        !is.na(decimalLongitude), !is.na(decimalLatitude),
        decimalLongitude >= -180, decimalLongitude <= 180,
        decimalLatitude >= -90, decimalLatitude <= 90
      )
    
    if (nrow(occ_clean) == 0) {
      warning("No valid coordinates after cleaning for AphiaID:", species_aphia_id, call. = FALSE)
      return(NULL)
    }
    
    # Return just the coordinate matrix
    return(as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")]))
    
  }, error = function(e) {
    warning("Error loading/cleaning occurrences for AphiaID ", species_aphia_id, ": ", e$message, call. = FALSE)
    return(NULL)
  })
}


#' Generate Background Points within the Raster Extent
#'
#' @param predictor_stack SpatRaster object defining the study area.
#' @param n_background Number of background points to generate.
#' @param config Configuration list.
#' @param seed Optional seed for reproducibility.
#' @return A data frame with 'x', 'y' columns for background points, or NULL on error.
generate_sdm_background <- function(predictor_stack, n_background, config, seed = NULL) {
  cat("    Generating", n_background, "background points...\n")
  if (is.null(predictor_stack)) {
    warning("Predictor stack is required to generate background points.", call. = FALSE); return(NULL)
  }
  if (!is.null(seed)) set.seed(seed)
  
  tryCatch({
    bg_points <- terra::spatSample(predictor_stack[[1]],
                                   size = n_background,
                                   method = "random",
                                   na.rm = TRUE, # Essential
                                   xy = TRUE,
                                   warn = FALSE)
    
    if (nrow(bg_points) < n_background) {
      warning("Could only sample ", nrow(bg_points), " background points (fewer than requested ", n_background, "). Check raster extent and NA values.", call. = FALSE)
    }
    if (nrow(bg_points) == 0) {
      warning("Failed to generate any background points.", call. = FALSE); return(NULL)
    }
    
    return(as.data.frame(bg_points[, c("x", "y")]))
    
  }, error = function(e) {
    warning("Error generating background points: ", e$message, call. = FALSE)
    return(NULL)
  })
}

#' Tune SDM Hyperparameters using SDMtune gridSearch with k-fold CV
#'
#' Prepares data, creates folds, and runs SDMtune::gridSearch to find the
#' best hyperparameters based on the mean cross-validation performance.
#'
#' @param occs_coords Matrix or data frame of occurrence coordinates (x, y).
#' @param predictor_stack SpatRaster object of predictors for the *tuning* scenario (usually 'current').
#' @param background_df Data frame with 'x', 'y' columns for background points.
#' @param config Configuration list (for SDM method, tuning grid, CV folds, metric).
#' @param species_name Character string for the species being modeled.
#' @return A list containing the best hyperparameters found (`best_hypers`) and the full SDMtune results object (`tuning_results`), or NULL on error.
run_sdm_tuning_kfold <- function(occs_coords, predictor_stack, background_df, config, species_name = "species") {
  cat("    Preparing data and tuning hyperparameters using", config$sdm_n_folds, "-fold CV...\n")
  
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df)) {
    warning("Invalid inputs for running SDM tuning.", call. = FALSE); return(NULL)
  }
  
  # 1. Prepare SWD object
  full_swd_data <- tryCatch({
    SDMtune::prepareSWD(
      species = species_name,
      p = occs_coords,
      a = background_df,
      env = predictor_stack,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("Failed to prepare SWD object for tuning: ", e$message, call. = FALSE); return(NULL)
  })
  if (is.null(full_swd_data)) return(NULL)
  
  # 2. Create k-folds
  folds <- tryCatch({
    SDMtune::randomFolds(
      data = full_swd_data,
      k = config$sdm_n_folds,
      only_presence = FALSE,
      seed = 123
    )
  }, error = function(e){
    warning("Failed to create k-folds: ", e$message, call.=FALSE); return(NULL)
  })
  if(is.null(folds)) return(NULL)
  
  # 3. Train initial CV base model
  cat("    Training initial CV base model for gridSearch input...\n")
  initial_cv_model <- tryCatch({
    SDMtune::train(method = config$sdm_method,
                   data = full_swd_data,
                   folds = folds,
                   progress = FALSE)
  }, error = function(e) {
    warning("Failed to train initial CV base model: ", e$message, call. = FALSE); return(NULL)
  })
  if (is.null(initial_cv_model) || !inherits(initial_cv_model, "SDMmodelCV")) {
    warning("Failed to create a valid SDMmodelCV object.", call. = FALSE); return(NULL)
  }
  
  # 4. Define hyperparameter grid
  hyper_grid <- config$sdm_tune_grid
  
  # 5. Run gridSearch
  tuning_results <- NULL
  tryCatch({
    cat("    Running SDMtune::gridSearch with k-fold CV...\n")
    tuning_results <- SDMtune::gridSearch(
      model = initial_cv_model,
      hypers = hyper_grid,
      metric = config$sdm_evaluation_metric, # Metric for CV evaluation
      save_models = TRUE,
      progress = FALSE,
      interactive = FALSE
    )
    cat("    SDMtune::gridSearch with k-fold CV completed.\n")
    
    # Access results using the @results slot
    res_df <- tuning_results@results
    
    if(is.null(res_df) || nrow(res_df) == 0) {
      warning("No results found in the SDMtune object after gridSearch.", call.=FALSE)
      return(NULL)
    }
    
    # *** CORRECTED METRIC COLUMN NAME IDENTIFICATION ***
    # When gridSearch uses an SDMmodelCV object, results columns are named
    # train_METRIC and test_METRIC (representing the mean across folds).
    metric_base_name <- toupper(config$sdm_evaluation_metric) # e.g., "AUC", "TSS"
    # AICc might be an exception, check both cases
    if (tolower(config$sdm_evaluation_metric) == "aicc") {
      metric_cv_name <- "AICc" # AICc column doesn't have train_/test_ prefix
    } else {
      metric_cv_name <- paste0("test_", metric_base_name) # e.g., test_AUC, test_TSS
    }
    
    # Check if the expected metric column exists
    if (!metric_cv_name %in% colnames(res_df)) {
      warning("Cannot find CV evaluation metric column '", metric_cv_name,
              "' in the results dataframe. Available columns: ",
              paste(colnames(res_df), collapse=", "),
              ". Cannot select best hypers.", call. = FALSE)
      print(head(res_df))
      return(NULL)
    }
    
    # Select best model (handle NAs, minimize AICc, maximize others)
    valid_metric_indices <- which(!is.na(res_df[[metric_cv_name]]))
    if(length(valid_metric_indices) == 0) {
      warning("All values for metric '", metric_cv_name, "' are NA. Cannot select best hypers.", call.=FALSE)
      return(NULL)
    }
    
    if (metric_col_name == "AICc") {
      best_row_relative_index <- which.min(res_df[[metric_cv_name]][valid_metric_indices])
    } else {
      best_row_relative_index <- which.max(res_df[[metric_cv_name]][valid_metric_indices])
    }
    
    best_row_index <- valid_metric_indices[best_row_relative_index]
    
    if(length(best_row_index) == 0 || best_row_index > nrow(res_df)) {
      warning("Could not determine best model index from CV results. Using first row as fallback.", call.=FALSE)
      best_row_index <- 1 # Fallback
    }
    
    # Extract the best hyperparameter combination from the results df
    best_hypers_df <- res_df[best_row_index, names(hyper_grid), drop = FALSE]
    
    cat("    Best hyperparameters found (Mean CV", metric_cv_name, "):",
        paste(names(best_hypers_df), best_hypers_df[1,], collapse=", "), "\n")
    
    return(list(best_hypers = best_hypers_df, tuning_results = tuning_results))
    
  }, error = function(e) {
    warning("SDMtune::gridSearch failed during CV tuning: ", e$message, call. = FALSE)
    return(NULL)
  })
}


#' Train Final SDM Model
#'
#' Trains a single SDM model using all provided occurrence and background data
#' with the specified best hyperparameters.
#'
#' @param occs_coords Matrix or data frame of ALL occurrence coordinates (x, y).
#' @param predictor_stack SpatRaster object of predictors for the training scenario (usually 'current').
#' @param background_df Data frame with ALL 'x', 'y' columns for background points.
#' @param best_hypers A data frame or list containing the best hyperparameter values.
#' @param config Configuration list.
#' @param species_name Character string for the species being modeled.
#' @return An SDMmodel object or NULL on error.
train_final_sdm <- function(occs_coords, predictor_stack, background_df, best_hypers, config, species_name = "species") {
  cat("    Training final model on full dataset with best hyperparameters...\n")
  
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df) || is.null(best_hypers)) {
    warning("Invalid inputs for training final SDM.", call. = FALSE); return(NULL)
  }
  
  # 1. Prepare SWD object with ALL data
  full_swd_data <- tryCatch({
    SDMtune::prepareSWD(
      species = species_name,
      p = occs_coords,
      a = background_df,
      env = predictor_stack,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("Failed to prepare SWD object for final model: ", e$message, call. = FALSE); return(NULL)
  })
  if (is.null(full_swd_data)) return(NULL)
  
  # 2. Prepare hyperparameter arguments for train()
  train_args <- list(
    method = config$sdm_method,
    data = full_swd_data,
    progress = FALSE
  )
  
  # Add hyperparameters dynamically
  for (h_name in names(best_hypers)) {
    h_value <- best_hypers[[h_name]]
    # Ensure correct type if needed (e.g., numeric conversion if stored as char)
    if (h_name == "reg") h_value <- as.numeric(h_value)
    # Feature classes should remain character
    train_args[[h_name]] <- h_value
  }
  cat("      Using hyperparameters:", paste(names(train_args)[-(1:2)], sapply(train_args[-(1:2)], as.character), collapse=", "), "\n")
  
  
  # 3. Train the final model using do.call
  final_model <- tryCatch({
    do.call(SDMtune::train, train_args)
  }, error = function(e) {
    warning("Failed to train final model: ", e$message, call. = FALSE); return(NULL)
  })
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) {
    warning("Final model training did not return a valid SDMmodel object.", call.=FALSE)
    return(NULL)
  }
  
  cat("    Final model trained successfully.\n")
  return(final_model)
}


#' Predict SDM Suitability
#'
#' Generates a prediction raster using a trained final SDMmodel object and a
#' predictor stack for the target scenario.
#'
#' @param final_sdm_model An SDMmodel object (output from `train_final_sdm`).
#' @param predictor_stack SpatRaster object for the prediction scenario.
#' @param config Configuration list.
#' @param output_type Character string for prediction type (e.g., "cloglog", "logistic"). Defaults to "cloglog".
#' @return A SpatRaster object with the prediction, or NULL on error.
predict_sdm_suitability <- function(final_sdm_model, predictor_stack, config, output_type = "cloglog") {
  cat("    Predicting suitability for scenario...\n")
  
  if (is.null(final_sdm_model) || !inherits(final_sdm_model, "SDMmodel")) {
    warning("Invalid final_sdm_model object provided.", call. = FALSE); return(NULL)
  }
  if (is.null(predictor_stack)) {
    warning("Predictor stack is required for prediction.", call. = FALSE); return(NULL)
  }
  
  prediction_raster <- NULL
  tryCatch({
    prediction_raster <- SDMtune::predict(
      object = final_sdm_model,
      data = predictor_stack,
      type = output_type,
      clamp = TRUE # Consistent with paper's likely Maxent approach
    )
    names(prediction_raster) <- "suitability"
    cat("    Prediction raster generated (type:", output_type, ").\n")
  }, error = function(e) {
    warning("SDMtune::predict failed: ", e$message, call. = FALSE)
    prediction_raster <- NULL
  })
  
  return(prediction_raster)
}

#-------------------------------------------------------------------------------