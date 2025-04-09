# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow v3 - Ensemble)
# - Trains multiple algorithms (Maxnet, RF, BRT).
# - Predicts using an ensemble mean.
# - Calculates VI based on the Maxnet model within the ensemble list.
# - Saves/Loads lists of models.
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr, future, furrr) # Added future/furrr

# --- (load_clean_individual_occ_coords remains the same) ---
#' Load, Clean, and Get Coordinates for Individual Species Occurrences
#' (Returns list with data and count)
#'
#' @param species_aphia_id AphiaID of the target species.
#' @param occurrence_dir Directory containing individual species CSV files.
#' @param config Project configuration list (needs `occurrence_crs`, `thinning_method`, `min_occurrences_sdm`, and optionally `predictor_stack_for_thinning` if thinning is 'cell').
#' @param logger A log4r logger object (can be NULL for less verbose output).
#' @param species_log_file Optional path to a species-specific log file for detailed messages.
#' @return A list containing cleaned coordinates `coords` (matrix) and the final `count`, or NULL on error.
load_clean_individual_occ_coords <- function(species_aphia_id, occurrence_dir, config, logger, species_log_file = NULL) {
  
  # Local logging function for this helper
  hlog <- function(level, ...) {
    msg <- paste(Sys.time(), paste0("[",level,"]"), "[OccLoadHelper]", paste0(..., collapse = " "))
    if (!is.null(species_log_file)) { # Prefer species log if provided
      cat(msg, "\n", file = species_log_file, append = TRUE)
    } else if (!is.null(logger)) { # Fallback to main logger
      log_level_func <- switch(level, DEBUG=log4r::debug, INFO=log4r::info, WARN=log4r::warn, log4r::error)
      tryCatch(log_level_func(logger, msg), error = function(e) cat(msg, "\n")) # Log to main logger, fallback to cat
    } else {
      # cat(msg, "\n") # Fallback to console if no loggers provided
    }
  }
  
  hlog("DEBUG", paste("Loading occurrences for AphiaID:", species_aphia_id))
  occ_file <- file.path(occurrence_dir, paste0(species_aphia_id, ".csv"))
  if (!file.exists(occ_file)) {
    hlog("WARN", paste("Occurrence file not found:", basename(occ_file)))
    return(list(coords = NULL, count = 0))
  }
  
  tryCatch({
    occ_df <- readr::read_csv(occ_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
    
    if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occ_df))) {
      hlog("WARN", paste("Missing coordinate columns in file:", basename(occ_file)))
      return(list(coords = NULL, count = 0))
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
    
    count_after_clean <- nrow(occ_clean)
    hlog("DEBUG", paste("  Retained", count_after_clean, "records after coordinate cleaning."))
    
    if (count_after_clean == 0) {
      hlog("WARN", "No valid coordinates after cleaning.")
      return(list(coords = NULL, count = 0))
    }
    
    # --- Spatial Thinning (Optional) ---
    if (!is.null(config$thinning_method) && config$thinning_method == "cell" && !is.null(config$predictor_stack_for_thinning)) {
      hlog("DEBUG", "  Applying cell-based thinning...")
      predictor_stack_thin <- config$predictor_stack_for_thinning
      if (is.null(predictor_stack_thin) || !inherits(predictor_stack_thin, "SpatRaster") || terra::nlyr(predictor_stack_thin) == 0) {
        hlog("WARN", " predictor_stack_for_thinning invalid/missing. Skipping thinning.")
        occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      } else {
        occs_sf <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs)
        target_crs <- terra::crs(predictor_stack_thin)
        if (!identical(sf::st_crs(occs_sf), sf::st_crs(target_crs))) {
          hlog("DEBUG", "Transforming occurrence CRS for thinning.")
          occs_sf <- sf::st_transform(occs_sf, crs = target_crs)
        }
        
        occs_spatvector <- terra::vect(occs_sf)
        occs_cells <- tryCatch(terra::extract(predictor_stack_thin[[1]], occs_spatvector, cells = TRUE), error=function(e){hlog("WARN",paste("Error extracting cells:",e$message));NULL})
        
        if(is.null(occs_cells)){
          hlog("WARN", " Cell extraction failed. Skipping thinning.")
          occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
        } else {
          valid_cells_indices <- which(!is.na(occs_cells$cell))
          if(length(valid_cells_indices) > 0) {
            # Ensure we keep the row associated with the *first* occurrence in a duplicated cell
            first_occurrence_in_cell_idx <- !duplicated(occs_cells$cell[valid_cells_indices])
            original_indices_to_keep <- valid_cells_indices[first_occurrence_in_cell_idx]
            
            occs_thinned_sf <- occs_sf[original_indices_to_keep, ]
            occ_thinned_coords <- sf::st_coordinates(occs_thinned_sf)
            colnames(occ_thinned_coords) <- c("decimalLongitude", "decimalLatitude") # Ensure names
            hlog("DEBUG", "Thinned from", nrow(occs_sf), "to", nrow(occ_thinned_coords))
          } else {
            hlog("WARN", " No valid cells found for occurrences on thinning raster. Skipping thinning.")
            occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
          }
        }
        rm(occs_sf, occs_spatvector, occs_cells); gc() # Clean up SF/terra objects
      }
      count_after_thin <- nrow(occ_thinned_coords)
      hlog("DEBUG", paste("  Retained", count_after_thin, "records after thinning."))
      if (count_after_thin == 0) {hlog("WARN", "No records left after thinning."); return(list(coords = NULL, count = 0))}
    } else {
      hlog("DEBUG", "  Skipping spatial thinning or method not 'cell'.")
      occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      count_after_thin <- nrow(occ_thinned_coords)
    }
    
    colnames(occ_thinned_coords) <- c("x", "y") # Ensure standard names for SWD
    return(list(coords = occ_thinned_coords, count = count_after_thin))
    
  }, error = function(e) {
    hlog("ERROR", paste("Error processing occurrences for AphiaID", species_aphia_id, ":", e$message))
    return(list(coords = NULL, count = 0))
  })
}


# --- (generate_sdm_background remains the same) ---
#' Generate Background Points within the Raster Extent
#' @param predictor_stack SpatRaster stack (used for extent/masking).
#' @param n_background Number of points to attempt generating.
#' @param config Project configuration list.
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @param seed Optional random seed for reproducibility.
#' @return A data frame of background point coordinates (x, y), or NULL on error.
generate_sdm_background <- function(predictor_stack, n_background, config, logger, species_log_file = NULL, seed = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[BgGenHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  hlog("DEBUG", paste("Generating", n_background, "background points..."))
  if (is.null(predictor_stack) || !inherits(predictor_stack, "SpatRaster") || terra::nlyr(predictor_stack) == 0) { hlog("ERROR", "Predictor stack is required."); return(NULL) }
  if (!is.null(seed)) set.seed(seed)
  
  tryCatch({
    # Use first layer as mask
    mask_layer <- predictor_stack[[1]]
    bg_points <- terra::spatSample(mask_layer, size = n_background, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE, replace=FALSE)
    if (nrow(bg_points) < n_background) { hlog("WARN", paste("Could only sample", nrow(bg_points), "background points (requested", n_background, ").")) }
    if (nrow(bg_points) == 0) { hlog("ERROR", "Failed to generate any background points."); return(NULL) }
    hlog("DEBUG", paste("Generated", nrow(bg_points), "background points."))
    return(as.data.frame(bg_points[, c("x", "y")]))
  }, error = function(e) { hlog("ERROR", paste("Error generating background points:", e$message)); return(NULL) })
}

# --- (run_sdm_tuning_kfold remains the same - tunes only Maxnet) ---
#' Tune SDM Hyperparameters using SDMtune gridSearch with k-fold CV
#' Returns the tuning object containing results and models.
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @return An `SDMtune` object containing tuning results, or NULL on error.
run_sdm_tuning_kfold <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species", species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[TuningHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}# cat(msg, "\n")} }
  
  hlog("INFO", paste("Tuning hyperparameters for", species_name, "using", config$sdm_n_folds, "-fold CV (focusing on", config$sdm_method, ")..."))
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df)) { hlog("ERROR", "Invalid inputs."); return(NULL) }
  
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE) }, error = function(e) { hlog("ERROR", paste("Failed prepareSWD:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  folds <- tryCatch({ SDMtune::randomFolds(data = full_swd_data, k = config$sdm_n_folds, only_presence = FALSE, seed = 123)}, error = function(e){ hlog("ERROR", paste("Failed create k-folds:", e$message)); return(NULL)})
  if(is.null(folds)) return(NULL)
  
  hlog("DEBUG", "Training initial CV base model (", config$sdm_method, ") for gridSearch...")
  initial_cv_model <- tryCatch({ SDMtune::train(method = config$sdm_method, data = full_swd_data, folds = folds, progress = FALSE)}, error = function(e) { hlog("ERROR", paste("Failed train CV base model:", e$message)); return(NULL)})
  if (is.null(initial_cv_model) || !inherits(initial_cv_model, "SDMmodelCV")) { hlog("ERROR", "Failed create valid SDMmodelCV object."); return(NULL)}
  
  hyper_grid <- config$sdm_tune_grid; tuning_results <- NULL
  tryCatch({
    hlog("DEBUG", "Running SDMtune::gridSearch with k-fold CV...")
    show_progress <- FALSE; if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::level(logger); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE})}
    
    tuning_results <- SDMtune::gridSearch(model = initial_cv_model, hypers = hyper_grid, metric = config$sdm_evaluation_metric, save_models = TRUE, progress = show_progress, interactive = FALSE)
    hlog("DEBUG", "SDMtune::gridSearch completed.")
    
    res_df <- tuning_results@results
    if(is.null(res_df) || nrow(res_df) == 0) { hlog("WARN", "No results found in tuning object."); return(NULL) }
    
    metric_base_upper <- toupper(config$sdm_evaluation_metric)
    target_metric_col <- if (tolower(config$sdm_evaluation_metric) == "aicc") "AICc" else paste0("test_", metric_base_upper)
    if (!target_metric_col %in% colnames(res_df)) { hlog("ERROR", paste("CV metric '", target_metric_col, "' not found. Available:", paste(colnames(res_df), collapse=", "))); print(head(res_df)); return(NULL) }
    
    valid_metric_indices <- which(!is.na(res_df[[target_metric_col]]))
    if(length(valid_metric_indices) == 0) { hlog("WARN", paste("All values for metric '", target_metric_col, "' are NA.")); return(NULL) }
    
    select_fun <- if (target_metric_col == "AICc") which.min else which.max
    best_row_relative_index <- select_fun(res_df[[target_metric_col]][valid_metric_indices])
    best_row_index <- valid_metric_indices[best_row_relative_index]
    if(length(best_row_index) == 0 || best_row_index > nrow(res_df)) { hlog("WARN", "Could not determine best model index, using first valid row."); best_row_index <- valid_metric_indices[1]}
    
    best_hypers_df <- res_df[best_row_index, names(hyper_grid), drop = FALSE]
    hlog("INFO", paste("  Best hypers (Mean CV", target_metric_col,"=", round(res_df[[target_metric_col]][best_row_index], 4),"): ", paste(names(best_hypers_df), best_hypers_df[1,], collapse=", ")))
    
    # Add best hypers to the object attributes for easy access
    attr(tuning_results, "best_hypers") <- best_hypers_df
    return(tuning_results) # Return the whole tuning object
    
  }, error = function(e) { hlog("ERROR", paste("SDMtune::gridSearch failed:", e$message)); return(NULL) })
}


# --- MODIFIED: train_final_sdm for Ensemble ---
#' Train Final SDM Models (Ensemble: Maxnet, RF, BRT)
#' Uses the same best hyperparameters derived from tuning the main method (Maxnet).
#' Returns a LIST of trained model objects.
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @return A named list of `SDMmodel` objects (e.g., list(Maxnet=..., RF=..., BRT=...)), or NULL on error.
train_final_sdm <- function(occs_coords, predictor_stack, background_df, best_hypers, config, logger, species_name = "species", species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[TrainHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  hlog("INFO", paste("Training final ENSEMBLE models for", species_name, "on full dataset..."))
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df) || is.null(best_hypers)) { hlog("ERROR", "Invalid inputs."); return(NULL) }
  
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE)}, error = function(e) { hlog("ERROR", paste("Failed prepareSWD:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  # Algorithms to include in the ensemble
  ensemble_methods <- c("Maxnet", "RF", "BRT") # Exclude ANN for now, may need absence data
  hlog("INFO", paste("  Ensemble algorithms:", paste(ensemble_methods, collapse=", ")))
  
  # Prepare arguments, applying best_hypers derived from Maxnet tuning to all
  base_args <- list(data = full_swd_data, progress = FALSE)
  for (h_name in names(best_hypers)) {
    h_value <- best_hypers[[h_name]]
    # Ensure correct type if needed (e.g., reg for Maxnet)
    if (h_name == "reg" && inherits(best_hypers[[h_name]], "factor")) {
      h_value <- as.numeric(as.character(h_value)) # Convert factor to numeric if needed
    } else if (h_name == "reg") {
      h_value <- as.numeric(h_value) # Ensure numeric
    }
    # Apply hyperparameter to base_args IF it's relevant for SDMtune::train, otherwise ignore
    # (Maxnet uses reg/fc, RF uses mtry/ntree/nodesize, BRT uses n.trees etc.)
    # For simplicity, we apply the Maxnet hypers (reg/fc) and use defaults for others.
    # More sophisticated: map relevant hypers if tuning involved them.
    if(h_name %in% c("reg", "fc")) { # Apply only Maxnet-tuned hypers for now
      base_args[[h_name]] <- h_value
    }
  }
  hlog("DEBUG", paste("  Using base hyperparameters for all models:", paste(names(best_hypers), sapply(best_hypers, as.character), collapse=", ")))
  
  # Use future_map to train in parallel if configured
  # Ensure plan is set before calling the main processing function
  show_progress <- FALSE # Progress bar inside future_map is tricky, disable for now
  
  trained_models_list <- tryCatch({
    furrr::future_map(ensemble_methods, function(method_name) {
      hlog("DEBUG", paste("    Training final model for method:", method_name))
      method_args <- base_args
      method_args$method <- method_name
      # Add default args if needed for RF/BRT (optional, train uses defaults)
      # if (method_name == "RF") { method_args$ntree <- 500 } # Example
      # if (method_name == "BRT") { method_args$n.trees <- 1000 } # Example
      
      # Train the model
      model_obj <- tryCatch({
        do.call(SDMtune::train, method_args)
      }, error = function(e) {
        hlog("ERROR", paste("    Failed train final model for method", method_name, ":", e$message)); NULL
      })
      
      # Basic validation
      if(is.null(model_obj) || !inherits(model_obj, "SDMmodel")) {
        hlog("ERROR", paste("    Training returned invalid object for method", method_name))
        return(NULL) # Return NULL on failure for this method
      } else {
        hlog("DEBUG", paste("    Successfully trained model for method:", method_name))
        return(model_obj) # Return the trained model object
      }
    }, .options = furrr_options(seed = TRUE)) # Ensure reproducibility if parallel
  }, error = function(e) {
    hlog("ERROR", paste("Error during parallel model training:", e$message))
    return(NULL)
  })
  
  # Check results and name the list
  if (is.null(trained_models_list) || length(trained_models_list) == 0 || all(sapply(trained_models_list, is.null))) {
    hlog("ERROR", "Failed to train ANY models for the ensemble.")
    return(NULL)
  }
  
  # Filter out NULLs (failed models) and name the list
  valid_models <- trained_models_list[!sapply(trained_models_list, is.null)]
  names(valid_models) <- ensemble_methods[!sapply(trained_models_list, is.null)]
  
  if(length(valid_models) == 0) {
    hlog("ERROR", "No models trained successfully.")
    return(NULL)
  }
  
  hlog("INFO", paste("  Final ensemble models trained successfully for:", species_name, "(Methods:", paste(names(valid_models), collapse=", "), ")"))
  return(valid_models) # Return the named list of models
}


# --- MODIFIED: predict_sdm_suitability for Ensemble ---
#' Predict SDM Suitability (Ensemble Mean)
#' Takes a list of models, predicts with each, returns the mean raster.
#' @param final_model_list A named list of trained `SDMmodel` objects.
#' @return A SpatRaster object with the mean prediction, or error message string.
predict_sdm_suitability <- function(final_model_list, predictor_stack, config, logger, species_log_file = NULL, output_type = "cloglog") {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[PredictHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  hlog("DEBUG", paste("Attempting ENSEMBLE prediction (type:", output_type, ")..."))
  if (is.null(final_model_list) || !is.list(final_model_list) || length(final_model_list) == 0 || !all(sapply(final_model_list, inherits, "SDMmodel"))) { msg <- "Invalid final_model_list provided."; hlog("ERROR", msg); return(msg) }
  if (is.null(predictor_stack)) { msg <- "Predictor stack required."; hlog("ERROR", msg); return(msg) }
  
  prediction_rasters <- list()
  errors_occurred <- FALSE
  
  for (model_name in names(final_model_list)) {
    model <- final_model_list[[model_name]]
    hlog("DEBUG", paste("  Predicting with model:", model_name))
    pred <- tryCatch({
      SDMtune::predict(object = model, data = predictor_stack, type = output_type, clamp = TRUE)
    }, error = function(e) {
      err_msg <- paste("  Prediction failed for model", model_name, ":", e$message); hlog("WARN", err_msg); errors_occurred <<- TRUE; NULL
    })
    if (!is.null(pred) && inherits(pred, "SpatRaster")) {
      prediction_rasters[[model_name]] <- pred
    } else {
      hlog("WARN", paste(" Prediction result invalid for", model_name))
      errors_occurred <<- TRUE
    }
  }
  
  if (length(prediction_rasters) == 0) {
    msg <- "No successful predictions generated for any model in the ensemble."; hlog("ERROR", msg); return(msg)
  }
  
  # Stack the successful predictions
  pred_stack <- tryCatch(terra::rast(prediction_rasters), error = function(e){
    msg <- paste("Failed to stack individual predictions:", e$message); hlog("ERROR", msg); return(msg)
  })
  if(inherits(pred_stack, "character")) return(pred_stack) # Return error message if stacking failed
  
  # Calculate the mean ensemble prediction
  ensemble_mean <- tryCatch({
    terra::mean(pred_stack, na.rm = TRUE)
  }, error = function(e) {
    msg <- paste("Failed to calculate ensemble mean:", e$message); hlog("ERROR", msg); return(msg)
  })
  if(inherits(ensemble_mean, "character")) return(ensemble_mean) # Return error message
  
  names(ensemble_mean) <- "suitability_ensemble_mean"
  hlog("INFO", paste("Ensemble mean prediction calculated using", length(prediction_rasters), "models."))
  if(errors_occurred) hlog("WARN", "Some models failed prediction during ensemble calculation.")
  
  # Optionally add SD raster here if needed
  # ensemble_sd <- tryCatch(terra::app(pred_stack, fun = sd, na.rm = TRUE), error=function(e) NULL)
  # if(!is.null(ensemble_sd)) names(ensemble_sd) <- "suitability_ensemble_sd"
  # return(c(ensemble_mean, ensemble_sd)) # If returning both mean and SD
  
  return(ensemble_mean)
}

# --- MODIFIED: save_tuning_results (No change needed, already saves based on tuning run) ---
# --- (Keep the existing save_tuning_results function) ---
#' Save Tuning Results (RDS and CSV)
#' Saves the full tuning object and the results table separately.
#' Constructs paths based on target structure in config.
#'
#' @param tuning_output The object returned by `run_sdm_tuning_kfold` (contains `@results`).
#' @param species_name_sanitized Sanitized species name (e.g., "Genus_species").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_combined_pca").
#' @param config The configuration list (needs target_results_base, model_output_subdir_map).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_tuning_results <- function(tuning_output, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveTuneHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  if (is.null(tuning_output) || !inherits(tuning_output, "SDMtune")) { hlog("ERROR", "Invalid tuning object provided."); return(FALSE) }
  
  # Determine target subdirectory using intermediate path
  target_subdir <- file.path(config$results_dir_intermediate, paste0(basename(config$occurrence_dir), predictor_type_suffix))
  dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
  
  target_base_name <- paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix)
  
  # 1. Save full tuning object (RDS)
  rds_file <- file.path(target_subdir, paste0(target_base_name, "_object.rds"))
  tryCatch({ saveRDS(tuning_output, file = rds_file); hlog("DEBUG", "Full tuning object saved (RDS):", basename(rds_file)) },
           error = function(e) { hlog("ERROR", paste("Failed save tuning RDS:", e$message)) })
  
  # 2. Save results table (CSV)
  csv_file <- file.path(target_subdir, paste0(target_base_name, "_results.csv"))
  results_df <- tuning_output@results
  if (is.null(results_df) || nrow(results_df) == 0) { hlog("WARN", "No results table found in tuning object."); return(TRUE) }
  
  tryCatch({ readr::write_csv(results_df, file = csv_file); hlog("DEBUG", "Tuning results table saved (CSV):", basename(csv_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save tuning CSV:", e$message)); FALSE })
}

# --- MODIFIED: save_final_model for Ensemble ---
#' Save Final SDM Model Objects (List)
#' Saves the LIST of trained SDMmodel objects as a single RDS file.
#' NOTE: Saves to *intermediate* location.
#'
#' @param final_model_list The named list of trained `SDMmodel` objects.
#' @return TRUE on success, FALSE on failure.
save_final_model <- function(final_model_list, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveModelHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  if (is.null(final_model_list) || !is.list(final_model_list) || length(final_model_list) == 0) { hlog("ERROR", "Invalid model list provided."); return(FALSE) }
  
  model_subdir <- file.path(config$models_dir_intermediate, paste0(basename(config$occurrence_dir), predictor_type_suffix))
  dir.create(model_subdir, recursive = TRUE, showWarnings = FALSE)
  model_file <- file.path(model_subdir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  tryCatch({ saveRDS(final_model_list, file = model_file); hlog("DEBUG", "Final model LIST saved (RDS):", basename(model_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save final model list RDS:", e$message)); FALSE })
}

# --- (construct_prediction_filename remains the same - handles suffix) ---
# --- (save_sdm_prediction remains the same - saves the raster passed to it) ---
#' Construct Target Prediction Filename
#' Builds the expected filename and path for a prediction raster based on target structure.
#' Used for checking if a prediction file already exists.
#'
#' @param species_name_sanitized Sanitized species name (e.g., "Genus_species").
#' @param scenario_name The name of the scenario (e.g., "current", "ssp119_2050").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_combined_pca").
#' @param config The configuration list.
#' @return Character string with the full path to the expected prediction file.
construct_prediction_filename <- function(species_name_sanitized, scenario_name, predictor_type_suffix, config) {
  
  # Construct base filename part (e.g., mean_pred_SPECIES)
  # Using sanitized name for now.
  base_filename <- paste0("mean_pred_", species_name_sanitized) # Using SDMtune mean convention implicitly
  
  target_dir <- NULL
  target_filename_stem <- NULL # Filename without extension
  
  if (scenario_name == "current") {
    target_dir <- config$target_predictions_current_dir
    # Filename format: mean_pred_SPECIES[_MODELTYPE?]
    target_filename_stem <- paste0(base_filename, predictor_type_suffix)
  } else {
    # Future scenario
    ssp_dir_name <- config$ssp_scenario_map[[scenario_name]]
    if (is.null(ssp_dir_name)) {
      warning("No SSP directory mapping for scenario: ", scenario_name, call. = FALSE)
      return(NULL) # Return NULL if mapping is missing
    }
    target_dir <- file.path(config$target_predictions_future_base, ssp_dir_name)
    
    # Extract time tag for filename
    time_tag_clean <- ifelse(grepl("2050$", scenario_name), "dec50", ifelse(grepl("2100$", scenario_name), "dec100", "UNKNOWN"))
    if(time_tag_clean == "UNKNOWN") {
      warning("Cannot extract time tag from scenario: ", scenario_name, call. = FALSE)
      return(NULL) # Return NULL if time tag is unknown
    }
    
    # Filename format: mean_pred_SPECIES_SSP_TIME[_MODELTYPE?]
    target_filename_stem <- paste0(base_filename, "_", ssp_dir_name, "_", time_tag_clean, predictor_type_suffix)
  }
  
  if(is.null(target_dir) || is.null(target_filename_stem)) {
    warning("Could not construct target path components.", call. = FALSE)
    return(NULL)
  }
  
  # Add extension
  target_filename <- paste0(target_filename_stem, ".tif")
  
  return(file.path(target_dir, target_filename))
}

#' Save SDM Prediction Raster
#' Saves the prediction SpatRaster to the target output directory structure.
#' Handles current vs future scenario paths.
#'
#' @param prediction_raster The predicted SpatRaster object.
#' @param species_name_sanitized Sanitized species name.
#' @param scenario_name The name of the scenario (e.g., "current", "ssp119_2050").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_combined_pca").
#' @param config The configuration list (needs target prediction dirs, ssp_scenario_map).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_sdm_prediction <- function(prediction_raster, species_name_sanitized, scenario_name, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SavePredHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  if (is.null(prediction_raster) || !inherits(prediction_raster, "SpatRaster")) { hlog("ERROR", "Invalid prediction raster provided."); return(FALSE) }
  
  pred_file_path <- construct_prediction_filename(species_name_sanitized, scenario_name, predictor_type_suffix, config)
  if(is.null(pred_file_path)) { hlog("ERROR", "Failed construct prediction path."); return(FALSE) }
  
  target_dir <- dirname(pred_file_path)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    terra::writeRaster(prediction_raster, filename = pred_file_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    hlog("DEBUG", paste("Prediction raster saved:", basename(pred_file_path), "to", target_dir))
    TRUE
  }, error = function(e) {
    hlog("ERROR", paste("Failed save prediction raster:", e$message))
    FALSE
  })
}

# --- MODIFIED: calculate_and_save_vi for Ensemble ---
#' Calculate and Save Variable Importance (v9 - Ensemble Adjusted)
#' Calculates permutation importance using SDMtune::varImp for the **Maxnet model**
#' within the provided list. Requires the model list and training SWD data.
#'
#' @param final_model_list A named list of `SDMmodel` objects.
#' @param training_swd An `SWD` object with full training data. REQUIRED.
#' @return TRUE on success, FALSE on failure.
calculate_and_save_vi <- function(final_model_list, training_swd, # Both are now required
                                  species_name_sanitized, group_name, predictor_type_suffix,
                                  config, logger, species_log_file = NULL) {
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[VarImpHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  # --- Input validation ---
  if (is.null(final_model_list) || !is.list(final_model_list) || length(final_model_list) == 0 || is.null(names(final_model_list))) {
    hlog("ERROR", "Invalid named list of SDMmodel objects provided for VI.")
    return(FALSE)
  }
  # Check specifically for Maxnet model
  maxnet_model <- final_model_list[["Maxnet"]]
  if (is.null(maxnet_model) || !inherits(maxnet_model, "SDMmodel")) {
    hlog("ERROR", "Maxnet model not found or invalid within the provided list for VI.")
    return(FALSE)
  }
  if (is.null(training_swd) || !inherits(training_swd, "SWD")) {
    hlog("ERROR", "Training SWD object required for varImp calculation.")
    return(FALSE)
  }
  # --- End Input Validation ---
  
  hlog("INFO", "Calculating variable importance (permutation for Maxnet model)...")
  vi_results <- NULL
  tryCatch({
    
    # Determine if progress should be shown
    show_progress <- FALSE
    if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::level(logger); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE}) }
    
    # --- Call varImp on the specific Maxnet model object ---
    hlog("DEBUG", "Calling varImp on the Maxnet model object...")
    
    vi_results <- SDMtune::varImp(
      model = maxnet_model, # Use the extracted Maxnet model
      permut = 10,          # Number of permutations
      progress = show_progress
    )
    # --- End varImp call ---
    
    if (is.null(vi_results) || nrow(vi_results) == 0) { hlog("WARN", "Variable importance calculation returned empty results."); return(TRUE)} # Still return TRUE if VI is empty but didn't error
    
    # --- Saving Logic (uses the overall predictor_type_suffix for the filename) ---
    vi_subdir_name <- paste0("vi_", group_name)
    target_subdir <- file.path(config$target_vi_base, vi_subdir_name)
    dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
    vi_filename <- paste0("vi_", species_name_sanitized, predictor_type_suffix, ".csv") # Uses the main suffix (e.g., _pca, _combined_pca)
    vi_file_path <- file.path(target_subdir, vi_filename)
    readr::write_csv(vi_results, vi_file_path)
    hlog("INFO", paste("Maxnet variable importance saved to:", vi_file_path))
    return(TRUE)
    
  }, error = function(e) {
    hlog("ERROR", paste("Variable importance calculation/saving failed:", e$message))
    return(FALSE)
  })
}

#-------------------------------------------------------------------------