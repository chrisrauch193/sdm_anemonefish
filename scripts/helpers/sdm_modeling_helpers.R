# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow v2)
# - Saving logic moved to dedicated helper functions.
# - Path construction uses new config structure for target output.
# - Added Variable Importance helper.
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr) # Added readr

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
      log_level_func(logger, msg) # Log to main logger
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
        if (sf::st_crs(occs_sf) != sf::st_crs(target_crs)) {
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
            unique_cell_indices_in_valid <- which(!duplicated(occs_cells$cell[valid_cells_indices]))
            original_indices_to_keep <- valid_cells_indices[unique_cell_indices_in_valid]
            occs_thinned_sf <- occs_sf[original_indices_to_keep, ]
            occ_thinned_coords <- sf::st_coordinates(occs_thinned_sf)
            colnames(occ_thinned_coords) <- c("decimalLongitude", "decimalLatitude") # Ensure names
          } else {
            hlog("WARN", " No valid cells found for occurrences on thinning raster. Skipping thinning.")
            occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
          }
        }
      }
      count_after_thin <- nrow(occ_thinned_coords)
      hlog("DEBUG", paste("  Retained", count_after_thin, "records after thinning."))
      if (count_after_thin == 0) {hlog("WARN", "No records left after thinning."); return(list(coords = NULL, count = 0))}
    } else {
      hlog("DEBUG", "  Skipping spatial thinning or method not 'cell'.")
      occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      count_after_thin <- nrow(occ_thinned_coords)
    }
    
    return(list(coords = occ_thinned_coords, count = count_after_thin))
    
  }, error = function(e) {
    hlog("ERROR", paste("Error processing occurrences for AphiaID", species_aphia_id, ":", e$message))
    return(list(coords = NULL, count = 0))
  })
}


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
  
  hlog("INFO", paste("Generating", n_background, "background points..."))
  if (is.null(predictor_stack) || !inherits(predictor_stack, "SpatRaster") || terra::nlyr(predictor_stack) == 0) { hlog("ERROR", "Predictor stack is required."); return(NULL) }
  if (!is.null(seed)) set.seed(seed)
  
  tryCatch({
    # Use first layer as mask
    mask_layer <- predictor_stack[[1]]
    bg_points <- terra::spatSample(mask_layer, size = n_background, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
    if (nrow(bg_points) < n_background) { hlog("WARN", paste("Could only sample", nrow(bg_points), "background points (requested", n_background, ").")) }
    if (nrow(bg_points) == 0) { hlog("ERROR", "Failed to generate any background points."); return(NULL) }
    hlog("DEBUG", paste("Generated", nrow(bg_points), "background points."))
    return(as.data.frame(bg_points[, c("x", "y")]))
  }, error = function(e) { hlog("ERROR", paste("Error generating background points:", e$message)); return(NULL) })
}


#' Tune SDM Hyperparameters using SDMtune gridSearch with k-fold CV
#' Returns the tuning object containing results and models.
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @return An `SDMtune` object containing tuning results, or NULL on error.
run_sdm_tuning_kfold <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species", species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[TuningHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}# cat(msg, "\n")} }
  
  hlog("INFO", paste("Tuning hyperparameters for", species_name, "using", config$sdm_n_folds, "-fold CV..."))
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df)) { hlog("ERROR", "Invalid inputs."); return(NULL) }
  
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE) }, error = function(e) { hlog("ERROR", paste("Failed prepareSWD:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  folds <- tryCatch({ SDMtune::randomFolds(data = full_swd_data, k = config$sdm_n_folds, only_presence = FALSE, seed = 123)}, error = function(e){ hlog("ERROR", paste("Failed create k-folds:", e$message)); return(NULL)})
  if(is.null(folds)) return(NULL)
  
  hlog("DEBUG", "Training initial CV base model for gridSearch...")
  initial_cv_model <- tryCatch({ SDMtune::train(method = config$sdm_method, data = full_swd_data, folds = folds, progress = FALSE)}, error = function(e) { hlog("ERROR", paste("Failed train CV base model:", e$message)); return(NULL)})
  if (is.null(initial_cv_model) || !inherits(initial_cv_model, "SDMmodelCV")) { hlog("ERROR", "Failed create valid SDMmodelCV object."); return(NULL)}
  
  hyper_grid <- config$sdm_tune_grid; tuning_results <- NULL
  tryCatch({
    hlog("DEBUG", "Running SDMtune::gridSearch with k-fold CV...")
    show_progress <- FALSE; if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::log_level(config$log_level); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE})}
    
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


#' Train Final SDM Model
#' Returns the trained model object.
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @return An `SDMmodel` object, or NULL on error.
train_final_sdm <- function(occs_coords, predictor_stack, background_df, best_hypers, config, logger, species_name = "species", species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[TrainHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  hlog("INFO", paste("Training final model for", species_name, "on full dataset..."))
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df) || is.null(best_hypers)) { hlog("ERROR", "Invalid inputs."); return(NULL) }
  
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE)}, error = function(e) { hlog("ERROR", paste("Failed prepareSWD:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  train_args <- list(method = config$sdm_method, data = full_swd_data, progress = FALSE)
  for (h_name in names(best_hypers)) { h_value <- best_hypers[[h_name]]; if (h_name == "reg") h_value <- as.numeric(h_value); train_args[[h_name]] <- h_value }
  hlog("DEBUG", paste("  Using final hyperparameters:", paste(names(train_args)[-(1:2)], sapply(train_args[-(1:2)], as.character), collapse=", ")))
  
  final_model <- tryCatch({ do.call(SDMtune::train, train_args)}, error = function(e) { hlog("ERROR", paste("Failed train final model:", e$message)); return(NULL)})
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Final model training returned invalid object."); return(NULL) }
  
  hlog("INFO", paste("  Final model trained successfully for", species_name))
  return(final_model)
}


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
  base_filename <- paste0("mean_pred_", species_name_sanitized)
  
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


#' Predict SDM Suitability (Returns Raster or Error Message)
#' @return A SpatRaster object with the prediction, or a character string with the error message on failure.
predict_sdm_suitability <- function(final_sdm_model, predictor_stack, config, logger, species_log_file = NULL, output_type = "cloglog") {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[PredictHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  hlog("DEBUG", paste("Attempting prediction (type:", output_type, ")..."))
  if (is.null(final_sdm_model) || !inherits(final_sdm_model, "SDMmodel")) { msg <- "Invalid final_sdm_model."; hlog("ERROR", msg); return(msg) }
  if (is.null(predictor_stack)) { msg <- "Predictor stack required."; hlog("ERROR", msg); return(msg) }
  
  prediction_result <- NULL
  tryCatch({
    prediction_result <- SDMtune::predict(object = final_sdm_model, data = predictor_stack, type = output_type, clamp = TRUE)
    if (is.null(prediction_result) || !inherits(prediction_result, "SpatRaster") || terra::nlyr(prediction_result) == 0) { msg <- "Prediction returned NULL or empty raster."; hlog("WARN", msg); return(msg) }
    names(prediction_result) <- "suitability"; hlog("DEBUG", "Prediction raster generated."); return(prediction_result)
  }, error = function(e) { err_msg <- paste("SDMtune::predict failed:", e$message); hlog("ERROR", err_msg); return(err_msg) })
}

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
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveTuneHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  if (is.null(tuning_output) || !inherits(tuning_output, "SDMtune")) { hlog("ERROR", "Invalid tuning object provided."); return(FALSE) }
  
  # Determine target subdirectory
  subdir_name <- config$model_output_subdir_map[[predictor_type_suffix]]
  if (is.null(subdir_name)) { hlog("ERROR", paste("No output subdirectory mapping found for suffix:", predictor_type_suffix)); return(FALSE) }
  target_subdir <- file.path(config$target_results_base, subdir_name)
  dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
  
  # Construct filenames matching the target repo style (using original species short codes if available, else sanitized name)
  # NOTE: Using sanitized name for now, post-analysis script will need to adapt.
  target_base_name <- paste0("CV_Results_", species_name_sanitized)
  
  # 1. Save full tuning object (RDS - for potential internal reuse)
  rds_file <- file.path(target_subdir, paste0(target_base_name, "_tuning_object.rds")) # Different name to avoid clash
  tryCatch({ saveRDS(tuning_output, file = rds_file); hlog("DEBUG", "Full tuning object saved (RDS):", basename(rds_file)) },
           error = function(e) { hlog("ERROR", paste("Failed save tuning RDS:", e$message)) }) # Log error but continue
  
  # 2. Save results table (CSV - for analysis mirroring target)
  csv_file <- file.path(target_subdir, paste0(target_base_name, ".csv"))
  results_df <- tuning_output@results
  if (is.null(results_df) || nrow(results_df) == 0) { hlog("WARN", "No results table found in tuning object."); return(TRUE) } # Return TRUE as RDS might have saved
  
  tryCatch({ readr::write_csv(results_df, file = csv_file); hlog("DEBUG", "Tuning results table saved (CSV):", basename(csv_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save tuning CSV:", e$message)); FALSE })
}


#' Save Final SDM Model Object
#' Saves the trained SDMmodel object as an RDS file.
#' NOTE: This saves to an *intermediate* location (`config$models_dir`), not the target analysis structure.
#'
#' @param final_model The trained `SDMmodel` object.
#' @param species_name_sanitized Sanitized species name.
#' @param predictor_type_suffix Suffix indicating model type.
#' @param config The configuration list (needs `models_dir`).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_final_model <- function(final_model, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveModelHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Invalid model object provided."); return(FALSE) }
  
  # Create model type subdirectory if needed
  model_subdir <- file.path(config$models_dir, paste0(basename(config$anemone_occurrence_dir), predictor_type_suffix)) # Example based on group occurrence dir name
  dir.create(model_subdir, recursive = TRUE, showWarnings = FALSE)
  
  model_file <- file.path(model_subdir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  tryCatch({ saveRDS(final_model, file = model_file); hlog("DEBUG", "Final model object saved (RDS):", basename(model_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save final model RDS:", e$message)); FALSE })
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
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SavePredHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  if (is.null(prediction_raster) || !inherits(prediction_raster, "SpatRaster")) { hlog("ERROR", "Invalid prediction raster provided."); return(FALSE) }
  
  # Determine target directory and filename structure
  target_dir <- NULL
  target_filename <- NULL
  
  # Construct base filename part (e.g., mean_pred_SPECIES)
  # NOTE: Using sanitized name. Adjust post-analysis or add mapping if target uses codes.
  base_filename <- paste0("mean_pred_", species_name_sanitized) # Using SDMtune mean convention implicitly
  
  if (scenario_name == "current") {
    target_dir <- config$target_predictions_current_dir
    # Filename format: mean_pred_SPECIES[_MODELTYPE?].tif
    target_filename <- paste0(base_filename, predictor_type_suffix, ".tif") # Add suffix to distinguish model types
  } else {
    # Future scenario
    ssp_dir_name <- config$ssp_scenario_map[[scenario_name]]
    if (is.null(ssp_dir_name)) { hlog("ERROR", paste("No SSP directory mapping found for scenario:", scenario_name)); return(FALSE) }
    target_dir <- file.path(config$target_predictions_future_base, ssp_dir_name)
    
    # Extract time tag for filename
    time_tag_clean <- ifelse(grepl("2050$", scenario_name), "dec50", ifelse(grepl("2100$", scenario_name), "dec100", "UNKNOWN"))
    if(time_tag_clean == "UNKNOWN") { hlog("ERROR", paste("Cannot extract time tag from scenario:", scenario_name)); return(FALSE) }
    
    # Filename format: mean_pred_SPECIES_SSP_TIME[_MODELTYPE?].tif
    target_filename <- paste0(base_filename, "_", ssp_dir_name, "_", time_tag_clean, predictor_type_suffix, ".tif") # Add suffix
  }
  
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  pred_file_path <- file.path(target_dir, target_filename)
  
  tryCatch({
    terra::writeRaster(prediction_raster, filename = pred_file_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    hlog("DEBUG", paste("Prediction raster saved:", basename(pred_file_path), "to", target_dir))
    TRUE
  }, error = function(e) {
    hlog("ERROR", paste("Failed save prediction raster:", e$message))
    FALSE
  })
}


#' Calculate and Save Variable Importance (v8 - Requires SDMmodel & SWD)
#' Calculates permutation importance using SDMtune::varImp.
#' REQUIRES a single trained SDMmodel object and the corresponding SWD object
#' used for training (passed to the 'test' argument for permutation).
#' Saves to the target structure: data/output/vi_[groupname]/vi_[species]_[suffix].csv
#'
#' @param final_model An `SDMmodel` object (output from `train_final_sdm` or `getBestModel`).
#' @param training_swd An `SWD` object containing the full training data (presence + background)
#'                     used for training the `final_model`. THIS IS REQUIRED.
#' @param species_name_sanitized Sanitized species name.
#' @param group_name Group name ("anemone" or "anemonefish").
#' @param predictor_type_suffix Suffix indicating model type. Used in filename.
#' @param config The configuration list.
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
#' @export
calculate_and_save_vi <- function(final_model, training_swd, # Both are now required
                                  species_name_sanitized, group_name, predictor_type_suffix,
                                  config, logger, species_log_file = NULL) {
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[VarImpHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  # --- Input validation ---
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) {
    hlog("ERROR", "Invalid SDMmodel object provided for VI.")
    return(FALSE)
  }
  if (is.null(training_swd) || !inherits(training_swd, "SWD")) {
    hlog("ERROR", "Training SWD object required for varImp calculation.")
    return(FALSE)
  }
  # --- End Input Validation ---
  
  hlog("INFO", "Calculating variable importance (permutation)...")
  vi_results <- NULL
  tryCatch({
    
    # Determine if progress should be shown
    show_progress <- FALSE
    if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::log_level(config$log_level %||% "INFO"); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE}) }
    
    # --- Call varImp on the SDMmodel, providing test data (training SWD) ---
    hlog("DEBUG", "Calling varImp on the provided SDMmodel object...")
    metric_vi <- config$sdm_evaluation_metric %||% "auc"
    if(tolower(metric_vi) == "aicc") metric_vi <- "auc" # Use AUC if AICc was tuning metric
    
    vi_results <- SDMtune::varImp(
      model = final_model,        # Use the single model object
      permut = 10,                # Number of permutations
      progress = show_progress
    )
    # --- End varImp call ---
    
    if (is.null(vi_results) || nrow(vi_results) == 0) { hlog("WARN", "Variable importance calculation returned empty results."); return(TRUE)}
    
    # --- Saving Logic ---
    vi_subdir_name <- paste0("vi_", group_name)
    target_subdir <- file.path(config$target_vi_base, vi_subdir_name)
    dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
    vi_filename <- paste0("vi_", species_name_sanitized, predictor_type_suffix, ".csv")
    vi_file_path <- file.path(target_subdir, vi_filename)
    readr::write_csv(vi_results, vi_file_path)
    hlog("INFO", paste("Variable importance saved to:", vi_file_path))
    return(TRUE)
    
  }, error = function(e) {
    hlog("ERROR", paste("Variable importance calculation/saving failed:", e$message))
    # print(rlang::trace_back()) # Uncomment for deeper debugging
    return(FALSE)
  })
}

#-------------------------------------------------------------------------