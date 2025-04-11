# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow v2)
# - Saving logic moved to dedicated helper functions.
# - Path construction uses new config structure for target output.
# - Added Variable Importance helper.
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr, blockCV) # Added readr

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


#' Generate Background Points within the Raster Extent (v2 - with Coral Masking)
#' Can optionally mask sampling to coral reef areas defined in config.
#' @param predictor_stack SpatRaster stack (used for extent/masking).
#' @param n_background Number of points to attempt generating.
#' @param config Project configuration list (needs `mask_background_points_to_coral`, `apply_coral_mask`, `coral_shapefile`). #<< ADDED config
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @param seed Optional random seed for reproducibility.
#' @return A data frame of background point coordinates (x, y), or NULL on error.
generate_sdm_background <- function(predictor_stack, n_background, config, logger, species_log_file = NULL, seed = NULL) { #<< ADDED config
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[BgGenHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  hlog("DEBUG", paste("Generating up to", n_background, "background points...")) # Changed log level
  if (is.null(predictor_stack) || !inherits(predictor_stack, "SpatRaster") || terra::nlyr(predictor_stack) == 0) { hlog("ERROR", "Predictor stack is required."); return(NULL) }
  if (!is.null(seed)) set.seed(seed)
  
  # --- Determine the layer to sample from ---
  sampling_layer <- predictor_stack[[1]] # Start with the first layer
  sampling_mask <- NULL # Initialize
  
  # Apply Coral Mask (if configured in config)
  if (config$mask_background_points_to_coral && config$apply_coral_mask) {
    hlog("DEBUG", "  Attempting coral reef mask for background point sampling...")
    if (!is.null(config$coral_shapefile) && file.exists(config$coral_shapefile)) {
      coral_areas_sf <- tryCatch({ sf::st_read(config$coral_shapefile, quiet = TRUE) }, error = function(e) {hlog("WARN",paste("   Failed load coral shapefile:", e$message)); NULL})
      if (!is.null(coral_areas_sf)) {
        coral_areas_vect <- tryCatch({ terra::vect(coral_areas_sf) }, error = function(e) {hlog("WARN",paste("   Failed convert coral sf to vect:", e$message)); NULL})
        if (!is.null(coral_areas_vect)) {
          if(terra::crs(coral_areas_vect) != terra::crs(sampling_layer)){
            hlog("DEBUG", "    Projecting coral shapefile CRS...")
            coral_areas_vect <- tryCatch(terra::project(coral_areas_vect, terra::crs(sampling_layer)), error = function(e){hlog("WARN",paste("     Failed project coral shapefile:", e$message)); NULL})
          }
          if(!is.null(coral_areas_vect)){
            sampling_mask <- tryCatch(terra::mask(sampling_layer, coral_areas_vect), error=function(e){hlog("WARN",paste("     Failed mask sampling layer:", e$message)); NULL})
            if(!is.null(sampling_mask)) { hlog("DEBUG", "  Coral reef mask applied for sampling.") }
            else { hlog("WARN", "  Failed to create mask. Sampling from unmasked layer.") }
          } else { hlog("WARN", "  CRS projection failed. Sampling from unmasked layer.") }
        } else { hlog("WARN", "  sf to vect conversion failed. Sampling from unmasked layer.") }
      } else { hlog("WARN", "  Coral shapefile loading failed. Sampling from unmasked layer.") }
    } else { hlog("WARN", "  Coral shapefile path missing/invalid. Sampling from unmasked layer.") }
  } else { hlog("DEBUG", "  Background point sampling not masked to coral reefs.") }
  
  # Use the masked layer if created, otherwise the original first layer
  layer_to_sample_from <- if(!is.null(sampling_mask)) sampling_mask else sampling_layer
  
  # --- Sample Points ---
  tryCatch({
    bg_points <- terra::spatSample(layer_to_sample_from, size = n_background, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
    if (nrow(bg_points) < n_background) { hlog("WARN", paste("Could only sample", nrow(bg_points), "background points (requested", n_background, ", likely due to mask/raster extent).")) }
    if (nrow(bg_points) == 0) { hlog("ERROR", "Failed to generate ANY background points from the sampling area."); return(NULL) }
    hlog("INFO", paste("Generated", nrow(bg_points), "background points.")) # Changed log level
    return(as.data.frame(bg_points[, c("x", "y")]))
  }, error = function(e) { hlog("ERROR", paste("Error generating background points:", e$message)); return(NULL) })
}


#' Tune SDM Hyperparameters using Spatial Cross-Validation (SCV)
#'
#' Uses blockCV::spatialBlock to create spatial folds and SDMtune::gridSearch
#' to tune hyperparameters based on a specified metric. Returns the tuning object.
#'
#' @param occs_coords Matrix or data frame of presence coordinates (decimalLongitude, decimalLatitude).
#' @param predictor_stack SpatRaster stack for the tuning scenario. MUST have a defined CRS.
#' @param background_df Data frame of background coordinates (x, y).
#' @param config Project configuration list (needs sdm_n_folds, sdm_tune_grid, sdm_evaluation_metric).
#' @param logger A log4r logger object (can be NULL).
#' @param species_name Character string for logging/output.
#' @param species_log_file Optional path to a species-specific log file.
#' @return An `SDMtune` object containing tuning results (including models trained during CV),
#'         with the best hyperparameters added as an attribute ("best_hypers"), or NULL on error.
#' @export
run_sdm_tuning_scv <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species", species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SCV_TuneHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}# cat(msg, "\n")} }
  
  hlog("INFO", paste("Tuning hyperparameters for", species_name, "using Spatial CV..."))
  
  # --- Input Checks ---
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm) { hlog("ERROR", "Insufficient occurrence points."); return(NULL) }
  if (is.null(predictor_stack) || !inherits(predictor_stack, "SpatRaster") || terra::nlyr(predictor_stack) == 0) { hlog("ERROR", "Invalid predictor stack."); return(NULL) }
  if (is.null(background_df) || nrow(background_df) == 0) { hlog("ERROR", "Invalid background points."); return(NULL) }
  target_crs_terra <- terra::crs(predictor_stack)
  if (target_crs_terra == "") { hlog("ERROR", "Predictor stack CRS is missing. Cannot perform spatial CV."); return(NULL) }
  target_crs_sf <- sf::st_crs(target_crs_terra)
  if (is.null(config$sdm_n_folds) || config$sdm_n_folds < 2) { hlog("ERROR", "Invalid k-folds number in config."); return(NULL) }
  # --- End Input Checks ---
  
  # --- Prepare Data for blockCV ---
  hlog("DEBUG", "Preparing data for spatial blocking...")
  tryCatch({
    # Presence points
    pres_sf <- sf::st_as_sf(as.data.frame(occs_coords), coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs)
    if(sf::st_crs(pres_sf) != target_crs_sf) { pres_sf <- sf::st_transform(pres_sf, crs = target_crs_sf)}
    pres_sf$presence <- 1
    
    # Background points
    back_sf <- sf::st_as_sf(background_df, coords = c("x", "y"), crs = target_crs_sf) # Assume background is already in target CRS or handle transformation if needed
    back_sf$presence <- 0
    
    # Combine
    combined_sf <- rbind(pres_sf, back_sf)
  }, error = function(e) { hlog("ERROR", paste("Failed preparing sf objects for blockCV:", e$message)); return(NULL)})
  
  # --- Create Spatial Blocks ---
  hlog("DEBUG", paste("Creating spatial blocks (k =", config$sdm_n_folds, ") using blockCV..."))
  scv_blocks <- NULL
  tryCatch({
    # Determine range based on points - may need adjustment for very clustered data
    # range_val <- blockCV::rangeSuggest(speciesData = combined_sf,
    #                                   species = "presence",
    #                                   rasterLayer = predictor_stack[[1]],
    #                                   sampleNumber = 5000) # Sample subset for speed
    # hlog("DEBUG", paste("Suggested block range:", range_val))
    # If range calculation fails or is unsuitable, use a default or calculate differently
    
    scv_blocks <- blockCV::spatialBlock(
      speciesData = combined_sf,
      species = "presence", # Name of the presence/absence column
      rasterLayer = predictor_stack[[1]], # Use first layer for spatial context
      theRange = NULL, # Let blockCV determine range (or set manually if needed)
      k = config$sdm_n_folds,
      selection = "random", # or "systematic" or "checkerboard"
      iteration = 100, # Number of attempts to find good blocks
      biomod2 = FALSE, # Output format compatible with SDMtune? Needs check - Assume list format okay.
      xOffset = 0,
      yOffset = 0,
      progress = FALSE, # Show progress bar?
      showBlocks = FALSE # Set to TRUE to plot blocks for debugging
    )
    
    # Optional: Save block plot for inspection
    # block_plot_file <- file.path(config$species_log_dir, paste0(species_name, predictor_type_suffix, "_scv_blocks.png"))
    # png(block_plot_file, width=7, height=7, units="in", res=150)
    # plot(predictor_stack[[1]], main="Spatial Blocks")
    # plot(scv_blocks$blocks, add=TRUE, border="blue", lwd=2)
    # plot(combined_sf[combined_sf$presence==1, ], pch=19, cex=0.5, col="red", add=TRUE)
    # plot(combined_sf[combined_sf$presence==0, ], pch=3, cex=0.3, col="black", add=TRUE)
    # dev.off()
    # hlog("DEBUG", paste("Spatial block plot saved to:", block_plot_file))
    
  }, error = function(e) { hlog("ERROR", paste("blockCV::spatialBlock failed:", e$message)); return(NULL)})
  
  if (is.null(scv_blocks) || is.null(scv_blocks$folds)) { hlog("ERROR", "Spatial block creation did not return valid folds."); return(NULL) }
  hlog("DEBUG", "Spatial blocks created.")
  
  # --- Prepare SWD ---
  hlog("DEBUG", "Preparing SWD object...")
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE) }, error = function(e) { hlog("ERROR", paste("Failed prepareSWD:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  # --- Run Grid Search with SCV Folds ---
  hyper_grid <- config$sdm_tune_grid; tuning_results <- NULL
  tryCatch({
    hlog("DEBUG", "Running SDMtune::gridSearch with spatial CV folds...")
    show_progress <- FALSE; if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::log_level(config$log_level %||% "INFO"); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE})}
    
    tuning_results <- SDMtune::gridSearch(
      data = full_swd_data,               # Pass SWD object
      hypers = hyper_grid,
      method = config$sdm_method,
      metric = config$sdm_evaluation_metric,
      folds = scv_blocks$folds,           # Pass the spatial folds list directly
      save_models = TRUE,                 # Important to save models trained on CV folds
      progress = show_progress,
      interactive = FALSE,
      parallel = FALSE                    # Parallel is usually handled outside this function
    )
    hlog("DEBUG", "SDMtune::gridSearch completed.")
    
    res_df <- tuning_results@results
    if(is.null(res_df) || nrow(res_df) == 0) { hlog("WARN", "No results found in tuning object."); return(NULL) }
    
    metric_base_upper <- toupper(config$sdm_evaluation_metric)
    target_metric_col <- if (tolower(config$sdm_evaluation_metric) == "aicc") "AICc" else paste0("test_", metric_base_upper) # Use test_ metric for CV
    if (!target_metric_col %in% colnames(res_df)) { hlog("ERROR", paste("CV metric '", target_metric_col, "' not found. Available:", paste(colnames(res_df), collapse=", "))); print(head(res_df)); return(NULL) }
    
    valid_metric_indices <- which(!is.na(res_df[[target_metric_col]]))
    if(length(valid_metric_indices) == 0) { hlog("WARN", paste("All values for metric '", target_metric_col, "' are NA.")); return(NULL) }
    
    select_fun <- if (target_metric_col == "AICc") which.min else which.max
    best_row_relative_index <- select_fun(res_df[[target_metric_col]][valid_metric_indices])
    best_row_index <- valid_metric_indices[best_row_relative_index]
    if(length(best_row_index) == 0 || best_row_index > nrow(res_df)) { hlog("WARN", "Could not determine best model index, using first valid row."); best_row_index <- valid_metric_indices[1]}
    
    best_hypers_df <- res_df[best_row_index, names(hyper_grid), drop = FALSE]
    hlog("INFO", paste("  Best hypers (Mean SCV", target_metric_col,"=", round(res_df[[target_metric_col]][best_row_index], 4),"): ", paste(names(best_hypers_df), best_hypers_df[1,], collapse=", ")))
    
    # Add best hypers to the object attributes for easy access later
    attr(tuning_results, "best_hypers") <- best_hypers_df
    return(tuning_results) # Return the whole tuning object
    
  }, error = function(e) { hlog("ERROR", paste("SDMtune::gridSearch with SCV failed:", e$message)); return(NULL) })
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


# scripts/helpers/sdm_modeling_helpers.R

# ... (keep other functions above) ...

#' Calculate and Save Variable Importance (v9 - Handles MaxNet Directly)
#' Calculates permutation importance. For MaxNet models, it extracts the
#' importance calculated during training. For other models, it uses SDMtune::varImp.
#' Saves results to the target structure.
#'
#' @param final_model An `SDMmodel` object (output from `train_final_sdm`).
#' @param training_swd An `SWD` object (potentially needed for varImp on non-MaxNet).
#' @param species_name_sanitized Sanitized species name.
#' @param group_name Group name ("anemone" or "anemonefish").
#' @param predictor_type_suffix Suffix indicating model type. Used in filename.
#' @param config The configuration list.
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
#' @export
calculate_and_save_vi <- function(final_model, training_swd, # training_swd kept for potential use by other methods
                                  species_name_sanitized, group_name, predictor_type_suffix,
                                  config, logger, species_log_file = NULL) {
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[VarImpHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  # --- Input validation ---
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) {
    hlog("ERROR", "Invalid SDMmodel object provided for VI.")
    return(FALSE)
  }
  # Training SWD validation might only be strictly needed if varImp is called below
  # if (is.null(training_swd) || !inherits(training_swd, "SWD")) {
  #   hlog("ERROR", "Training SWD object required for varImp calculation.")
  #   return(FALSE)
  # }
  
  hlog("INFO", "Calculating/Extracting variable importance...")
  vi_results_df <- NULL
  
  tryCatch({
    model_method <- final_model@method # Get the method used (e.g., "Maxnet")
    
    if (model_method == "Maxnet") {
      hlog("DEBUG", "Maxnet model detected. Extracting pre-calculated permutation importance.")
      # Access the raw maxnet model object stored by SDMtune
      raw_maxnet_model <- final_model@model
      if (!is.null(raw_maxnet_model) && !is.null(raw_maxnet_model$variable.importance)) {
        # The importance is stored as a named numeric vector
        importance_vector <- raw_maxnet_model$variable.importance
        # Convert to a standard data frame
        vi_results_df <- data.frame(
          Variable = names(importance_vector),
          Importance = as.numeric(importance_vector),
          stringsAsFactors = FALSE
        )
        vi_results_df <- vi_results_df[order(-vi_results_df$Importance), ] # Order descending
        hlog("DEBUG", "Successfully extracted MaxNet variable importance.")
      } else {
        hlog("WARN", "Could not find pre-calculated variable importance in the MaxNet model object.")
      }
    } else {
      # Fallback to SDMtune::varImp for other methods (e.g., RF, BRT) - may need testing
      hlog("DEBUG", paste("Model method is", model_method, ". Attempting SDMtune::varImp..."))
      if (is.null(training_swd) || !inherits(training_swd, "SWD")) {
        hlog("ERROR", "Training SWD object required for varImp calculation with non-MaxNet models.")
        return(FALSE)
      }
      show_progress <- FALSE
      if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::log_level(config$log_level %||% "INFO"); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE}) }
      
      vi_results <- SDMtune::varImp(
        model = final_model,
        permut = 10,
        progress = show_progress,
        test = training_swd # Provide test data for permutation
      )
      if (!is.null(vi_results) && nrow(vi_results) > 0) {
        vi_results_df <- vi_results # Already a data frame
        hlog("DEBUG", "SDMtune::varImp successful for non-MaxNet model.")
      } else {
        hlog("WARN", "SDMtune::varImp returned empty results for non-MaxNet model.")
      }
    }
    
    # --- Saving Logic ---
    if (is.null(vi_results_df)) {
      hlog("ERROR", "Variable importance could not be obtained.")
      return(FALSE)
    }
    
    vi_subdir_name <- paste0("vi_", group_name)
    # Use the correct target base directory from config
    target_subdir <- file.path(config$target_vi_base, vi_subdir_name) # Use config$target_vi_base
    dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
    vi_filename <- paste0("vi_", species_name_sanitized, predictor_type_suffix, ".csv")
    vi_file_path <- file.path(target_subdir, vi_filename)
    
    readr::write_csv(vi_results_df, vi_file_path)
    hlog("INFO", paste("Variable importance saved to:", vi_file_path))
    return(TRUE)
    
  }, error = function(e) {
    hlog("ERROR", paste("Variable importance processing/saving failed:", e$message))
    # print(rlang::trace_back()) # Uncomment for deeper debugging
    return(FALSE)
  })
}

#-------------------------------------------------------------------------