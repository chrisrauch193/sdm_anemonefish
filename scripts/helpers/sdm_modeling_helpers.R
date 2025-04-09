# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow)
# v3: Added slog argument passing for species-specific logging within helpers.
#-------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(SDMtune, terra, sf, dplyr, readr, tools, stringr, stats, log4r, Hmisc, corrplot, ggplot2) # Ensure all needed packages are listed

#' Load, Clean, and Get Coordinates for Individual Species Occurrences
#' (Returns list with data and count)
#'
#' @param species_aphia_id AphiaID of the target species.
#' @param occurrence_dir Directory containing individual species CSV files.
#' @param config Project configuration list (must include `occurrence_crs`, `min_occurrences_sdm`, `thinning_method`, `predictor_stack_for_thinning` if thinning).
#' @param logger A log4r logger object (for general script progress, can be NULL).
#' @param slog A function for species-specific logging (can be NULL).
#' @return A list containing cleaned coordinates `coords` (matrix) and the `count`, or NULL on error.
load_clean_individual_occ_coords <- function(species_aphia_id, occurrence_dir, config, logger, slog = NULL) {
  # --- Safe Logging ---
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  slog_warn <- function(...) if (!is.null(slog) && is.function(slog)) slog("WARN", ...) else log_warn(...)
  slog_info <- function(...) if (!is.null(slog) && is.function(slog)) slog("INFO", ...) else log_info(...)
  # --- End Safe Logging ---
  
  slog_debug(paste("Loading occurrences for AphiaID:", species_aphia_id))
  occ_file <- file.path(occurrence_dir, paste0(species_aphia_id, ".csv"))
  if (!file.exists(occ_file)) {
    slog_warn(paste("Occurrence file not found:", basename(occ_file)))
    return(NULL)
  }
  
  tryCatch({
    occ_df <- readr::read_csv(occ_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
    
    if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occ_df))) {
      slog_warn(paste("Missing coordinate columns in:", basename(occ_file)))
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
    
    count_after_clean <- nrow(occ_clean)
    slog_debug(paste("  Retained", count_after_clean, "records after coordinate cleaning."))
    
    if (count_after_clean == 0) {
      slog_warn("No valid coordinates after cleaning.")
      return(NULL)
    }
    
    # --- Spatial Thinning (Optional) ---
    occ_thinned_coords <- NULL
    count_after_thin <- 0
    if (!is.null(config$thinning_method) && config$thinning_method == "cell" && !is.null(config$predictor_stack_for_thinning)) {
      slog_debug("  Applying cell-based thinning...")
      predictor_stack_thin <- config$predictor_stack_for_thinning
      if(is.null(predictor_stack_thin) || !inherits(predictor_stack_thin, "SpatRaster")){
        slog_warn(" Predictor stack for thinning invalid. Skipping thinning.")
        occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      } else {
        tryCatch({
          occs_sf <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs)
          target_crs <- terra::crs(predictor_stack_thin)
          if (sf::st_crs(occs_sf) != sf::st_crs(target_crs)) {
            occs_sf <- sf::st_transform(occs_sf, crs = target_crs)
          }
          occs_spatvector <- terra::vect(occs_sf)
          if(terra::nlyr(predictor_stack_thin) > 0){
            occs_cells <- terra::extract(predictor_stack_thin[[1]], occs_spatvector, cells = TRUE)
            valid_cells_indices <- which(!is.na(occs_cells$cell))
            if(length(valid_cells_indices) > 0) {
              unique_cell_indices_in_valid <- which(!duplicated(occs_cells$cell[valid_cells_indices]))
              original_indices_to_keep <- valid_cells_indices[unique_cell_indices_in_valid]
              occs_thinned_sf <- occs_sf[original_indices_to_keep, ]
              occ_thinned_coords <- sf::st_coordinates(occs_thinned_sf)
              colnames(occ_thinned_coords) <- c("decimalLongitude", "decimalLatitude") # Ensure names
            } else {
              slog_warn(" No valid cells found for occurrences on thinning raster. Skipping thinning.")
              occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
            }
          } else {
            slog_warn(" Thinning predictor stack has no layers. Skipping thinning.")
            occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
          }
        }, error = function(e) {
          slog_error(" Error during cell-based thinning:", e$message);
          occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")]) # Fallback
        })
      }
      count_after_thin <- nrow(occ_thinned_coords)
      slog_info(paste("  Retained", count_after_thin, "records after thinning.")) # Log to species log
      if (count_after_thin == 0) { slog_warn("No records left after thinning."); return(NULL) }
      
    } else {
      if (!is.null(config$thinning_method) && config$thinning_method != "none") {
        slog_debug("  Thinning method '", config$thinning_method, "' not implemented or thinning disabled. Using unthinned data.")
      } else {
        slog_debug("  Skipping spatial thinning.")
      }
      occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      count_after_thin <- nrow(occ_thinned_coords) # Same as count_after_clean here
    }
    
    return(list(coords = occ_thinned_coords, count = count_after_thin))
    
  }, error = function(e) {
    slog_error(paste("Error loading/cleaning occurrences:", e$message))
    return(NULL)
  })
}


#' Generate Background Points within the Raster Extent
generate_sdm_background <- function(predictor_stack, n_background, config, logger, seed = NULL, slog = NULL) {
  # --- Safe Logging ---
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  slog_warn <- function(...) if (!is.null(slog) && is.function(slog)) slog("WARN", ...) else log_warn(...)
  # --- End Safe Logging ---
  
  log_info(paste("Generating", n_background, "background points...")) # Main log
  if (is.null(predictor_stack)) { slog_error("Predictor stack is required."); return(NULL) }
  if (!is.null(seed)) set.seed(seed)
  
  tryCatch({
    layer_to_sample <- if ("bathymetry_mean" %in% names(predictor_stack)) predictor_stack[["bathymetry_mean"]] else predictor_stack[[1]]
    bg_points <- terra::spatSample(layer_to_sample, size = n_background, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
    if (nrow(bg_points) < n_background) { slog_warn(paste("Sampled only", nrow(bg_points), "background points (requested", n_background, ").")) }
    if (nrow(bg_points) == 0) { slog_error("Failed to generate any background points."); return(NULL) }
    slog_debug(paste("Generated", nrow(bg_points), "background points.")) # Species log
    return(as.data.frame(bg_points[, c("x", "y")]))
  }, error = function(e) { slog_error(paste("Error generating background points:", e$message)); return(NULL) })
}


#' Tune SDM Hyperparameters using SDMtune gridSearch with k-fold CV
#' Returns a list containing best_hypers (dataframe) and tuning_results_object (SDMtune object), or NULL on error.
run_sdm_tuning_kfold <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species", slog = NULL) {
  # --- Safe Logging ---
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  slog_warn <- function(...) if (!is.null(slog) && is.function(slog)) slog("WARN", ...) else log_warn(...)
  slog_info <- function(...) if (!is.null(slog) && is.function(slog)) slog("INFO", ...) else log_info(...)
  # --- End Safe Logging ---
  
  slog_info(paste("Starting hyperparameter tuning using", config$sdm_n_folds, "-fold CV..."))
  
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df) || is.null(config$sdm_tune_grid)) {
    msg <- "Invalid inputs for SDM tuning (occs, stack, background, or tune_grid missing/invalid)."
    slog_error(msg); return(NULL)
  }
  
  slog_debug("Preparing SWD data for tuning...")
  full_swd_data <- tryCatch({
    SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE)
  }, error = function(e) { err_msg <- paste("Failed prepare SWD for tuning:", e$message); slog_error(err_msg); return(NULL) })
  if (is.null(full_swd_data)) return(NULL)
  slog_debug("SWD data prepared.")
  
  slog_debug("Creating k-folds...")
  folds <- tryCatch({
    SDMtune::randomFolds(data = full_swd_data, k = config$sdm_n_folds, only_presence = FALSE, seed = 123)
  }, error = function(e){ err_msg <- paste("Failed create k-folds:", e$message); slog_error(err_msg); return(NULL) })
  if(is.null(folds)) return(NULL)
  slog_debug(paste(config$sdm_n_folds, "folds created."))
  
  hyper_grid <- config$sdm_tune_grid
  tuning_results <- NULL
  tryCatch({
    slog_debug("Running SDMtune::gridSearch...")
    show_progress <- FALSE
    if (!is.null(logger)) {
      current_log_level_num <- tryCatch(log4r::log_level(config$log_level), error=function(e) log4r::INFO)
      debug_level_num <- tryCatch(log4r::log_level("DEBUG"), error=function(e) log4r::DEBUG)
      show_progress <- current_log_level_num <= debug_level_num
    }
    
    tuning_results <- SDMtune::gridSearch(
      data = full_swd_data, hypers = hyper_grid, metric = config$sdm_evaluation_metric,
      folds = folds, save_models = TRUE, progress = show_progress,
      interactive = FALSE, method = config$sdm_method
    )
    slog_debug("SDMtune::gridSearch completed.")
    
    res_df <- tuning_results@results
    if(is.null(res_df) || nrow(res_df) == 0) { msg <- "No results in tuning object."; slog_warn(msg); return(NULL) }
    
    metric_base_upper <- toupper(config$sdm_evaluation_metric)
    target_metric_col <- if (tolower(config$sdm_evaluation_metric) == "aicc") "AICc" else paste0("test_", metric_base_upper)
    
    if (!target_metric_col %in% colnames(res_df)) { msg <- paste("Metric '", target_metric_col, "' not found."); slog_error(msg); print(head(res_df)); return(NULL) }
    
    valid_metric_indices <- which(!is.na(res_df[[target_metric_col]]))
    if(length(valid_metric_indices) == 0) { msg <- paste("All values for metric '", target_metric_col, "' are NA."); slog_warn(msg); return(NULL) }
    
    select_fun <- if (target_metric_col == "AICc") which.min else which.max
    best_row_relative_index <- select_fun(res_df[[target_metric_col]][valid_metric_indices])
    best_row_index <- valid_metric_indices[best_row_relative_index]
    
    if(length(best_row_index) == 0 || best_row_index > nrow(res_df)) { slog_warn("Could not determine best model index, using first valid."); best_row_index <- valid_metric_indices[1] }
    
    best_hypers_df <- res_df[best_row_index, names(hyper_grid), drop = FALSE]
    msg <- paste("Best hypers (Mean CV", target_metric_col,"=", round(res_df[[target_metric_col]][best_row_index], 4),"): ", paste(names(best_hypers_df), best_hypers_df[1,], collapse=", "))
    slog_info(msg)
    
    return(list(best_hypers = best_hypers_df, tuning_results_object = tuning_results))
    
  }, error = function(e) { err_msg <- paste("SDMtune::gridSearch failed:", e$message); slog_error(err_msg); return(NULL) })
}


#' Train Final SDM Model (Returns Model Object)
train_final_sdm <- function(occs_coords, predictor_stack, background_df, best_hypers, config, logger, species_name = "species", slog = NULL) {
  # --- Safe Logging ---
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  slog_info <- function(...) if (!is.null(slog) && is.function(slog)) slog("INFO", ...) else log_info(...)
  # --- End Safe Logging ---
  
  slog_info("Training final model...")
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df) || is.null(best_hypers)) {
    msg <- "Invalid inputs for final training."; slog_error(msg); return(NULL)
  }
  
  slog_debug("Preparing SWD for final training...")
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE)},
                            error = function(e) { msg <- paste("Failed prepare SWD for final model:", e$message); slog_error(msg); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  slog_debug("SWD prepared for final training.")
  
  train_args <- list(method = config$sdm_method, data = full_swd_data, progress = FALSE)
  for (h_name in names(best_hypers)) {
    h_value <- best_hypers[[h_name]]; if (h_name == "reg") h_value <- as.numeric(h_value); train_args[[h_name]] <- h_value
  }
  slog_debug(paste("  Using final hyperparameters:", paste(names(train_args)[-(1:2)], sapply(train_args[-(1:2)], as.character), collapse=", ")))
  
  final_model <- tryCatch({ do.call(SDMtune::train, train_args) },
                          error = function(e) { msg <- paste("Failed train final model:", e$message); slog_error(msg); return(NULL) })
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { msg <- "Final model training returned invalid object."; slog_error(msg); return(NULL) }
  
  slog_info("Final model trained successfully.")
  return(final_model) # Return the model object
}


#' Predict SDM Suitability
#' Returns raster on success, or the error message string on failure.
predict_sdm_suitability <- function(final_sdm_model, predictor_stack, config, logger, output_type = "cloglog", slog = NULL) {
  # --- Safe Logging ---
  log_debug <- function(...) if (!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if (!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if (!is.null(logger)) log4r::warn(logger, ...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  slog_warn <- function(...) if (!is.null(slog) && is.function(slog)) slog("WARN", ...) else log_warn(...)
  # --- End Safe Logging ---
  
  slog_debug(paste("Attempting prediction (type:", output_type, ")..."))
  
  if (is.null(final_sdm_model) || !inherits(final_sdm_model, "SDMmodel")) { msg <- "Invalid final_sdm_model object."; slog_error(msg); return(msg) }
  if (is.null(predictor_stack) || !inherits(predictor_stack, "SpatRaster") || terra::nlyr(predictor_stack) == 0) { msg <- "Invalid predictor_stack provided."; slog_error(msg); return(msg) }
  
  # Ensure predictor names match the model's expected names
  expected_vars <- names(final_sdm_model@data@data) |> setdiff(c("species", "pa"))
  current_vars <- names(predictor_stack)
  if(!identical(sort(expected_vars), sort(current_vars))) {
    msg <- paste0("Predictor mismatch. Model expects: ", paste(sort(expected_vars), collapse=", "), ". Stack has: ", paste(sort(current_vars), collapse=", "))
    slog_error(msg)
    return(msg)
  }
  # Reorder stack if necessary
  if(!identical(expected_vars, current_vars)) {
    predictor_stack <- predictor_stack[[expected_vars]]
    slog_debug("Reordered prediction stack layers to match model.")
  }
  
  prediction_result <- NULL
  tryCatch({
    prediction_result <- SDMtune::predict(
      object = final_sdm_model,
      data = predictor_stack,
      type = output_type,
      clamp = TRUE
    )
    if (is.null(prediction_result) || !inherits(prediction_result, "SpatRaster") || terra::nlyr(prediction_result) == 0) { msg <- "Prediction returned NULL/empty."; slog_warn(msg); return(msg) }
    # names(prediction_result) <- "suitability" # Naming handled by writeRaster
    slog_debug("Prediction successful.")
    return(prediction_result) # Return the SpatRaster
    
  }, error = function(e) { err_msg <- paste("Prediction failed:", e$message); slog_error(err_msg); return(err_msg) })
}


#' Calculate and Save Variable Importance using SDMtune::varImp
calculate_and_save_vi <- function(sdm_model, test_data, permutations, save_path, logger, slog = NULL) {
  # --- Safe Logging ---
  log_debug <- function(...) if (!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if (!is.null(logger)) log4r::error(logger, ...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  # --- End Safe Logging ---
  
  if (is.null(sdm_model) || !inherits(sdm_model, "SDMmodel")) { slog_error("Invalid SDM model for varImp."); return(NULL) }
  if (is.null(test_data) || !inherits(test_data, "SWD")) { slog_error("Invalid SWD test data for varImp."); return(NULL) }
  if(nrow(test_data@data) == 0) { slog_error("SWD test data has zero rows."); return(NULL) }
  if(is.null(test_data@pa)) { slog_error("SWD test data missing @pa vector."); return(NULL) }
  
  model_vars <- names(sdm_model@data@data) |> setdiff(c("species", "pa"))
  test_vars <- names(test_data@data) |> setdiff(c("species", "pa"))
  if (!all(model_vars %in% test_vars)) { slog_error("Test data variables mismatch model variables for varImp."); return(NULL) }
  
  slog_debug(paste("Calculating variable importance with", permutations, "permutations..."))
  
  vi_results_df <- NULL
  tryCatch({
    vi_results <- SDMtune::varImp(sdm_model, perm = permutations, test = test_data, progress = FALSE)
    vi_results_df <- as.data.frame(vi_results)
    vi_results_df$Variable <- rownames(vi_results_df)
    rownames(vi_results_df) <- NULL
    vi_results_df <- vi_results_df[, c("Variable", "Permutation_importance")]
    
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(vi_results_df, save_path)
    slog_debug(paste("Variable importance saved to:", basename(save_path)))
    
  }, error = function(e) { slog_error(paste("VI calculation failed:", e$message)); vi_results_df <- NULL })
  return(invisible(vi_results_df))
}


# --- Plotting Functions (Kept from previous version, updated to use slog_safe/config) ---
plot_vif_results_original <- function(vif_result, save_path = NULL, config, display_lookup = NULL, logger = NULL, slog = NULL) {
  # --- Safe Logging ---
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  slog_warn <- function(...) if (!is.null(slog) && is.function(slog)) slog("WARN", ...) else log_warn(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  # --- End Safe Logging ---
  if(!is.numeric(vif_result) || is.null(names(vif_result)) || length(vif_result) == 0){slog_warn("plot_vif_results_original: invalid input."); return(NULL)}
  if(is.null(config$core_var_display_names) || is.null(config$get_display_name) || !is.function(config$get_display_name)){slog_warn("Config missing display names/function for VIF plot."); get_display_name_func <- function(name, ...) name; lookup_to_use <- NULL} else {get_display_name_func <- config$get_display_name; lookup_to_use <- display_lookup %||% config$core_var_display_names}
  
  original_names <- names(vif_result); display_names <- sapply(original_names, get_display_name_func, lookup = lookup_to_use, USE.NAMES = FALSE)
  df <- data.frame(OriginalName = original_names, DisplayName = display_names, VIF = as.numeric(vif_result)); num_vars <- nrow(df)
  if (num_vars == 0) return(NULL); df$VIF <- pmax(0, df$VIF); df$DisplayName <- factor(df$DisplayName, levels = df$DisplayName[order(df$VIF)])
  vif_threshold_val <- config$vif_threshold %||% 5 # Default to 5 if not in config
  p <- ggplot2::ggplot(df, aes(x = DisplayName, y = VIF)) + ggplot2::geom_bar(stat = "identity", aes(fill = OriginalName), show.legend = FALSE) + ggplot2::scale_fill_viridis_d(option = "plasma") + ggplot2::labs(title = "VIF Analysis Results (car::vif)", x = "Environmental Variables", y = "VIF Value") + ggplot2::scale_y_continuous(limits = c(0, max(vif_threshold_val, ceiling(max(df$VIF, na.rm = TRUE)))), breaks = scales::pretty_breaks(n = 5)) + ggplot2::geom_hline(yintercept = vif_threshold_val, linetype = "dashed", color = "red") + ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(save_path)) {tryCatch({ggplot2::ggsave(save_path, plot = p, width = max(8, 4 + 0.2 * num_vars), height = 6, limitsize = FALSE); slog_debug("VIF plot saved:", basename(save_path))}, error = function(e) { slog_warn("Failed save VIF plot:", e$message)})}
  return(p)
}

plot_correlation_results_original <- function(env_extract, save_path = NULL, config, display_lookup = NULL, logger = NULL, slog = NULL) {
  # --- Safe Logging ---
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  slog_warn <- function(...) if (!is.null(slog) && is.function(slog)) slog("WARN", ...) else log_warn(...)
  slog_error <- function(...) if (!is.null(slog) && is.function(slog)) slog("ERROR", ...) else log_error(...)
  slog_debug <- function(...) if (!is.null(slog) && is.function(slog)) slog("DEBUG", ...) else log_debug(...)
  # --- End Safe Logging ---
  if(!is.data.frame(env_extract) && !is.matrix(env_extract) || ncol(env_extract) < 2) {slog_warn("plot_correlation_results_original: requires df/matrix with >= 2 columns."); return(invisible(NULL))}
  if(is.null(config$core_var_display_names) || is.null(config$get_display_name) || !is.function(config$get_display_name)){slog_warn("Config missing display names/function for Corr plot."); get_display_name_func <- function(name, ...) name; lookup_to_use <- NULL} else {get_display_name_func <- config$get_display_name; lookup_to_use <- display_lookup %||% config$core_var_display_names}
  
  env.cor <- tryCatch(stats::cor(env_extract, method = "pearson", use = "pairwise.complete.obs"), error = function(e) {slog_error("Corr calc failed: ", e$message); NULL})
  if(is.null(env.cor)) return(invisible(NULL))
  env.p <- NULL
  if(requireNamespace("Hmisc", quietly = TRUE)){ numeric_cols <- sapply(env_extract, is.numeric); if(sum(numeric_cols) < 2) { slog_warn("Need >= 2 numeric cols for Hmisc::rcorr.") } else { env_extract_num <- as.matrix(env_extract[, numeric_cols, drop = FALSE]); n_obs <- crossprod(!is.na(env_extract_num)); if(any(n_obs < 3)) { slog_warn("Insufficient pairs, p-values may be unreliable.") }; hmisc_result <- tryCatch(Hmisc::rcorr(env_extract_num, type="pearson"), error=function(e){slog_warn("Hmisc::rcorr failed: ",e$message); NULL}); if(!is.null(hmisc_result)) env.p <- hmisc_result$P }} else { slog_warn("Hmisc not found, p-values missing.") }
  original_names <- colnames(env.cor); display_names <- sapply(original_names, get_display_name_func, lookup = lookup_to_use, USE.NAMES = FALSE)
  rownames(env.cor) <- display_names; colnames(env.cor) <- display_names
  if (!is.null(env.p)) { if ( all(dim(env.p) == sum(numeric_cols)) ) { p_original_names <- colnames(env_extract)[numeric_cols]; p_display_names <- sapply(p_original_names, get_display_name_func, lookup = lookup_to_use, USE.NAMES = FALSE); rownames(env.p) <- p_display_names; colnames(env.p) <- p_display_names; env.p <- env.p[display_names, display_names] } else { slog_warn("Dim mismatch corr/p-value matrix. Disabling p-values."); env.p <- NULL }}
  
  tryCatch({ plot_obj <- function() { corrplot::corrplot(corr = env.cor, method="color", type = "upper", order = "hclust", p.mat = env.p, sig.level = c(.01, .05), insig = "label_sig", pch.cex = 1.5, pch.col = "grey", tl.col = "black", tl.srt = 45, diag = FALSE, na.label = "NA", mar=c(0,0,1,0), title = "Pearson Correlation Matrix") }
  if (!is.null(save_path)) { plot_dim <- max(8, 4 + 0.3 * ncol(env.cor)); grDevices::png(filename = save_path, width = plot_dim, height = plot_dim, units = "in", res = 300); plot_obj(); grDevices::dev.off(); slog_debug("Correlation plot saved:", basename(save_path)) } else { plot_obj() }
  return(invisible(NULL))
  }, error = function(e){slog_warn("Failed create/save correlation plot: ", e$message); return(NULL)})
}
#-------------------------------------------------------------------------------