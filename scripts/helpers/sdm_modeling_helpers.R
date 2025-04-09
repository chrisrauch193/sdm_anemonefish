# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow)
# - Added logger argument to functions (and safe handling)
# - Modified occurrence loading to return counts
# - Corrected progress bar logic based on logger level strings
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stats, log4r) # Added log4r

#' Load, Clean, and Get Coordinates for Individual Species Occurrences
#' (Returns list with data and count)
#'
#' @param species_aphia_id AphiaID of the target species.
#' @param occurrence_dir Directory containing individual species CSV files.
#' @param config Project configuration list (potentially including `predictor_stack_for_thinning`).
#' @param logger A log4r logger object (can be NULL).
#' @return A list containing cleaned coordinates `coords` (matrix) and the `count`, or NULL on error.
load_clean_individual_occ_coords <- function(species_aphia_id, occurrence_dir, config, logger) {
  # Use logger safely
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  
  log_debug(paste("Loading occurrences for AphiaID:", species_aphia_id))
  occ_file <- file.path(occurrence_dir, paste0(species_aphia_id, ".csv"))
  if (!file.exists(occ_file)) {
    log_warn(paste("Occurrence file not found for AphiaID:", species_aphia_id))
    return(NULL)
  }
  
  tryCatch({
    occ_df <- readr::read_csv(occ_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
    
    if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occ_df))) {
      log_warn(paste("Missing coordinate columns in file:", basename(occ_file)))
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
    log_debug(paste("  Retained", count_after_clean, "records after coordinate cleaning."))
    
    if (count_after_clean == 0) {
      log_warn(paste("No valid coordinates after cleaning for AphiaID:", species_aphia_id))
      return(NULL)
    }
    
    # --- Spatial Thinning (Optional) ---
    # Note: predictor_stack_for_thinning must be passed via the config object passed to this function
    if (!is.null(config$thinning_method) && config$thinning_method == "cell" && !is.null(config$predictor_stack_for_thinning)) {
      log_debug("  Applying cell-based thinning...")
      predictor_stack_thin <- config$predictor_stack_for_thinning
      if(is.null(predictor_stack_thin)){
        log_warn(" Predictor stack for thinning not provided in config. Skipping thinning.")
        occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      } else {
        occs_sf <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs)
        target_crs <- terra::crs(predictor_stack_thin)
        if (sf::st_crs(occs_sf) != sf::st_crs(target_crs)) {
          occs_sf <- sf::st_transform(occs_sf, crs = target_crs)
        }
        occs_spatvector <- terra::vect(occs_sf)
        # Use the first layer of the stack for cell extraction
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
            log_warn(" No valid cells found for occurrences on thinning raster. Skipping thinning.")
            occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
          }
        } else {
          log_warn(" Thinning predictor stack has no layers. Skipping thinning.")
          occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
        }
      }
      count_after_thin <- nrow(occ_thinned_coords)
      # Use main logger if available, otherwise print (for parallel context)
      msg <- paste("  Retained", count_after_thin, "records after thinning for AphiaID:", species_aphia_id)
      if(!is.null(logger)) log_info(msg) else cat(msg, "\n") # Print if no logger
      
      if (count_after_thin == 0) {
        log_warn("No records left after thinning.")
        return(NULL)
      }
    } else {
      log_debug("  Skipping spatial thinning or method not 'cell'.")
      occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      count_after_thin <- nrow(occ_thinned_coords) # Same as count_after_clean here
    }
    
    return(list(coords = occ_thinned_coords, count = count_after_thin))
    
  }, error = function(e) {
    log_error(paste("Error loading/cleaning occurrences for AphiaID", species_aphia_id, ":", e$message))
    return(NULL)
  })
}


#' Generate Background Points within the Raster Extent
generate_sdm_background <- function(predictor_stack, n_background, config, logger, seed = NULL) {
  # Use logger safely
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  
  log_info(paste("Generating", n_background, "background points..."))
  if (is.null(predictor_stack)) {
    log_error("Predictor stack is required to generate background points."); return(NULL)
  }
  if (!is.null(seed)) set.seed(seed)
  
  tryCatch({
    # Ensure we sample from a layer that likely has fewer NAs if possible, e.g., bathymetry
    layer_to_sample <- if ("bathymetry_mean" %in% names(predictor_stack)) predictor_stack[["bathymetry_mean"]] else predictor_stack[[1]]
    
    bg_points <- terra::spatSample(layer_to_sample, size = n_background, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
    
    if (nrow(bg_points) < n_background) { log_warn(paste("Could only sample", nrow(bg_points), "background points (requested", n_background, ")."))}
    if (nrow(bg_points) == 0) {log_error("Failed to generate any background points."); return(NULL)}
    log_debug(paste("Generated", nrow(bg_points), "background points."))
    return(as.data.frame(bg_points[, c("x", "y")]))
  }, error = function(e) {
    log_error(paste("Error generating background points:", e$message)); return(NULL)
  })
}


#' Calculate and Save Variable Importance using SDMtune::varImp
#'
#' @param sdm_model The trained SDMmodel object (from `train` or loaded from RDS).
#' @param test_data The SWD object used for testing/evaluation during tuning (if available), or the full SWD used for training.
#' @param permutations Integer, number of permutations for `varImp`.
#' @param save_path Character, the full path to save the results CSV.
#' @param logger A log4r logger object (can be NULL).
#' @return The variable importance data frame (invisibly), or NULL on error.
calculate_and_save_vi <- function(sdm_model, test_data, permutations, save_path, logger) {
  log_info <- function(...) if (!is.null(logger)) log4r::info(logger, ...)
  log_debug <- function(...) if (!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if (!is.null(logger)) log4r::error(logger, ...)
  
  if (is.null(sdm_model) || !inherits(sdm_model, "SDMmodel")) {
    log_error("Invalid SDM model provided for varImp.")
    return(NULL)
  }
  if (is.null(test_data) || !inherits(test_data, "SWD")) {
    log_error("Invalid SWD test data provided for varImp.")
    return(NULL)
  }
  if(nrow(test_data@data) == 0) {
    log_error("SWD test data has zero rows for varImp.")
    return(NULL)
  }
  # Ensure test_data has presence/absence
  if(is.null(test_data@pa)) {
    log_error("SWD test data is missing presence/absence vector (@pa).")
    return(NULL)
  }
  # Ensure test data matches model variables (basic check)
  model_vars <- sdm_model@data@data |> colnames() |> setdiff(c("species", "pa"))
  test_vars <- test_data@data |> colnames() |> setdiff(c("species", "pa"))
  if (!all(model_vars %in% test_vars)) {
    log_error("Test data variables do not match model variables for varImp.")
    return(NULL)
  }
  
  log_debug(paste("Calculating variable importance with", permutations, "permutations..."))
  
  vi_results_df <- NULL
  tryCatch({
    vi_results <- SDMtune::varImp(sdm_model, perm = permutations, test = test_data, progress = FALSE)
    vi_results_df <- as.data.frame(vi_results)
    vi_results_df$Variable <- rownames(vi_results_df)
    rownames(vi_results_df) <- NULL
    vi_results_df <- vi_results_df[, c("Variable", "Permutation_importance")] # Keep relevant columns
    
    # Save the results
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(vi_results_df, save_path)
    log_debug(paste("Variable importance saved to:", basename(save_path)))
    
  }, error = function(e) {
    log_error(paste("Variable importance calculation failed:", e$message))
    vi_results_df <- NULL # Ensure NULL is returned on error
  })
  
  return(invisible(vi_results_df))
}


#' Tune SDM Hyperparameters using SDMtune gridSearch with k-fold CV
run_sdm_tuning_kfold <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species") {
  # ... (setup and SWD prep as before) ...
  
  # --- ADDED: Check if tuning results file exists ---
  tuning_results_file <- file.path(config$results_dir, # Use YOUR results dir
                                   paste0(species_name, "_", config$sdm_method, "_tuning.rds")) # Consistent naming
  if (!config$force_rerun$run_standard_sdms && file.exists(tuning_results_file)) {
    log_info("Loading existing tuning results from:", basename(tuning_results_file))
    tuning_results <- tryCatch(readRDS(tuning_results_file), error = function(e){
      log_warn(paste("Could not load existing tuning results, re-running tuning:", e$message)); NULL
    })
    if (!is.null(tuning_results) && inherits(tuning_results, "SDMtune")) {
      # Extract best hypers from loaded object
      res_df <- tuning_results@results
      # ... (rest of best hyper extraction logic from original function) ...
      if(length(best_row_index) > 0) {
        best_hypers_df <- res_df[best_row_index, names(config$sdm_tune_grid), drop = FALSE]
        log_info(paste("  Loaded Best hypers:", paste(names(best_hypers_df), best_hypers_df[1,], collapse=", ")))
        # Return BOTH the loaded object AND the best hypers
        return(list(best_hypers = best_hypers_df, tuning_results_object = tuning_results))
      } else {
        log_warn("Could not determine best model from loaded tuning results. Re-running tuning.")
      }
    }
  }
  # --- END ADDED Check ---
  
  # ... (run folds, initial model, gridSearch as before) ...
  
  tryCatch({
    # ... (gridSearch call) ...
    
    # ... (best hyper extraction logic) ...
    
    # --- MODIFIED: Save the full tuning object ---
    save_tuning_result <- tryCatch({
      saveRDS(tuning_results, tuning_results_file) # Save the whole object
      slog("DEBUG", "Tuning results object saved to:", basename(tuning_results_file)) # Use slog if available
      NULL
    }, error = function(e){
      slog("ERROR", "Failed save tuning results object:", e$message)
      return(list(status = "error_saving_tuning_results", ...)) # Return list with status
    })
    if (!is.null(save_tuning_result)) return(save_tuning_result)
    # --- END MODIFIED Save ---
    
    # Return BOTH the best hypers and the full tuning object
    return(list(best_hypers = best_hypers_df, tuning_results_object = tuning_results))
    
  }, error = function(e) { ... ; return(NULL) }) # Return NULL on gridSearch error
}

# --- Modify train_final_sdm to RETURN the model ---
train_final_sdm <- function(...) {
  # ... (preparation as before) ...
  
  final_model <- tryCatch({ do.call(SDMtune::train, train_args) },
                          error = function(e) { log_error(paste("Failed train final model:", e$message)); return(NULL) })
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) {
    log_error("Final model training returned invalid object."); return(NULL)
  }
  
  log_info(paste("  Final model trained successfully for", species_name))
  return(final_model) # <<< RETURN the model object
}



#' Train Final SDM Model
train_final_sdm <- function(occs_coords, predictor_stack, background_df, best_hypers, config, logger, species_name = "species") {
  # Use logger safely
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  
  log_info(paste("Training final model for", species_name, "on full dataset..."))
  
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df) || is.null(best_hypers)) {
    log_error("Invalid inputs for training final SDM."); return(NULL)
  }
  
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE)},
                            error = function(e) { log_error(paste("Failed prepare SWD for final model:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  train_args <- list(method = config$sdm_method, data = full_swd_data, progress = FALSE)
  for (h_name in names(best_hypers)) {
    h_value <- best_hypers[[h_name]]; if (h_name == "reg") h_value <- as.numeric(h_value); train_args[[h_name]] <- h_value
  }
  log_debug(paste("  Using final hyperparameters:", paste(names(train_args)[-(1:2)], sapply(train_args[-(1:2)], as.character), collapse=", ")))
  
  final_model <- tryCatch({ do.call(SDMtune::train, train_args)},
                          error = function(e) { log_error(paste("Failed train final model:", e$message)); return(NULL)})
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { log_error("Final model training returned invalid object."); return(NULL) }
  
  log_info(paste("  Final model trained successfully for", species_name))
  return(final_model)
}


#' Predict SDM Suitability (Reverted to SDMtune::predict with error handling)
#' Returns raster on success, or the error message string on failure.
#'
#' @param final_sdm_model An SDMmodel object (output from `train_final_sdm`).
#' @param predictor_stack SpatRaster object for the prediction scenario.
#' @param config Configuration list.
#' @param logger A log4r logger object (can be NULL).
#' @param output_type Character string for prediction type (e.g., "cloglog", "logistic"). Defaults to "cloglog".
#' @return A SpatRaster object with the prediction, or a character string with the error message on failure.
predict_sdm_suitability <- function(final_sdm_model, predictor_stack, config, logger, output_type = "cloglog") {
  log_debug <- function(...) if (!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if (!is.null(logger)) log4r::error(logger, ...)
  
  log_debug(paste("Attempting prediction using SDMtune::predict (type:", output_type, ")..."))
  
  if (is.null(final_sdm_model) || !inherits(final_sdm_model, "SDMmodel")) {
    msg <- "Invalid final_sdm_model object."
    if(!is.null(logger)) log_error(msg); return(msg)
  }
  if (is.null(predictor_stack) || !inherits(predictor_stack, "SpatRaster") || terra::nlyr(predictor_stack) == 0) {
    msg <- "Invalid predictor_stack provided."
    if(!is.null(logger)) log_error(msg); return(msg)
  }
  
  prediction_result <- NULL
  tryCatch({
    prediction_result <- SDMtune::predict(
      object = final_sdm_model,
      data = predictor_stack,
      type = output_type,
      clamp = TRUE # Or get from config: config$sdm_clamp
    )
    
    if (is.null(prediction_result) || !inherits(prediction_result, "SpatRaster") || terra::nlyr(prediction_result) == 0) {
      msg <- "SDMtune::predict returned NULL or empty SpatRaster."
      if(!is.null(logger)) log_warn(msg); return(msg)
    }
    
    # Rename layer for consistency if needed, but might be overwritten by writeRaster filename
    # names(prediction_result) <- "suitability"
    log_debug("  SDMtune::predict prediction raster generated successfully.")
    return(prediction_result) # Return the SpatRaster
    
  }, error = function(e) {
    err_msg <- paste("SDMtune::predict failed:", e$message)
    if(!is.null(logger)) log_error(err_msg)
    return(err_msg) # Return the error message string
  })
}


#' Plot Pearson Correlation using corrplot (Corrected Argument Handling v2)
plot_correlation_results_original <- function(env_extract, save_path = NULL, config, display_lookup = NULL) {
  if(!is.data.frame(env_extract) && !is.matrix(env_extract) || ncol(env_extract) < 2) {
    warning("plot_correlation_results_original requires df/matrix with >= 2 columns.", call.=FALSE); return(invisible(NULL))
  }
  
  # Use provided lookup first, then fallback to config lookup
  lookup_to_use <- if (!is.null(display_lookup)) display_lookup else config$core_var_display_names
  get_display_name_func <- if (!is.null(display_lookup)) config$get_display_name else config$get_display_name
  
  if(is.null(lookup_to_use) || is.null(get_display_name_func) || !is.function(get_display_name_func)){
    warning("Display name lookup or function missing/invalid for Corr plot.", call.=FALSE)
    get_display_name_func <- function(name, ...) name # Fallback
    lookup_to_use <- NULL
  }
  
  env.cor <- tryCatch(stats::cor(env_extract, method = "pearson", use = "pairwise.complete.obs"), error = function(e) {warning("Correlation calculation failed: ", e$message, call.=FALSE); NULL})
  if(is.null(env.cor)) return(invisible(NULL))
  
  env.p <- NULL
  if(requireNamespace("Hmisc", quietly = TRUE)){
    numeric_cols <- sapply(env_extract, is.numeric)
    if(sum(numeric_cols) < 2) { warning("Need >= 2 numeric columns for Hmisc::rcorr.", call.=FALSE) }
    else {
      env_extract_num <- as.matrix(env_extract[, numeric_cols, drop = FALSE])
      n_obs <- crossprod(!is.na(env_extract_num))
      if(any(n_obs < 3)) { warning("Insufficient non-NA pairs, p-values from Hmisc may be unreliable.", call.=FALSE) }
      hmisc_result <- tryCatch(Hmisc::rcorr(env_extract_num, type="pearson"), error=function(e){warning("Hmisc::rcorr failed: ",e$message, call.=FALSE); NULL})
      if(!is.null(hmisc_result)) env.p <- hmisc_result$P
    }
  } else { warning("Hmisc package not found, p-values missing.", call.=FALSE) }
  
  original_names <- colnames(env.cor)
  display_names <- sapply(original_names, get_display_name_func, lookup = lookup_to_use, USE.NAMES = FALSE)
  rownames(env.cor) <- display_names
  colnames(env.cor) <- display_names
  
  if (!is.null(env.p)) {
    if ( all(dim(env.p) == sum(numeric_cols)) ) {
      p_original_names <- colnames(env_extract)[numeric_cols]
      p_display_names <- sapply(p_original_names, get_display_name_func, lookup = lookup_to_use, USE.NAMES = FALSE)
      rownames(env.p) <- p_display_names
      colnames(env.p) <- p_display_names
      env.p <- env.p[display_names, display_names]
    } else {
      warning("Dimension mismatch between correlation matrix and p-value matrix. Disabling p-values.", call.=FALSE)
      env.p <- NULL
    }
  }
  
  tryCatch({
    plot_obj <- function() {
      corrplot::corrplot(
        corr = env.cor, method = "color", type = "upper", order = "hclust",
        p.mat = env.p, sig.level = c(.01, .05), insig = "label_sig", pch.cex = 1.5,
        pch.col = "grey", tl.col = "black", tl.srt = 45, diag = FALSE,
        na.label = "NA", mar = c(0, 0, 1.5, 0), # Adjusted margin slightly
        title = "Pearson Correlation Matrix"
      )
    }
    if (!is.null(save_path)) {
      plot_dim <- max(8, 4 + 0.3 * ncol(env.cor))
      grDevices::png(filename = save_path, width = plot_dim, height = plot_dim, units = "in", res = 300)
      plot_obj()
      grDevices::dev.off()
      cat("  Correlation plot saved to", save_path, "\n")
    } else { plot_obj() }
    return(invisible(NULL))
  }, error = function(e){ warning("Failed to create or save correlation plot: ", e$message, call.=FALSE); return(NULL)})
}
#-------------------------------------------------------------------------------