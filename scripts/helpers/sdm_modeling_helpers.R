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

#' Tune SDM Hyperparameters using SDMtune gridSearch with k-fold CV
run_sdm_tuning_kfold <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species") {
  # Use logger safely
  log_info <- function(...) if(!is.null(logger)) log4r::info(logger, ...)
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  
  log_info(paste("Tuning hyperparameters for", species_name, "using", config$sdm_n_folds, "-fold CV..."))
  
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df)) {
    log_error("Invalid inputs for running SDM tuning."); return(NULL)
  }
  
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE) },
                            error = function(e) { log_error(paste("Failed prepare SWD for tuning:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  folds <- tryCatch({ SDMtune::randomFolds(data = full_swd_data, k = config$sdm_n_folds, only_presence = FALSE, seed = 123)},
                    error = function(e){ log_error(paste("Failed create k-folds:", e$message)); return(NULL)})
  if(is.null(folds)) return(NULL)
  
  log_debug("  Training initial CV base model for gridSearch...")
  initial_cv_model <- tryCatch({ SDMtune::train(method = config$sdm_method, data = full_swd_data, folds = folds, progress = FALSE)},
                               error = function(e) { log_error(paste("Failed train initial CV base model:", e$message)); return(NULL)})
  if (is.null(initial_cv_model) || !inherits(initial_cv_model, "SDMmodelCV")) { log_error("Failed create valid SDMmodelCV object."); return(NULL)}
  
  hyper_grid <- config$sdm_tune_grid; tuning_results <- NULL
  tryCatch({
    log_debug("  Running SDMtune::gridSearch with k-fold CV...")
    
    # Determine if progress should be shown safely
    show_progress <- FALSE
    # Only check log level if logger is valid
    if (!is.null(logger)) {
      # Use log4r::log_level to convert string from config to numeric level
      current_log_level_num <- tryCatch(log4r::log_level(config$log_level), error = function(e) log4r::INFO) # Default to INFO on error
      debug_level_num <- tryCatch(log4r::log_level("DEBUG"), error = function(e) log4r::DEBUG) # Get numeric for DEBUG
      show_progress <- current_log_level_num <= debug_level_num
    }
    
    tuning_results <- SDMtune::gridSearch(
      model = initial_cv_model,
      hypers = hyper_grid,
      metric = config$sdm_evaluation_metric,
      save_models = TRUE,
      progress = show_progress, # Use the safe value
      interactive = FALSE
    )
    log_debug("  SDMtune::gridSearch completed.")
    
    res_df <- tuning_results@results
    if(is.null(res_df) || nrow(res_df) == 0) { log_warn("No results found in tuning object."); return(NULL) }
    
    metric_base_upper <- toupper(config$sdm_evaluation_metric)
    target_metric_col <- if (tolower(config$sdm_evaluation_metric) == "aicc") "AICc" else paste0("test_", metric_base_upper)
    
    if (!target_metric_col %in% colnames(res_df)) { log_error(paste("CV metric '", target_metric_col, "' not found. Available:", paste(colnames(res_df), collapse=", "))); print(head(res_df)); return(NULL) }
    
    valid_metric_indices <- which(!is.na(res_df[[target_metric_col]]))
    if(length(valid_metric_indices) == 0) { log_warn(paste("All values for metric '", target_metric_col, "' are NA.")); return(NULL) }
    
    select_fun <- if (target_metric_col == "AICc") which.min else which.max
    best_row_relative_index <- select_fun(res_df[[target_metric_col]][valid_metric_indices])
    best_row_index <- valid_metric_indices[best_row_relative_index]
    
    if(length(best_row_index) == 0 || best_row_index > nrow(res_df)) { log_warn("Could not determine best model index, using first row."); best_row_index <- 1}
    
    best_hypers_df <- res_df[best_row_index, names(hyper_grid), drop = FALSE]
    
    log_info(paste("  Best hypers (Mean CV", target_metric_col,"=", round(res_df[[target_metric_col]][best_row_index], 4),"): ", paste(names(best_hypers_df), best_hypers_df[1,], collapse=", ")))
    
    return(list(best_hypers = best_hypers_df, tuning_results = tuning_results))
    
  }, error = function(e) {
    log_error(paste("SDMtune::gridSearch failed during CV tuning:", e$message))
    return(NULL)
  })
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


#' Predict SDM Suitability
predict_sdm_suitability <- function(final_sdm_model, predictor_stack, config, logger, output_type = "cloglog") {
  # Use logger safely
  log_debug <- function(...) if(!is.null(logger)) log4r::debug(logger, ...)
  log_error <- function(...) if(!is.null(logger)) log4r::error(logger, ...)
  log_warn <- function(...) if(!is.null(logger)) log4r::warn(logger, ...)
  
  log_debug(paste("Predicting suitability (type:", output_type, ")..."))
  
  if (is.null(final_sdm_model) || !inherits(final_sdm_model, "SDMmodel")) { log_error("Invalid final_sdm_model object provided."); return(NULL) }
  if (is.null(predictor_stack)) { log_error("Predictor stack required for prediction."); return(NULL) }
  
  prediction_raster <- NULL
  tryCatch({
    # Determine if progress should be shown safely for predict
    show_progress_predict <- FALSE
    if (!is.null(logger)) {
      current_log_level_num <- tryCatch(log4r::log_level(config$log_level), error = function(e) log4r::INFO)
      debug_level_num <- tryCatch(log4r::log_level("DEBUG"), error = function(e) log4r::DEBUG)
      show_progress_predict <- current_log_level_num <= debug_level_num
    }
    
    prediction_raster <- SDMtune::predict(
      object = final_sdm_model,
      data = predictor_stack,
      type = output_type,
      clamp = TRUE,
      progress = show_progress_predict # Use the safe value
    )
    names(prediction_raster) <- "suitability"
    log_debug("  Prediction raster generated.")
  }, error = function(e) { log_error(paste("SDMtune::predict failed:", e$message)); prediction_raster <- NULL })
  
  return(prediction_raster)
}

#-------------------------------------------------------------------------------
# Plotting helpers (Unchanged from previous version - rely on config for display names)
#-------------------------------------------------------------------------------
#' Plot VIF results using ggplot (Corrected Argument Handling v2)
plot_vif_results_original <- function(vif_result, save_path = NULL, config, display_lookup = NULL) {
  if(!is.numeric(vif_result) || is.null(names(vif_result)) || length(vif_result) == 0){
    warning("plot_vif_results_original expects non-empty named numeric vector.", call.=FALSE); return(NULL)
  }
  
  # Use provided lookup first, then fallback to config lookup
  lookup_to_use <- if (!is.null(display_lookup)) display_lookup else config$core_var_display_names
  get_display_name_func <- if (!is.null(display_lookup)) config$get_display_name else config$get_display_name
  
  if(is.null(lookup_to_use) || is.null(get_display_name_func) || !is.function(get_display_name_func)){
    warning("Display name lookup or function missing/invalid for VIF plot.", call.=FALSE)
    get_display_name_func <- function(name, ...) name # Fallback
    lookup_to_use <- NULL
  }
  
  original_names <- names(vif_result)
  display_names <- sapply(original_names, get_display_name_func, lookup = lookup_to_use, USE.NAMES = FALSE)
  
  df <- data.frame(OriginalName = original_names, DisplayName = display_names, VIF = as.numeric(vif_result))
  num_vars <- nrow(df); if (num_vars == 0) return(NULL)
  df$VIF <- pmax(0, df$VIF)
  df$DisplayName <- factor(df$DisplayName, levels = df$DisplayName[order(df$VIF)])
  
  p <- ggplot2::ggplot(df, aes(x = DisplayName, y = VIF)) +
    ggplot2::geom_bar(stat = "identity", aes(fill = OriginalName), show.legend = FALSE) +
    ggplot2::scale_fill_viridis_d(option = "plasma") +
    ggplot2::labs(title = "VIF Analysis Results (car::vif)", x = "Environmental Variables", y = "VIF Value") +
    ggplot2::scale_y_continuous(limits = c(0, max(config$vif_threshold, ceiling(max(df$VIF, na.rm = TRUE)))), breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::geom_hline(yintercept = config$vif_threshold, linetype = "dashed", color = "red") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (!is.null(save_path)) {
    tryCatch({
      ggplot2::ggsave(save_path, plot = p, width = max(8, 4 + 0.2 * num_vars), height = 6, limitsize = FALSE)
      # Use cat for simple messages outside logger context or if logger is NULL
      cat("  VIF plot saved to", save_path, "\n")
    }, error = function(e) { warning("Failed to save VIF plot: ", e$message, call.=FALSE)})
  }
  return(p)
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