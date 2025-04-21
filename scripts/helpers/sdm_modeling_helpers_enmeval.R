# scripts/helpers/sdm_modeling_helpers_enmeval.R
#-------------------------------------------------------------------------------
# Helper functions specifically for the ENMeval SDM Workflow (v1.1 - Internal Partitions)
#-------------------------------------------------------------------------------
pacman::p_load(ENMeval, sf, terra, dplyr, readr, tools, stringr, log4r, predicts)

#' Select Optimal Model from ENMevaluation Results
#'
#' Selects the best model tuning settings based on criteria specified in the config.
#'
#' @param enmeval_results An ENMevaluation object.
#' @param config Project configuration list.
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return The tune.args name (character string) of the selected best model, or NULL.
select_enmeval_model <- function(enmeval_results, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[ModelSelectHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  if (!inherits(enmeval_results, "ENMevaluation")) { hlog("ERROR", "Invalid ENMevaluation object."); return(NULL) }
  
  results_df <- enmeval_results@results
  if (is.null(results_df) || nrow(results_df) == 0) { hlog("ERROR", "Results table is empty."); return(NULL) }
  
  selection_metric <- config$enmeval_selection_metric
  hlog("INFO", paste("Selecting model based on criteria:", selection_metric))
  
  # Optionally remove rows with NA results for key metrics before selection
  if (config$enmeval_omit_na_models_from_selection %||% TRUE) { # Default to TRUE if missing
    key_metrics <- c("AICc", "auc.val.avg", "cbi.val.avg", "or.mtp.avg", "or.10p.avg")
    cols_to_check <- intersect(key_metrics, names(results_df))
    if (length(cols_to_check) > 0) {
      rows_before <- nrow(results_df)
      results_df <- results_df[stats::complete.cases(results_df[, cols_to_check, drop = FALSE]), ]
      rows_after <- nrow(results_df)
      if (rows_after < rows_before) {
        hlog("WARN", paste("Removed", rows_before - rows_after, "models with NA evaluation metrics before selection."))
      }
    }
    if (nrow(results_df) == 0) { hlog("ERROR", "No models remaining after removing NAs."); return(NULL) }
  }
  
  best_model_tune_args <- NULL
  
  tryCatch({
    if (selection_metric == "AICc") {
      if (!"AICc" %in% names(results_df)) stop("AICc column not found in results.")
      min_aicc <- min(results_df$AICc, na.rm = TRUE)
      if (!is.finite(min_aicc)) stop("Minimum AICc is not finite.")
      
      # Select models within the delta AICc threshold
      selected_models <- results_df[which(results_df$AICc <= min_aicc + config$enmeval_delta_aicc_threshold), ]
      
      # Among those, choose the one with the fewest coefficients (simplest)
      if (nrow(selected_models) > 0) {
        if ("ncoef" %in% names(selected_models)) {
          best_row_index <- which.min(selected_models$ncoef)
        } else {
          hlog("WARN", "ncoef column not found for AICc tie-breaking, selecting first model within threshold.")
          best_row_index <- 1
        }
        best_model_tune_args <- selected_models$tune.args[best_row_index]
        hlog("INFO", paste("Selected model:", best_model_tune_args, "(Lowest AICc within threshold, fewest coefficients if available)"))
      }
    } else if (selection_metric == "AUC_MTPO") { # Example sequential criteria
      if (!"or.mtp.avg" %in% names(results_df)) stop("or.mtp.avg column not found.")
      if (!"auc.val.avg" %in% names(results_df)) stop("auc.val.avg column not found.")
      
      # Ensure metrics are numeric before filtering
      results_df$or.mtp.avg <- as.numeric(results_df$or.mtp.avg)
      results_df$auc.val.avg <- as.numeric(results_df$auc.val.avg)
      
      min_omission <- min(results_df$or.mtp.avg, na.rm = TRUE)
      if (!is.finite(min_omission)) stop("Minimum MTP Omission Rate is not finite.")
      
      selected_models <- results_df %>%
        dplyr::filter(or.mtp.avg <= min_omission + 1e-6) # Add tolerance for float comparison
      
      if (nrow(selected_models) > 0) {
        max_auc <- max(selected_models$auc.val.avg, na.rm = TRUE)
        if (!is.finite(max_auc)) stop("Maximum AUC is not finite among models with lowest MTP Omission.")
        selected_models <- selected_models %>%
          dplyr::filter(auc.val.avg >= max_auc - 1e-6) # Add tolerance
      }
      
      if (nrow(selected_models) == 1) {
        best_model_tune_args <- selected_models$tune.args[1]
        hlog("INFO", paste("Selected model:", best_model_tune_args, "(Lowest MTP Omission, highest AUC)"))
      } else if (nrow(selected_models) > 1) {
        # Tie-breaking: choose simplest (fewest coefficients) if available
        if ("ncoef" %in% names(selected_models)) {
          best_row_index <- which.min(selected_models$ncoef)
          best_model_tune_args <- selected_models$tune.args[best_row_index]
          hlog("INFO", paste("Selected model:", best_model_tune_args, "(Lowest MTP Omission, highest AUC, fewest coefficients tiebreak)"))
        } else {
          hlog("WARN", "ncoef column not found for tie-breaking, selecting first model.")
          best_model_tune_args <- selected_models$tune.args[1]
        }
      }
    } else if (selection_metric %in% names(results_df)) {
      # Simple selection based on a single metric (assuming higher is better unless AICc)
      metric_values <- as.numeric(results_df[[selection_metric]])
      if (all(is.na(metric_values))) stop(paste("All values for metric", selection_metric, "are NA."))
      
      select_func <- if (selection_metric == "AICc") which.min else which.max
      best_row_index <- select_func(metric_values)
      
      if (length(best_row_index) == 1) {
        best_model_tune_args <- results_df$tune.args[best_row_index]
        hlog("INFO", paste("Selected model:", best_model_tune_args, "(Optimized for", selection_metric, ")"))
      } else if (length(best_row_index) > 1) {
        # Tie-breaking if multiple models have the exact same best metric value
        hlog("WARN", paste("Multiple models tied for best", selection_metric, ". Selecting first one:", results_df$tune.args[best_row_index[1]]))
        best_model_tune_args <- results_df$tune.args[best_row_index[1]]
      }
    } else {
      stop(paste("Unsupported or missing selection metric:", selection_metric))
    }
    
    if (is.null(best_model_tune_args)) {
      hlog("WARN", "Could not determine a single best model based on criteria. Returning NULL.")
    }
    
    return(best_model_tune_args)
    
  }, error = function(e) {
    hlog("ERROR", paste("Failed during model selection:", e$message))
    return(NULL)
  })
}


#' Save ENMeval Results (Evaluation Table CSV and Full Object RDS)
#'
#' Saves the results table to the target directory and the full object
#' to the intermediate directory.
#'
#' @param enmeval_results The ENMevaluation object.
#' @param species_name_sanitized Sanitized species name.
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_enmeval_pca"). Ensure this is specific for ENMeval runs.
#' @param config Project configuration list.
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure of saving the CSV.
save_enmeval_results <- function(enmeval_results, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveENMResHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  if (!inherits(enmeval_results, "ENMevaluation")) { hlog("ERROR", "Invalid ENMevaluation object provided."); return(FALSE) }
  
  # --- Save Full ENMevaluation Object (Intermediate) ---
  # Ensure group_name is accessible within config or passed explicitly
  group_name <- config$group_name %||% "unknown_group"
  intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix))
  dir.create(intermediate_results_dir, recursive = TRUE, showWarnings = FALSE)
  rds_filename <- paste0("enmevaluation_object_", species_name_sanitized, predictor_type_suffix, ".rds")
  rds_file_path <- file.path(intermediate_results_dir, rds_filename)
  tryCatch({
    saveRDS(enmeval_results, file = rds_file_path)
    hlog("DEBUG", "Full ENMevaluation object saved (RDS):", basename(rds_file_path))
  }, error = function(e) {
    hlog("WARN", paste("Failed to save full ENMevaluation RDS:", e$message))
  })
  
  # --- Save Results Table (Target) ---
  results_df <- enmeval_results@results
  if (is.null(results_df) || nrow(results_df) == 0) {
    hlog("WARN", "No results table found in ENMevaluation object to save as CSV.")
    return(TRUE) # Return TRUE as the main object might have saved
  }
  
  # Construct target path using helper
  csv_file_path <- construct_results_filename(species_name_sanitized, predictor_type_suffix, config)
  if(is.null(csv_file_path)){
    hlog("ERROR", "Failed to construct target path for results CSV.")
    return(FALSE)
  }
  
  # Ensure target directory exists
  target_subdir <- dirname(csv_file_path)
  dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    readr::write_csv(results_df, file = csv_file_path)
    hlog("DEBUG", "ENMeval results table saved (CSV):", basename(csv_file_path))
    TRUE
  }, error = function(e) {
    hlog("ERROR", paste("Failed save ENMeval results CSV:", e$message))
    FALSE
  })
}

#' Calculate and Save Variable Importance from ENMeval Model (MaxNet specific)
#'
#' Extracts pre-calculated permutation importance from a MaxNet model object
#' obtained via ENMeval. Saves results to the target structure.
#'
#' @param optimal_model An `SDMmodel` object (the selected one from ENMeval).
#' @param species_name_sanitized Sanitized species name.
#' @param group_name Group name ("anemone" or "anemonefish").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_enmeval_pca"). Ensure this is specific for ENMeval runs.
#' @param config The configuration list.
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
calculate_and_save_vi_enmeval <- function(optimal_model, species_name_sanitized, group_name, predictor_type_suffix,
                                          config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[VarImpENMHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  if (is.null(optimal_model) || !inherits(optimal_model, "SDMmodel")) {
    hlog("ERROR", "Invalid optimal_model object provided for VI.")
    return(FALSE)
  }
  
  hlog("INFO", "Extracting variable importance...")
  vi_results_df <- NULL
  
  tryCatch({
    model_method <- optimal_model@algorithm # ENMeval stores algorithm name here
    
    if (model_method == "maxnet") {
      hlog("DEBUG", "Maxnet model detected. Extracting pre-calculated permutation importance.")
      raw_maxnet_model <- optimal_model@model # Access the underlying maxnet object
      if (!is.null(raw_maxnet_model) && !is.null(raw_maxnet_model$variable.importance)) {
        importance_vector <- raw_maxnet_model$variable.importance
        vi_results_df <- data.frame(
          Variable = names(importance_vector),
          Importance = as.numeric(importance_vector),
          stringsAsFactors = FALSE
        )
        vi_results_df <- vi_results_df[order(-vi_results_df$Importance), ] # Order descending
        hlog("DEBUG", "Successfully extracted MaxNet variable importance.")
      } else {
        hlog("WARN", "Could not find pre-calculated VI in the MaxNet model object.")
        # NOTE: ENMeval doesn't easily support recalculating VI without the original SWD data easily accessible
        # We will rely on the built-in importance. If this fails, VI cannot be calculated here.
      }
    } else {
      # VI extraction for other ENMeval algorithms might differ
      hlog("WARN", paste("Variable importance extraction not explicitly implemented for algorithm:", model_method, "in this helper. Check model object structure."))
      # Attempt to find a standard 'variable.importance' slot if it exists
      if ("variable.importance" %in% slotNames(optimal_model@model)) {
        importance_data <- slot(optimal_model@model, "variable.importance")
        if(is.numeric(importance_data) && !is.null(names(importance_data))){
          vi_results_df <- data.frame(Variable = names(importance_data), Importance = as.numeric(importance_data))
          vi_results_df <- vi_results_df[order(-vi_results_df$Importance), ]
          hlog("DEBUG", "Extracted VI from 'variable.importance' slot.")
        } else if (is.data.frame(importance_data) && all(c("Variable", "Importance") %in% names(importance_data))){
          vi_results_df <- importance_data[, c("Variable", "Importance")]
          vi_results_df <- vi_results_df[order(-vi_results_df$Importance), ]
          hlog("DEBUG", "Extracted VI from 'variable.importance' data frame.")
        } else {
          hlog("WARN", "Found 'variable.importance' slot, but format not recognized.")
        }
      } else {
        hlog("WARN", "No standard 'variable.importance' found for this algorithm.")
      }
    }
    
    # --- Saving Logic ---
    if (is.null(vi_results_df)) {
      hlog("ERROR", "Variable importance could not be obtained.")
      return(FALSE)
    }
    
    # Construct target path using helper
    vi_file_path <- construct_vi_filename(species_name_sanitized, group_name, predictor_type_suffix, config)
    if(is.null(vi_file_path)){
      hlog("ERROR", "Failed to construct target path for VI CSV.")
      return(FALSE)
    }
    
    target_subdir <- dirname(vi_file_path)
    dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
    
    readr::write_csv(vi_results_df, vi_file_path)
    hlog("INFO", paste("Variable importance saved to:", vi_file_path))
    return(TRUE)
    
  }, error = function(e) {
    hlog("ERROR", paste("Variable importance processing/saving failed:", e$message))
    return(FALSE)
  })
}


#' Predict SDM Suitability using Optimal Model from ENMeval
#' Uses terra::predict.
#' @param optimal_sdm_model The selected optimal `SDMmodel` object from ENMeval.
#' @return A SpatRaster object with the prediction, or a character string with the error message on failure.
predict_sdm_suitability <- function(optimal_sdm_model, predictor_stack, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[PredictHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  output_type <- config$enmeval_pred_type %||% "cloglog"
  do_clamp <- config$enmeval_clamp %||% TRUE
  
  hlog("DEBUG", paste("Attempting prediction (type:", output_type, ", clamp:", do_clamp, ")..."))
  
  if (is.null(optimal_sdm_model) || !inherits(optimal_sdm_model, "SDMmodel")) { msg <- "Invalid optimal_sdm_model."; hlog("ERROR", msg); return(msg) }
  if (is.null(predictor_stack)) { msg <- "Predictor stack required."; hlog("ERROR", msg); return(msg) }
  
  prediction_result <- NULL
  tryCatch({
    # Use terra::predict with the internal model object
    # We need the actual prediction function associated with the algorithm
    # For maxnet, we can use predicts::predict
    model_algorithm <- optimal_sdm_model@algorithm
    
    predict_func <- NULL
    if(model_algorithm == "maxnet"){
      predict_func <- predicts::predict # Function from predicts package for maxnet/glmnet
    } else if (model_algorithm == "bioclim") {
      # Need to confirm predicts::predict works for bioclim or find the right func
      predict_func <- predicts::predict # Assuming it works, may need adjustment
      # BIOCLIM might need specific arguments like 'tails'
      # predict_args <- list(tails = ...)
    } else {
      msg <- paste("Prediction function mapping not defined for algorithm:", model_algorithm); hlog("ERROR", msg); return(msg)
    }
    
    hlog("DEBUG", paste("Using prediction function for", model_algorithm))
    
    prediction_result <- terra::predict(
      object = predictor_stack,       # Raster data
      model = optimal_sdm_model@model, # The internal model object (e.g., maxnet object)
      fun = predict_func,             # The prediction function
      type = output_type,             # Passed to the predict_func
      clamp = do_clamp,               # Passed to the predict_func
      na.rm = TRUE
    )
    
    if (is.null(prediction_result) || !inherits(prediction_result, "SpatRaster") || terra::nlyr(prediction_result) == 0) {
      msg <- "Prediction returned NULL or empty raster."; hlog("WARN", msg); return(msg)
    }
    
    min_pred <- global(prediction_result, "min", na.rm=TRUE)$min
    max_pred <- global(prediction_result, "max", na.rm=TRUE)$max
    if (!is.na(min_pred) && !is.na(max_pred) && abs(min_pred - max_pred) < 1e-9) { # Check for near-constant values
      hlog("WARN", paste("Prediction resulted in a near-constant value:", round(min_pred, 4)))
    }
    
    names(prediction_result) <- "suitability"
    hlog("DEBUG", "Prediction raster generated.")
    return(prediction_result)
    
  }, error = function(e) {
    err_msg <- paste("terra::predict failed:", e$message)
    hlog("ERROR", err_msg)
    hlog("ERROR", " Predictor stack names: ", paste(names(predictor_stack), collapse=", "))
    # Depending on the model object structure, add more debug info if possible
    # e.g., if optimal_sdm_model@model$betas exists for maxnet:
    # if(optimal_sdm_model@algorithm == "maxnet" && !is.null(optimal_sdm_model@model$betas)) {
    #    hlog("ERROR", " Model coefficient names: ", paste(names(optimal_sdm_model@model$betas), collapse=", "))
    # }
    return(err_msg)
  })
}

# --- Include necessary helper functions from the original file ---
# Make sure these functions are defined either here or loaded via source()
# construct_results_filename(...)
# construct_vi_filename(...)
# construct_prediction_filename(...)  # From original helpers
# save_final_model(...)            # From original helpers (saves INTERMEDIATE model)
# save_sdm_prediction(...)         # From original helpers (saves TARGET prediction)

# Define here if not sourced:
construct_results_filename <- function(species_name_sanitized, predictor_type_suffix, config) {
  subdir_name <- config$model_output_subdir_map[[predictor_type_suffix]]
  if (is.null(subdir_name)) {
    warning("No output subdirectory mapping found for suffix: ", predictor_type_suffix, call. = FALSE)
    return(NULL)
  }
  target_subdir <- file.path(config$target_results_base, subdir_name)
  target_base_name <- paste0("CV_Results_", species_name_sanitized)
  target_filename <- paste0(target_base_name, ".csv")
  return(file.path(target_subdir, target_filename))
}

construct_vi_filename <- function(species_name_sanitized, group_name, predictor_type_suffix, config) {
  vi_subdir_name <- paste0("vi_", group_name)
  target_subdir <- file.path(config$target_vi_base, vi_subdir_name)
  vi_filename <- paste0("vi_", species_name_sanitized, predictor_type_suffix, ".csv")
  return(file.path(target_subdir, vi_filename))
}

# --- Make sure these are loaded/defined before use in the main script ---
# source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))
# --- Or copy the functions directly: ---
construct_prediction_filename <- function(species_name_sanitized, scenario_name, predictor_type_suffix, config) {
  base_filename <- paste0("mean_pred_", species_name_sanitized)
  target_dir <- NULL; target_filename_stem <- NULL
  if (scenario_name == "current") {
    target_dir <- config$target_predictions_current_dir
    target_filename_stem <- paste0(base_filename, predictor_type_suffix)
  } else {
    ssp_dir_name <- config$ssp_scenario_map[[scenario_name]]
    if (is.null(ssp_dir_name)) { warning("No SSP directory mapping for scenario: ", scenario_name, call. = FALSE); return(NULL) }
    target_dir <- file.path(config$target_predictions_future_base, ssp_dir_name)
    time_tag_clean <- ifelse(grepl("2050$", scenario_name), "dec50", ifelse(grepl("2100$", scenario_name), "dec100", "UNKNOWN"))
    if(time_tag_clean == "UNKNOWN") { warning("Cannot extract time tag from scenario: ", scenario_name, call. = FALSE); return(NULL) }
    target_filename_stem <- paste0(base_filename, "_", ssp_dir_name, "_", time_tag_clean, predictor_type_suffix)
  }
  if(is.null(target_dir) || is.null(target_filename_stem)) { warning("Could not construct target path components.", call. = FALSE); return(NULL) }
  target_filename <- paste0(target_filename_stem, ".tif")
  return(file.path(target_dir, target_filename))
}

save_final_model <- function(final_model, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveModelHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Invalid model object provided."); return(FALSE) }
  group_name <- config$group_name %||% "unknown_group" # Get group name
  model_subdir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
  dir.create(model_subdir, recursive = TRUE, showWarnings = FALSE)
  model_file <- file.path(model_subdir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  tryCatch({ saveRDS(final_model, file = model_file); hlog("DEBUG", "Optimal model object saved (RDS - Intermediate):", basename(model_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save optimal model RDS:", e$message)); FALSE })
}

save_sdm_prediction <- function(prediction_raster, species_name_sanitized, scenario_name, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SavePredHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  if (is.null(prediction_raster) || !inherits(prediction_raster, "SpatRaster")) { hlog("ERROR", "Invalid prediction raster provided."); return(FALSE) }
  pred_file_path <- construct_prediction_filename(species_name_sanitized, scenario_name, predictor_type_suffix, config)
  if (is.null(pred_file_path)) { hlog("ERROR", "Failed to construct prediction filename."); return(FALSE) }
  target_dir <- dirname(pred_file_path)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    terra::writeRaster(prediction_raster, filename = pred_file_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    hlog("DEBUG", paste("Prediction raster saved (Target):", basename(pred_file_path), "to", target_dir))
    TRUE
  }, error = function(e) { hlog("ERROR", paste("Failed save prediction raster:", e$message)); FALSE })
}