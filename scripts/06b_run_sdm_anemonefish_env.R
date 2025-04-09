# scripts/06b_run_sdm_anemonefish_env.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemonefish Species using ONLY Environmental
# Predictors (PCA or VIF-selected) (SDMtune Workflow) with Parallel, Logging, Progress
# Includes species-specific detailed log files.
#-------------------------------------------------------------------------------
cat("--- Running Script 06b: Run Standard Anemonefish Env SDMs (Parallel & Logged) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr, log4r, future, furrr, progressr)

# Source helpers
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))

# --- 2. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path,
                       log_level = config$log_level,
                       append = config$log_append,
                       log_to_console = config$log_to_console,
                       console_level = config$log_console_level)

log4r::info(logger, "--- Starting Script 06b: Run Standard Anemonefish Env SDMs ---")

# --- 3. Define Group Specifics & Predictor Type ---
group_name <- "anemonefish" # << CHANGED
species_list_file <- config$anemonefish_species_list_file # << CHANGED
occurrence_dir <- config$anemonefish_occurrence_dir   # << CHANGED
use_pca <- config$use_pca_predictors
predictor_type_suffix <- ifelse(use_pca, "_pca", "_vif")
log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---"))

# --- 4. Load Predictor Information ---
predictor_paths_or_list <- NULL
if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (!file.exists(pca_paths_rds)) {
    log4r::fatal(logger, paste("PCA raster paths file not found:", pca_paths_rds))
    stop("PCA paths missing.")
  }
  predictor_paths_or_list <- readRDS(pca_paths_rds)
  if (!is.list(predictor_paths_or_list) || length(predictor_paths_or_list) == 0) {
    log4r::fatal(logger, "PCA raster paths list empty/invalid.")
    stop("PCA paths list invalid.")
  }
  log4r::info(logger, paste("Loaded PCA raster paths for scenarios:", paste(names(predictor_paths_or_list), collapse=", ")))
} else {
  # Load the VIF list specific to anemonefish environment-only models
  predictor_paths_or_list <- config$final_vars_vif_anemonefish_env # << CHANGED
  if(is.null(predictor_paths_or_list) || length(predictor_paths_or_list) < 2) {
    log4r::fatal(logger, "Final VIF var list `final_vars_vif_anemonefish_env` missing/short.") # << CHANGED message
    stop("VIF vars missing.")
  }
  log4r::info(logger, paste("Using VIF-selected core variables:", paste(predictor_paths_or_list, collapse=", ")))
}

# --- 5. Create Output Directories ---
group_pred_dir <- file.path(config$predictions_dir, paste0(group_name, predictor_type_suffix))
group_results_dir <- file.path(config$results_dir, paste0(group_name, predictor_type_suffix))
group_models_dir <- file.path(config$models_dir, paste0(group_name, predictor_type_suffix))
species_log_dir <- config$species_log_dir # Get species log dir from config

dir.create(group_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(species_log_dir, recursive = TRUE, showWarnings = FALSE)

log4r::debug(logger, paste("Output directories created/checked for:", group_name, predictor_type_suffix))
log4r::debug(logger, paste("Species log directory:", species_log_dir))

# --- 6. Load Species List ---
tryCatch({
  species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) {
  log4r::fatal(logger, paste("Failed load species list:", e$message))
  stop("Species list loading failed.")
})
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))

# --- 7. Define Function to Process Single Species (for parallel execution) ---
# This function is identical to the one in 06a, just needs the correct inputs passed below.
# (Keep the function definition exactly as it was in 06a)
process_species_sdm <- function(species_row, config, predictor_paths_or_list, group_pred_dir,
                                group_results_dir, group_models_dir, species_log_dir, # Added species_log_dir
                                predictor_type_suffix, use_pca, occurrence_dir,
                                tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  log_prefix <- paste0("[", species_name, "] ")
  
  species_log_file <- file.path(species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) {
    msg <- paste(Sys.time(), paste0("[",level,"]"), log_prefix, paste0(..., collapse = " "))
    cat(msg, "\n", file = species_log_file, append = TRUE)
  }
  slog("INFO", "--- Starting processing ---")
  
  tuning_results_file <- file.path(group_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, ".rds"))
  final_model_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  slog("DEBUG", "Loading tuning predictors for scenario:", tuning_scenario)
  tuning_predictor_stack <- NULL
  if(use_pca){
    tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
    if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) {
      msg <- paste0(log_prefix, "Skipping: PCA stack for tuning scenario '", tuning_scenario, "' not found.")
      slog("ERROR", msg)
      return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg))
    }
    tuning_predictor_stack <- tryCatch(terra::rast(tuning_predictor_path), error = function(e) {slog("ERROR", "Failed to load tuning PCA stack:", e$message); NULL})
  } else {
    current_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, tuning_scenario, config)
    tuning_predictor_stack <- load_selected_env_data(tuning_scenario, current_vif_vars, config)
  }
  if(is.null(tuning_predictor_stack)) {
    msg <- paste0(log_prefix, "Skipping: Failed load predictor stack for tuning scenario '", tuning_scenario, "'.")
    slog("ERROR", msg)
    return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg))
  }
  slog("DEBUG", "Tuning predictor stack loaded.")
  
  config_for_occ_load <- config
  config_for_occ_load$predictor_stack_for_thinning <- tuning_predictor_stack
  
  slog("DEBUG", "Loading/cleaning/thinning occurrences using tuning stack.")
  occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL)
  if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) {
    msg <- paste0(log_prefix, "Skipping: Insufficient valid/thinned occurrences (found ", occ_data_result$count %||% 0, ", need ", config$min_occurrences_sdm, ").")
    slog("WARN", msg)
    return(list(status = "skipped_occurrences", species = species_name, occurrence_count = occ_data_result$count %||% 0, message = msg))
  }
  occs_coords <- occ_data_result$coords
  occurrence_count_after_thinning <- occ_data_result$count
  slog("INFO", "Occurrence count after clean/thin:", occurrence_count_after_thinning)
  
  if (!config$force_rerun$run_standard_sdms && file.exists(final_model_file)) {
    slog("INFO", "Skipping tuning/training: Final model exists.")
    load_existing_model <- TRUE
    final_model <- tryCatch(readRDS(final_model_file), error=function(e) { slog("ERROR", "Failed to load existing model:", e$message); NULL})
    if(is.null(final_model) || !inherits(final_model, "SDMmodel")){
      msg <- paste0(log_prefix, "Existing final model file is invalid or failed to load. Will attempt re-run.")
      slog("WARN", msg)
      load_existing_model <- FALSE
    }
  } else {
    load_existing_model <- FALSE
    slog("INFO", "Final model not found or rerun forced. Proceeding with tuning/training.")
    
    slog("DEBUG", "Generating background points.")
    background_points <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config, logger = NULL, seed = species_aphia_id)
    if (is.null(background_points)) {
      msg <- paste0(log_prefix, "Skipping: Failed background point generation.")
      slog("ERROR", msg)
      return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("DEBUG", "Background points generated.")
    
    slog("INFO", "Starting hyperparameter tuning.")
    tuning_output <- run_sdm_tuning_kfold(occs_coords, tuning_predictor_stack, background_points, config, logger = NULL, species_name)
    if (is.null(tuning_output) || is.null(tuning_output$best_hypers)) {
      msg <- paste0(log_prefix, "Skipping: Hyperparameter tuning failed.")
      slog("ERROR", msg)
      return(list(status = "error_tuning", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    best_hypers <- tuning_output$best_hypers
    slog("INFO", "Tuning complete. Best hypers:", paste(names(best_hypers), best_hypers[1,], collapse=", "))
    save_tuning_result <- tryCatch({
      saveRDS(tuning_output, tuning_results_file)
      slog("DEBUG", "Tuning results saved to:", basename(tuning_results_file))
      NULL
    }, error = function(e){
      slog("ERROR", "Failed save tuning results:", e$message)
      return(list(status = "error_saving_tuning_results", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0(log_prefix, "Failed save tuning results: ", e$message)))
    })
    if (!is.null(save_tuning_result)) return(save_tuning_result)
    
    slog("INFO", "Starting final model training.")
    final_model <- train_final_sdm(occs_coords, tuning_predictor_stack, background_points, best_hypers, config, logger = NULL, species_name)
    if(is.null(final_model)){
      msg <- paste0(log_prefix, "Skipping: Final model training failed.")
      slog("ERROR", msg)
      return(list(status = "error_training", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("INFO", "Final model training complete.")
    save_model_result <- tryCatch({
      saveRDS(final_model, final_model_file)
      slog("DEBUG", "Final model saved to:", basename(final_model_file))
      NULL
    }, error = function(e){
      slog("ERROR", "Failed save final model:", e$message)
      return(list(status = "error_saving_model", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0(log_prefix, "Failed save final model: ", e$message)))
    })
    if (!is.null(save_model_result)) return(save_model_result)
    
    rm(tuning_predictor_stack, background_points, tuning_output); gc()
  }
  
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- if(use_pca) names(predictor_paths_or_list) else config$env_scenarios
  slog("INFO", paste("Starting predictions for", length(scenarios_to_predict), "scenarios."))
  
  if (!is.null(final_model)) {
    for (pred_scenario in scenarios_to_predict) {
      slog("DEBUG", paste("  Predicting scenario:", pred_scenario))
      pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", pred_scenario, predictor_type_suffix, ".tif"))
      
      if (!config$force_rerun$run_standard_sdms && file.exists(pred_file)) {
        slog("DEBUG", "  Prediction exists. Skipping.")
        next
      }
      
      pred_predictor_stack <- NULL
      if(use_pca) {
        pred_predictor_path <- predictor_paths_or_list[[pred_scenario]]
        if (is.null(pred_predictor_path) || !file.exists(pred_predictor_path)) {
          slog("WARN", paste("  PCA stack not found for scenario:", pred_scenario))
          prediction_errors <- prediction_errors + 1; next
        }
        pred_predictor_stack <- tryCatch(terra::rast(pred_predictor_path), error = function(e) {slog("WARN", paste("  Error loading PCA stack:", e$message)); NULL})
      } else {
        scenario_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, pred_scenario, config)
        if(length(scenario_vif_vars) < 1) {
          slog("WARN", paste("  No VIF variables generated for scenario:", pred_scenario))
          prediction_errors <- prediction_errors + 1; next
        }
        pred_predictor_stack <- load_selected_env_data(pred_scenario, scenario_vif_vars, config)
      }
      if(is.null(pred_predictor_stack)) {
        slog("WARN", paste("  Failed to load predictor stack for scenario:", pred_scenario))
        prediction_errors <- prediction_errors + 1; next
      }
      
      prediction_output <- predict_sdm_suitability(final_model, pred_predictor_stack, config, logger = NULL)
      
      if (inherits(prediction_output, "SpatRaster")) {
        if (!is.null(prediction_output) && terra::nlyr(prediction_output) > 0) {
          save_success <- tryCatch({
            terra::writeRaster(prediction_output, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
            slog("DEBUG", paste("  Prediction raster saved:", basename(pred_file)))
            TRUE
          }, error = function(e) {
            slog("ERROR", paste("  Failed to save prediction raster:", e$message))
            FALSE
          })
          if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
        } else { slog("WARN", paste("  Prediction output was SpatRaster but NULL or empty for scenario:", pred_scenario)); prediction_errors <- prediction_errors + 1 }
        rm(prediction_output); gc()
      } else if (is.character(prediction_output)) { slog("WARN", paste("  Prediction failed for scenario:", pred_scenario, "Reason:", prediction_output)); prediction_errors <- prediction_errors + 1
      } else { slog("ERROR", paste("  Unexpected output type from predict_sdm_suitability for scenario:", pred_scenario, "- Class:", class(prediction_output))); prediction_errors <- prediction_errors + 1 }
      rm(pred_predictor_stack); gc()
    }
  } else {
    prediction_errors <- length(scenarios_to_predict)
    msg <- paste0(log_prefix, "Skipping all predictions: Final model was not available.")
    slog("ERROR", msg)
    status <- if(load_existing_model) "error_loading_model" else "error_training"
    return(list(status = status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
  }
  
  final_status <- "success"
  status_message <- paste0(log_prefix, "Finished. Occurrences (thin/clean): ", occurrence_count_after_thinning, ". Predictions attempted: ", length(scenarios_to_predict), ". Predictions made: ", predictions_made, ". Errors/Skipped preds: ", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) {
    final_status <- "error_prediction_all"
    status_message <- paste0(log_prefix, "Finished with training. Occurrences (thin/clean): ", occurrence_count_after_thinning, ". ALL ", length(scenarios_to_predict), " prediction steps failed/skipped.")
    slog("ERROR", status_message)
  } else if (prediction_errors > 0) {
    final_status <- "success_with_pred_errors"
    status_message <- paste0(log_prefix, "Finished. Occurrences (thin/clean): ", occurrence_count_after_thinning, ". Predictions made: ", predictions_made, ". Errors/Skipped preds: ", prediction_errors, ". Check logs.")
    slog("WARN", status_message)
  } else { slog("INFO", status_message) }
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
}

# --- 8. Setup Parallel Backend & Run ---
if (config$use_parallel && config$num_cores > 1) {
  log4r::info(logger, paste("Setting up parallel backend with", config$num_cores, "cores (multisession)."))
  future::plan(future::multisession, workers = config$num_cores)
} else {
  log4r::info(logger, "Running sequentially.")
  future::plan(future::sequential)
}

# Enable progressr handlers
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

log4r::info(logger, paste("Starting SDM processing for", nrow(species_df), "anemonefish..."))

# Run in parallel with progress
results_list <- progressr::with_progress({
  furrr::future_map(1:nrow(species_df), ~{
    process_species_sdm(
      species_row = species_df[.x, ],
      config = config,
      predictor_paths_or_list = predictor_paths_or_list,
      group_pred_dir = group_pred_dir,
      group_results_dir = group_results_dir,
      group_models_dir = group_models_dir,
      species_log_dir = config$species_log_dir, # Pass the path
      predictor_type_suffix = predictor_type_suffix,
      use_pca = use_pca,
      occurrence_dir = occurrence_dir # Pass correct occurrence dir
    )
  }, .options = furrr_options(seed = TRUE))
})

log4r::info(logger, "Parallel/sequential processing complete.")

# --- 9. Process Results ---
# (Keep the result processing and summary logging exactly as in 06a)
occurrence_counts <- list()
success_count <- 0
error_count <- 0
skipped_count <- 0
partial_success_count <- 0

log4r::info(logger, "--- Processing Results Summary (Anemonefish Env Models) ---")
for (res in results_list) {
  if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "Received NULL result from a species process."); next }
  log_level_func <- switch(res$status,
                           success = log4r::info,
                           success_with_pred_errors = log4r::warn,
                           error_prediction_all = log4r::error,
                           skipped_occurrences = log4r::warn,
                           log4r::error) # Default to error for other error statuses
  log_level_func(logger, res$message)
  occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
  if (res$status == "success") success_count <- success_count + 1
  else if (res$status == "success_with_pred_errors") partial_success_count <- partial_success_count + 1
  else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
  else error_count <- error_count + 1
}

log4r::info(logger, paste("--- Overall Summary (Anemonefish Env Models) ---"))
log4r::info(logger, paste("Total Species Targets:", length(results_list)))
log4r::info(logger, paste("Fully Successful Runs (Train + All Preds):", success_count))
log4r::info(logger, paste("Partially Successful Runs (Train + Some Pred Errors):", partial_success_count))
log4r::info(logger, paste("Skipped (Occurrences/Thinning/etc.):", skipped_count))
log4r::info(logger, paste("Errors during processing (Train/Tune/Setup/Save):", error_count))

if (length(occurrence_counts) > 0) {
  occ_count_df <- data.frame(
    Species = names(occurrence_counts),
    OccurrenceCountAfterThinning = unlist(occurrence_counts)
  ) %>% dplyr::arrange(Species)
  occ_count_df_print <- occ_count_df
  occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
  occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, ".csv"))
  tryCatch({ readr::write_csv(occ_count_df, occ_count_file); log4r::info(logger, paste("Occurrence counts saved to:", occ_count_file)) }, error = function(e) { log4r::error(logger, paste("Failed to save occurrence counts:", e$message)) })
  cat("\n--- Final Occurrence Counts (Anemonefish, after cleaning/thinning) ---\n")
  print(occ_count_df_print)
} else { log4r::warn(logger, "No occurrence counts were recorded.") }

# Shut down parallel workers
future::plan(future::sequential)

log4r::info(logger, "--- Script 06b finished. ---")
#-------------------------------------------------------------------------------