# scripts/ENMeval/06a_run_sdm_anemone_enmeval.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using ENMeval for tuning/evaluation.
# Uses Environmental Predictors (PCA or VIF-selected).
# Includes Parallel Processing, Logging, Progress, and Species-Specific Logs.
# Saves outputs to target structure defined in config.R.
# Based on ENMeval vignette methodology (v1.1 - Internal Partitions).
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (ENMeval Workflow v1.1) ---\n")

# --- 1. Setup: Load Config FIRST ---
# Assuming config.R is in the parent 'scripts' directory relative to base_dir
if (file.exists("../config.R")) { # Relative path from within ENMeval dir
  source("../config.R")
  if (!exists("config") || !is.list(config)) {
    stop("FATAL: 'config' list object not found or invalid after sourcing config.R")
  }
} else {
  # Fallback if run from base_dir
  if (file.exists("scripts/config.R")) {
    source("scripts/config.R")
    if (!exists("config") || !is.list(config)) {
      stop("FATAL: 'config' list object not found or invalid after sourcing config.R")
    }
  } else {
    stop("FATAL: Configuration file 'scripts/config.R' not found.")
  }
}


# --- 2. Load Required Packages ---
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
# Ensure all necessary packages are listed in config/01_install
pacman::p_load(terra, sf, dplyr, readr, tools, stringr, log4r, future, furrr, progressr, ENMeval, predicts)

# --- 3. Source Helper Functions ---
# Paths relative to the base_dir defined in config
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R")) # Contains occ loading, bg gen, etc.
source(file.path(config$helpers_dir, "sdm_modeling_helpers_enmeval.R")) # Contains ENMeval specific logic

# --- 4. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path,
                       log_level = config$log_level,
                       append = config$log_append,
                       log_to_console = config$log_to_console,
                       console_level = config$log_console_level)

log4r::info(logger, "--- Starting Script 06a: Run Standard Anemone SDMs (ENMeval Workflow v1.1) ---")

# --- 5. Define Group Specifics & Predictor Type ---
group_name <- "anemone"
config$group_name <- group_name # Add to config dynamically for helpers
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors
# Adjust suffix to clearly indicate ENMeval run
predictor_type_suffix <- ifelse(use_pca, "_enmeval_pca", "_enmeval_vif")

log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---"))
log4r::info(logger, paste("--- ENMeval Algorithm(s):", paste(config$enmeval_algorithms, collapse=", "), "---"))
log4r::info(logger, paste("--- ENMeval Partition Method:", config$enmeval_partitions, "---"))
log4r::info(logger, paste("--- ENMeval Model Selection Metric:", config$enmeval_selection_metric, "---"))

# --- 6. Load Predictor Information ---
predictor_paths_or_list <- NULL
if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (is.null(pca_paths_rds) || !file.exists(pca_paths_rds)) { log4r::fatal(logger, "PCA paths RDS missing/invalid."); stop("PCA paths RDS invalid.") }
  predictor_paths_or_list <- tryCatch(readRDS(pca_paths_rds), error = function(e) { log4r::fatal(logger, paste("Failed load PCA paths RDS:", e$message)); NULL })
  if (!is.list(predictor_paths_or_list) || length(predictor_paths_or_list) == 0) { log4r::fatal(logger, "PCA paths list empty/invalid."); stop("PCA paths list invalid.") }
  log4r::info(logger, paste("Loaded PCA raster paths for scenarios:", paste(names(predictor_paths_or_list), collapse=", ")))
} else {
  predictor_paths_or_list <- config$final_vars_vif_anemone # Use VIF list defined for anemone
  if(is.null(predictor_paths_or_list) || !is.character(predictor_paths_or_list) || length(predictor_paths_or_list) < 2) { log4r::fatal(logger, "Final VIF var list `final_vars_vif_anemone` invalid."); stop("VIF vars missing.") }
  log4r::info(logger, paste("Using VIF-selected core variables:", paste(predictor_paths_or_list, collapse=", ")))
}

# --- 7. Create Intermediate Output Dirs ---
base_intermediate_model_path <- config$models_dir_intermediate
base_intermediate_results_path <- config$results_dir_intermediate
base_species_log_path <- config$species_log_dir
intermediate_models_dir <- file.path(base_intermediate_model_path, paste0(group_name, predictor_type_suffix))
intermediate_results_dir <- file.path(base_intermediate_results_path, paste0(group_name, predictor_type_suffix))
tryCatch({
  dir.create(intermediate_models_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(intermediate_results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_species_log_path, recursive = TRUE, showWarnings = FALSE)
  log4r::debug(logger, paste("Intermediate dirs created/checked:", group_name, predictor_type_suffix))
}, error = function(e) { log4r::fatal(logger, paste("Failed create intermediate dirs:", e$message)); stop("Directory creation failed.") })

# --- 8. Load Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE) }, error = function(e) { log4r::fatal(logger, paste("Failed load species list:", e$message)); stop("Species list failed.") })
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))

# --- 9. Define Function to Process Single Species using ENMeval ---
process_species_sdm_enmeval <- function(species_row, config, predictor_paths_or_list, group_name, predictor_type_suffix, use_pca, occurrence_dir,
                                        tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  species_log_file <- file.path(config$species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), paste0("[", species_name, "]"), paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
  slog("INFO", paste0("--- Starting ENMeval processing (", predictor_type_suffix, ") ---"))
  
  # --- Define File Paths ---
  intermediate_models_dir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
  intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix))
  optimal_model_file <- file.path(intermediate_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  enmeval_rds_file <- file.path(intermediate_results_dir, paste0("enmevaluation_object_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  # --- Load Predictors for Tuning Scenario ---
  slog("DEBUG", "Loading tuning predictors for scenario:", tuning_scenario)
  tuning_predictor_stack <- NULL
  if(use_pca){
    tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
    if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) { msg <- paste0("Skipping: PCA stack path for tuning scenario '", tuning_scenario, "' not found."); slog("ERROR", msg); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
    tuning_predictor_stack <- tryCatch(terra::rast(tuning_predictor_path), error = function(e) {slog("ERROR", "Failed to load tuning PCA stack:", e$message); NULL})
  } else {
    current_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, tuning_scenario, config)
    tuning_predictor_stack <- load_selected_env_data(tuning_scenario, current_vif_vars, config)
  }
  if(is.null(tuning_predictor_stack)) { msg <- paste0("Skipping: Failed load predictor stack for tuning scenario '", tuning_scenario, "'."); slog("ERROR", msg); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
  if(terra::crs(tuning_predictor_stack) == "") { msg <- paste0("Skipping: Tuning predictor stack has no CRS assigned."); slog("ERROR", msg); return(list(status = "error_tuning_predictors_crs", species = species_name, occurrence_count = NA, message = msg)) }
  slog("DEBUG", "Tuning predictor stack loaded. Names:", paste(names(tuning_predictor_stack), collapse=", "))
  
  # --- Load, Clean, and Thin Occurrences ---
  config_for_occ_load <- config; config_for_occ_load$predictor_stack_for_thinning <- tuning_predictor_stack
  slog("DEBUG", "Loading/cleaning/thinning occurrences.")
  occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger=NULL, species_log_file=species_log_file)
  if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) { msg <- paste0("Skipping: Insufficient valid/thinned occurrences (found ", occ_data_result$count %||% 0, "). Required: ", config$min_occurrences_sdm); slog("WARN", msg); return(list(status = "skipped_occurrences", species = species_name, occurrence_count = occ_data_result$count %||% 0, message = msg)) }
  occs_coords_df <- as.data.frame(occ_data_result$coords); colnames(occs_coords_df) <- c("longitude", "latitude")
  occurrence_count_after_thinning <- occ_data_result$count
  slog("INFO", "Occurrence count after clean/thin:", occurrence_count_after_thinning)
  
  # --- Generate Background Points ---
  slog("DEBUG", "Generating background points.")
  background_points_df <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config, logger=NULL, species_log_file=species_log_file, seed = species_aphia_id)
  if (is.null(background_points_df)) { msg <- paste0("Skipping: Failed background point generation."); slog("ERROR", msg); return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg)) }
  colnames(background_points_df) <- c("longitude", "latitude")
  slog("DEBUG", "Background points generated:", nrow(background_points_df))
  
  # --- Check if Final Model Already Exists ---
  optimal_model <- NULL; enmeval_results <- NULL; load_existing_model <- FALSE
  if (!config$force_rerun$run_enmeval_sdms && file.exists(optimal_model_file)) { # Check ENMeval specific flag
    slog("INFO", "Optimal model file exists. Attempting to load...")
    optimal_model <- tryCatch(readRDS(optimal_model_file), error=function(e) { slog("ERROR", "Failed to load existing optimal model:", e$message); NULL})
    if(!is.null(optimal_model) && inherits(optimal_model, "SDMmodel")){ # Check class? ENMeval models are nested
      # Check if the loaded model is compatible (e.g., algorithm matches config) - needs more robust check
      if(optimal_model@algorithm %in% config$enmeval_algorithms){
        slog("INFO", "Successfully loaded existing optimal model. Skipping tuning.")
        load_existing_model <- TRUE
        if(file.exists(enmeval_rds_file)) { enmeval_results <- tryCatch(readRDS(enmeval_rds_file), error=function(e){slog("WARN","Could not load full ENMevaluation RDS."); NULL}) }
      } else {
        msg <- paste0("Existing model algorithm (", optimal_model@algorithm, ") does not match config. Will re-run."); slog("WARN", msg)
        optimal_model <- NULL; load_existing_model <- FALSE
      }
    } else {
      msg <- paste0("Existing optimal model file invalid/failed load. Will re-run tuning."); slog("WARN", msg)
      optimal_model <- NULL; load_existing_model <- FALSE
    }
  }
  
  # --- Run ENMeval Tuning if needed ---
  if (!load_existing_model) {
    slog("INFO", "Optimal model not found/invalid or rerun forced. Proceeding with ENMeval tuning.")
    
    # Prepare arguments for ENMevaluate
    enmeval_args <- list(
      occs = occs_coords_df,
      envs = tuning_predictor_stack,
      bg = background_points_df,
      tune.args = config$enmeval_tuning_settings,
      algorithm = config$enmeval_algorithms, # Vector of algorithms
      partitions = config$enmeval_partitions,
      other.settings = list(pred.type = config$enmeval_pred_type, abs.auc.diff = FALSE, validation.bg = "partition"),
      partition.settings = config$enmeval_partition_settings, # List of settings like ag. factor
      clamp = config$enmeval_clamp,
      parallel = config$use_parallel && config$enmeval_num_cores > 1,
      numCores = if(config$use_parallel && config$enmeval_num_cores > 1) config$enmeval_num_cores else 1,
      quiet = TRUE # Set to FALSE for more verbose output
    )
    
    # Run ENMevaluate
    slog("INFO", paste("Running ENMevaluate for algorithm(s):", paste(enmeval_args$algorithm, collapse=", ")))
    enmeval_results <- tryCatch({
      do.call(ENMeval::ENMevaluate, enmeval_args)
    }, error = function(e) {
      slog("ERROR", paste("ENMevaluate failed:", e$message))
      NULL
    })
    
    if (is.null(enmeval_results) || !inherits(enmeval_results, "ENMevaluation")) {
      msg <- paste0("Skipping: ENMevaluate did not return a valid result object."); slog("ERROR", msg)
      return(list(status = "error_enmeval_run", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("INFO", "ENMevaluate completed.")
    
    # Save ENMeval Results (CSV to target, RDS to intermediate)
    if (!save_enmeval_results(enmeval_results, species_name_sanitized, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)) {
      slog("WARN", "Failed to save ENMeval results table (CSV) to target directory.")
    }
    
    # Select the best model (assuming only one algorithm run for now)
    # If multiple algorithms were run, selection logic needs adjustment
    if(length(config$enmeval_algorithms) > 1) {
      slog("WARN", "Multiple algorithms run, selection logic currently selects best overall based on metric. Consider algorithm-specific selection.")
      # Might need to filter results_df by algorithm before calling select_enmeval_model
    }
    selected_model_name <- select_enmeval_model(enmeval_results, config, logger=NULL, species_log_file=species_log_file)
    if (is.null(selected_model_name)) {
      msg <- paste0("Skipping: Failed to select an optimal model from ENMeval results."); slog("ERROR", msg)
      return(list(status = "error_model_selection", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    
    # Get the optimal model object
    optimal_model <- enmeval_results@models[[selected_model_name]]
    if (is.null(optimal_model)) {
      msg <- paste0("Skipping: Selected optimal model '", selected_model_name,"' is NULL."); slog("ERROR", msg)
      return(list(status = "error_optimal_model_null", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("INFO", paste("Optimal model selected:", selected_model_name))
    
    # Save the selected optimal model object (intermediate)
    if (!save_final_model(optimal_model, species_name_sanitized, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)) {
      slog("WARN", "Proceeding without saved intermediate model file, but predictions might fail if run separately.")
    }
  } # End if !load_existing_model
  
  # --- Variable Importance ---
  if (!is.null(optimal_model)) {
    slog("INFO", "Calculating/Extracting and saving variable importance...")
    vi_success <- calculate_and_save_vi_enmeval( # Use the ENMeval specific VI helper
      optimal_model = optimal_model,
      species_name_sanitized = species_name_sanitized,
      group_name = group_name,
      predictor_type_suffix = predictor_type_suffix,
      config = config,
      logger = NULL,
      species_log_file = species_log_file)
    if(!vi_success) slog("WARN", "Variable importance calculation/saving failed.")
  } else {
    slog("WARN", "Skipping variable importance: optimal model object is unavailable.")
  }
  
  # --- Prediction Loop ---
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- if(use_pca) names(predictor_paths_or_list) else config$env_scenarios
  slog("INFO", paste("Starting predictions for", length(scenarios_to_predict), "scenarios using optimal model."))
  
  if (!is.null(optimal_model)) {
    for (pred_scenario in scenarios_to_predict) {
      slog("DEBUG", paste("  Predicting scenario:", pred_scenario))
      target_pred_file_path <- construct_prediction_filename(species_name_sanitized, pred_scenario, predictor_type_suffix, config)
      if (is.null(target_pred_file_path)) { slog("WARN", "  Could not construct prediction filename. Skipping."); prediction_errors <- prediction_errors + 1; next }
      
      if (!config$force_rerun$run_enmeval_sdms && file.exists(target_pred_file_path)) { slog("DEBUG", "  Prediction exists in target location. Skipping."); next } # Check ENMeval flag
      
      # Load predictor stack for the specific scenario
      pred_predictor_stack <- NULL
      if(use_pca) {
        pred_predictor_path <- predictor_paths_or_list[[pred_scenario]]
        if (is.null(pred_predictor_path) || !file.exists(pred_predictor_path)) { slog("WARN", paste("  PCA stack path missing for prediction:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
        pred_predictor_stack <- tryCatch(terra::rast(pred_predictor_path), error = function(e) {slog("WARN", paste("  Error loading PCA stack for prediction:", e$message)); NULL})
      } else {
        scenario_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, pred_scenario, config)
        if(length(scenario_vif_vars) < 1) { slog("WARN", paste("  No VIF vars found for prediction scenario:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
        pred_predictor_stack <- load_selected_env_data(pred_scenario, scenario_vif_vars, config)
      }
      if(is.null(pred_predictor_stack)) { slog("WARN", paste("  Failed load predictor stack for prediction:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      
      # Check if predictor names match model expectations (basic check)
      expected_vars <- tryCatch(names(optimal_model@variable.importance), error = function(e) NULL) # MaxNet specific
      if(is.null(expected_vars)) expected_vars <- tryCatch(optimal_model@model$var.names, error = function(e) NULL) # Fallback attempt
      if(!is.null(expected_vars) && !all(expected_vars %in% names(pred_predictor_stack))){
        missing_pred_vars <- setdiff(expected_vars, names(pred_predictor_stack))
        slog("ERROR", paste(" Predictor variables missing in stack for scenario", pred_scenario, "that are expected by the model:", paste(missing_pred_vars, collapse=", ")))
        prediction_errors <- prediction_errors + 1
        rm(pred_predictor_stack); gc(); next
      } else if (!is.null(expected_vars)) {
        # Reorder stack if necessary
        pred_predictor_stack <- pred_predictor_stack[[expected_vars]]
      }
      
      # Predict using the optimal model
      prediction_output <- predict_sdm_suitability(optimal_model, pred_predictor_stack, config, logger=NULL, species_log_file=species_log_file)
      
      # Save prediction
      if (inherits(prediction_output, "SpatRaster")) {
        save_success <- save_sdm_prediction(prediction_output, species_name_sanitized, pred_scenario, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)
        if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
        rm(prediction_output); gc()
      } else {
        slog("WARN", paste("  Prediction failed for scenario", pred_scenario, ":", prediction_output))
        prediction_errors <- prediction_errors + 1
      }
      rm(pred_predictor_stack); gc()
    } # End prediction scenario loop
  } else {
    prediction_errors <- length(scenarios_to_predict); msg <- paste0("Skipping all predictions: Optimal model object was not available."); slog("ERROR", msg)
    status <- if(load_existing_model) "error_loading_model" else "error_optimal_model_null"; return(list(status = status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
  }
  
  # --- Prepare return status ---
  final_status <- "success"
  status_message <- paste0("Finished ENMeval SDM. Occs:", occurrence_count_after_thinning, ". Preds attempted:", length(scenarios_to_predict), ". Made:", predictions_made, ". Errors/Skipped:", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) { final_status <- "error_prediction_all"; slog("ERROR", status_message) }
  else if (prediction_errors > 0) { final_status <- "success_with_pred_errors"; slog("WARN", status_message) }
  else { slog("INFO", status_message) }
  
  # --- Clean up ---
  if (!is.null(tuning_predictor_stack)) rm(tuning_predictor_stack)
  if (!is.null(occs_coords_df)) rm(occs_coords_df)
  if (!is.null(background_points_df)) rm(background_points_df)
  if (!load_existing_model && !is.null(enmeval_results)) rm(enmeval_results)
  gc()
  
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm_enmeval function


# --- 10. Setup Parallel Backend & Run ---
if (config$use_parallel && config$num_cores > 1 && !(config$enmeval_num_cores > 1) ) {
  log4r::info(logger, paste("Setting up parallel backend (furrr) with", config$num_cores, "cores."))
  gc(full=TRUE); future::plan(future::multisession, workers = config$num_cores, gc = TRUE)
} else if (config$use_parallel && config$enmeval_num_cores > 1) {
  log4r::info(logger, paste("Using internal ENMeval parallel with", config$enmeval_num_cores, "cores. Running species sequentially."))
  future::plan(future::sequential)
} else {
  log4r::info(logger, "Running sequentially.")
  future::plan(future::sequential)
}
progressr::handlers(global = TRUE); progressr::handlers("txtprogressbar")

log4r::info(logger, paste("Starting ENMeval SDM processing for", nrow(species_df), group_name, "species..."))

run_map_function <- if (inherits(future::plan(), "sequential") || inherits(future::plan(), "uniprocess")) {
  purrr::map
} else {
  furrr::future_map
}

results_list <- progressr::with_progress({
  run_map_function(1:nrow(species_df), ~{
    # Add tryCatch around the whole species processing for robustness in parallel
    tryCatch({
      process_species_sdm_enmeval(
        species_row = species_df[.x, ],
        config = config,
        predictor_paths_or_list = predictor_paths_or_list,
        group_name = group_name,
        predictor_type_suffix = predictor_type_suffix,
        use_pca = use_pca,
        occurrence_dir = occurrence_dir,
        tuning_scenario = "current"
      )
    }, error = function(e) {
      # Log the error and return an error status
      species_name_err <- species_df$scientificName[.x] %||% paste("Row", .x)
      error_msg <- paste("FATAL error processing species", species_name_err, ":", e$message)
      log4r::error(logger, error_msg) # Log to main logger
      # Also attempt to log to species file if possible
      species_log_file_err <- file.path(config$species_log_dir, paste0(gsub(" ", "_", species_name_err), predictor_type_suffix, "_detail.log"))
      cat(paste(Sys.time(), "[ERROR]", error_msg, "\n"), file = species_log_file_err, append = TRUE)
      # Return a consistent error structure
      return(list(status = "error_fatal_process", species = species_name_err, occurrence_count = NA, message = error_msg))
    })
  }, .options = furrr_options(seed = TRUE, scheduling = 1.0))
})

log4r::info(logger, "Processing loop complete.")

# --- 11. Process Results ---
# (Result processing and summary logic remains the same)
occurrence_counts <- list()
success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0
log4r::info(logger, paste("--- Processing Results Summary (", group_name, predictor_type_suffix, "Models - ENMeval) ---"))
for (res in results_list) {
  if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "Received NULL result."); next }
  log_level_func <- switch(res$status,
                           success = log4r::info,
                           success_with_pred_errors = log4r::warn,
                           skipped_occurrences = log4r::warn,
                           log4r::error) # Default to ERROR
  
  log_level_func(logger, paste0("[", res$species %||% "UNKNOWN", "] ", res$message %||% "No message.")) # Add species name to main log summary
  occurrence_counts[[res$species %||% paste("Unknown_",error_count)]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
  if (res$status == "success") success_count <- success_count + 1
  else if (res$status == "success_with_pred_errors") partial_success_count <- partial_success_count + 1
  else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
  else error_count <- error_count + 1
}
log4r::info(logger, paste("--- Overall Summary (", group_name, predictor_type_suffix, "Models - ENMeval) ---"))
log4r::info(logger, paste("Total Species Targets:", length(results_list)))
log4r::info(logger, paste("Fully Successful Runs:", success_count))
log4r::info(logger, paste("Partially Successful Runs (Prediction Errors):", partial_success_count))
log4r::info(logger, paste("Skipped (Occs/Input):", skipped_count))
log4r::info(logger, paste("Errors (Tuning/Processing/Fatal):", error_count))

if (length(occurrence_counts) > 0) {
  occ_count_df <- data.frame( Species = names(occurrence_counts), OccurrenceCountAfterThinning = unlist(occurrence_counts)) %>% dplyr::arrange(Species)
  occ_count_df_print <- occ_count_df; occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
  occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, ".csv")) # Suffix reflects ENMeval run
  tryCatch({ readr::write_csv(occ_count_df, occ_count_file); log4r::info(logger, paste("Occurrence counts saved:", occ_count_file))}, error = function(e) { log4r::error(logger, paste("Failed save counts:", e$message)) })
  cat("\n--- Final Occurrence Counts (", group_name, predictor_type_suffix, ", after cleaning/thinning) ---\n"); print(occ_count_df_print)
} else { log4r::warn(logger, "No occurrence counts were recorded.") }

future::plan(future::sequential); gc(full=TRUE) # Reset parallel plan
log4r::info(logger, "--- Script 06a (ENMeval Workflow v1.1) finished. ---")
#-------------------------------------------------------------------------------