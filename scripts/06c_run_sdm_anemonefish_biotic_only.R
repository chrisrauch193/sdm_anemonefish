# scripts/06c_run_sdm_anemonefish_biotic_pc1.R
#-------------------------------------------------------------------------------
# Run SDMs for Anemonefish Species using PC1 + Max Host Suitability
# (SDMtune Workflow). Aims to primarily capture biotic influence + main env gradient.
# Includes Parallel Processing, Logging, Progress, and Species-Specific Logs.
# Assumes 06a produced host predictions for ALL scenarios.
# V4: True Biotic + PC1 Implementation.
#-------------------------------------------------------------------------------
cat("--- Running Script 06c: Run Anemonefish Biotic+PC1 SDMs (v4) ---\n")

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

log4r::info(logger, "--- Starting Script 06c: Run Anemonefish Biotic+PC1 SDMs (v4) ---")

# --- 3. Define Group Specifics & Predictor Type ---
group_name <- "anemonefish"
species_list_file <- config$anemonefish_species_list_file
occurrence_dir <- config$anemonefish_occurrence_dir
host_group_name <- "anemone"

# Verify PCA is being used
if(!config$use_pca_predictors) {
  log4r::fatal(logger, "Script 06c (Biotic+PC1) requires config$use_pca_predictors=TRUE.")
  stop("This script requires PCA environmental predictors.")
}
host_predictor_type_suffix <- "_pca" # Suffix of host predictions to load

# Suffix for THIS script's outputs
predictor_type_suffix <- "_biotic_pc1" # <--- CHANGED SUFFIX

log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors: Biotic (Max Host Suitability", host_predictor_type_suffix, ") + PC1 ---")) # <--- UPDATED MESSAGE

# --- 4. Load Environmental Predictor Paths (PCA) ---
# Need paths to get PC1 layer for each scenario
log4r::debug(logger, "Loading Environmental (PCA) paths...")
pca_paths_rds <- config$pca_raster_paths_rds_path
if (!file.exists(pca_paths_rds)) {
  log4r::fatal(logger, paste("PCA raster paths file not found:", pca_paths_rds)); stop("PCA paths missing.")
}
env_predictor_paths_list <- tryCatch(readRDS(pca_paths_rds), error = function(e) {
  log4r::fatal(logger, paste("Failed load PCA paths:", e$message)); NULL
})
if (!is.list(env_predictor_paths_list) || length(env_predictor_paths_list) == 0) {
  log4r::fatal(logger, "PCA paths list empty/invalid."); stop("PCA paths list invalid.")
}
scenarios_to_process <- names(env_predictor_paths_list)
log4r::debug(logger, paste("Loaded Env PCA paths for scenarios:", paste(scenarios_to_process, collapse=", ")))

# --- 5. Load Processed Host Association Data ---
association_file_path <- file.path(config$data_dir, "processed_anemonefish_host_associations.csv")
if(!file.exists(association_file_path)) {
  log4r::fatal(logger, "Association file missing."); stop("Association file missing.")
}
associations_df <- tryCatch(readr::read_csv(association_file_path, show_col_types = FALSE), error = function(e){
  log4r::fatal(logger, paste("Failed load associations:", e$message)); NULL
})
if(is.null(associations_df)) stop("Failed to load associations file.")
log4r::info(logger, paste("Loaded", nrow(associations_df), "fish-host associations."))

# --- 6. Get File Paths for ALL Host Anemone Predictions (ALL Scenarios) ---
log4r::info(logger, "Gathering file paths for ALL potential host anemone predictions...")
host_predictions_dir <- file.path(config$predictions_dir, paste0(host_group_name, host_predictor_type_suffix))
if(!dir.exists(host_predictions_dir)) {
  log4r::fatal(logger, paste("Host predictions directory not found:", host_predictions_dir)); stop("Host predictions missing.")
}
host_species_list_file <- config$anemone_species_list_file
if (!file.exists(host_species_list_file)) stop("Anemone list missing.")
host_species_df <- readr::read_csv(host_species_list_file, show_col_types = FALSE)

all_host_prediction_paths <- list()
for (scenario in scenarios_to_process) {
  all_host_prediction_paths[[scenario]] <- list()
  log4r::debug(logger, paste("Checking host prediction paths for scenario:", scenario))
  missing_files_scenario = 0
  for (i in 1:nrow(host_species_df)) {
    host_sci_name <- host_species_df$scientificName[i]
    host_sci_name_sanitized <- gsub(" ", "_", host_sci_name)
    pred_file_path <- file.path(host_predictions_dir, paste0("sdm_prediction_", host_sci_name_sanitized, "_", scenario, host_predictor_type_suffix, ".tif"))
    if (file.exists(pred_file_path)) {
      all_host_prediction_paths[[scenario]][[host_sci_name_sanitized]] <- pred_file_path
    } else {
      log4r::warn(logger, paste("Host prediction file missing:", basename(pred_file_path)))
      missing_files_scenario = missing_files_scenario + 1
    }
  }
  if(length(all_host_prediction_paths[[scenario]]) == 0 && missing_files_scenario > 0){
    log4r::error(logger, paste("FATAL: No host predictions found for ANY species in scenario:", scenario))
  } else if (missing_files_scenario > 0) {
    log4r::warn(logger, paste("Found", length(all_host_prediction_paths[[scenario]]), "host prediction files, but", missing_files_scenario, "were missing for scenario:", scenario))
  } else {
    log4r::info(logger, paste("Found all", length(all_host_prediction_paths[[scenario]]), "host prediction files for scenario:", scenario))
  }
}

# --- 7. Create Output Directories (Using _biotic_pc1 suffix) ---
group_pred_dir <- file.path(config$predictions_dir, paste0(group_name, predictor_type_suffix))
group_results_dir <- file.path(config$results_dir, paste0(group_name, predictor_type_suffix))
group_models_dir <- file.path(config$models_dir, paste0(group_name, predictor_type_suffix))
species_log_dir <- config$species_log_dir

dir.create(group_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_models_dir, recursive = TRUE, showWarnings = FALSE)

log4r::debug(logger, paste("Output directories created/checked for:", group_name, predictor_type_suffix))

# --- 8. Load Anemonefish Species List ---
tryCatch({
  species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) {
  log4r::fatal(logger, paste("Failed load species list:", e$message)); stop("Species list load failed.")
})
log4r::info(logger, paste("Loaded", nrow(species_df), "anemonefish species."))

# --- 9. Define Function to Process Single Species (Biotic + PC1) ---
process_species_sdm_biotic_pc1 <- function(species_row, config, env_predictor_paths_list, all_host_prediction_paths, associations_df,
                                           group_pred_dir, group_results_dir, group_models_dir,
                                           species_log_dir, predictor_type_suffix, occurrence_dir,
                                           tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  log_prefix <- paste0("[", species_name, "] ")
  species_log_file <- file.path(species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), log_prefix, paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
  slog("INFO", "--- Starting Biotic+PC1 processing ---")
  
  # --- Get Associated Hosts ---
  associated_hosts <- associations_df %>%
    filter(AnemonefishScientificName == species_name) %>%
    pull(AssociatedAnemoneScientificName)
  associated_hosts_sanitized <- gsub(" ", "_", associated_hosts)
  if (length(associated_hosts_sanitized) == 0) {
    msg <- paste0(log_prefix, "Skipping: No associated hosts found."); slog("WARN", msg)
    return(list(status = "skipped_no_hosts", species = species_name, occurrence_count = NA, message = msg))
  }
  slog("DEBUG", paste("Associated hosts:", paste(associated_hosts_sanitized, collapse=", ")))
  
  # --- Prepare Biotic+PC1 Predictor Stack for Tuning Scenario ---
  slog("DEBUG", "Preparing Biotic+PC1 predictor stack for tuning scenario:", tuning_scenario)
  
  # 1. Load PC1 layer from Env (PCA) predictors for tuning
  tuning_env_path <- env_predictor_paths_list[[tuning_scenario]]
  if (is.null(tuning_env_path) || !file.exists(tuning_env_path)) {
    msg <- paste0(log_prefix, "Skipping: Env (PCA) stack path missing for tuning."); slog("ERROR", msg)
    return(list(status = "error_tuning_env_pca_path", species = species_name, occurrence_count = NA, message = msg))
  }
  # Load only the first layer (PC1)
  pc1_tuning_layer <- tryCatch(terra::rast(tuning_env_path)[[1]], error = function(e) { slog("ERROR", "Failed load tuning PC1 layer:", e$message); NULL })
  if (is.null(pc1_tuning_layer)) { msg <- paste0(log_prefix, "Skipping: Failed load tuning PC1 layer."); slog("ERROR", msg); return(list(status = "error_tuning_env_pc1_load", species = species_name, occurrence_count = NA, message = msg)) }
  reference_geom_tuning <- pc1_tuning_layer # Use PC1 for geometry reference
  
  # 2. Get host paths for tuning
  host_paths_for_tuning <- all_host_prediction_paths[[tuning_scenario]][names(all_host_prediction_paths[[tuning_scenario]]) %in% associated_hosts_sanitized]
  if (length(host_paths_for_tuning) == 0) {
    msg <- paste0(log_prefix, "Skipping: No host prediction paths found for tuning."); slog("ERROR", msg)
    return(list(status = "error_tuning_biotic_paths", species = species_name, occurrence_count = NA, message = msg))
  }
  if (length(host_paths_for_tuning) < length(associated_hosts_sanitized)) { slog("WARN", "Missing paths for some hosts in tuning.") }
  
  # 3. Load, Resample, and Stack Host predictors for tuning
  host_rasters_tuning_list <- list()
  for(host_name in names(host_paths_for_tuning)){
    host_path <- host_paths_for_tuning[[host_name]]
    host_rast <- tryCatch(terra::rast(host_path), error=function(e){slog("WARN", paste("Failed load host", basename(host_path),":", e$message)); NULL})
    if(!is.null(host_rast)){
      if (!terra::compareGeom(reference_geom_tuning, host_rast, stopOnError=FALSE)) {
        slog("DEBUG", paste("Resampling host", basename(host_path), "for tuning."))
        host_rast <- tryCatch(terra::resample(host_rast, reference_geom_tuning, method="bilinear"), error = function(e) { slog("ERROR", paste("Failed resample host", basename(host_path),":", e$message)); NULL })
      }
      if(!is.null(host_rast)){ names(host_rast) <- host_name; host_rasters_tuning_list[[host_name]] <- host_rast }
    }
  }
  if (length(host_rasters_tuning_list) == 0) {
    msg <- paste0(log_prefix, "Skipping: Failed load/resample ANY hosts for tuning."); slog("ERROR", msg)
    return(list(status = "error_tuning_biotic_load", species = species_name, occurrence_count = NA, message = msg))
  }
  host_stack_tuning <- tryCatch(terra::rast(host_rasters_tuning_list), error = function(e) { slog("ERROR", paste("Failed stack host predictors:", e$message)); NULL })
  if (is.null(host_stack_tuning) || terra::nlyr(host_stack_tuning) == 0) {
    msg <- paste0(log_prefix, "Skipping: Failed create valid host stack for tuning."); slog("ERROR", msg)
    return(list(status = "error_tuning_biotic_stack", species = species_name, occurrence_count = NA, message = msg))
  }
  
  # 4. Calculate Max Host Suitability for tuning
  max_host_suitability_tuning <- tryCatch(terra::app(host_stack_tuning, fun = "max", na.rm = TRUE), error = function(e) { slog("ERROR", paste("Failed calc max host suitability:", e$message)); NULL })
  if(is.null(max_host_suitability_tuning)) { msg <- paste0(log_prefix, "Skipping: Failed calc max host."); slog("ERROR", msg); return(list(status = "error_tuning_biotic_max", species = species_name, occurrence_count = NA, message = msg)) }
  names(max_host_suitability_tuning) <- "host_suitability_max"
  
  # 5. Combine PC1 + Biotic for tuning stack
  tuning_predictor_stack <- tryCatch(c(pc1_tuning_layer, max_host_suitability_tuning), error = function(e) { slog("ERROR", paste("Failed combine PC1+Biotic for tuning:", e$message)); NULL })
  if(is.null(tuning_predictor_stack)) { msg <- paste0(log_prefix, "Skipping: Failed combine stacks for tuning."); slog("ERROR", msg); return(list(status = "error_tuning_combine_stack", species = species_name, occurrence_count = NA, message = msg)) }
  # Ensure names are correct (might get reset by c())
  names(tuning_predictor_stack) <- c("PC1", "host_suitability_max")
  slog("INFO", paste("Biotic+PC1 tuning stack created with layers:", paste(names(tuning_predictor_stack), collapse=", ")))
  
  # --- Load, Clean, Thin Occurrences (using Biotic+PC1 tuning stack) ---
  config_for_occ_load <- config
  config_for_occ_load$predictor_stack_for_thinning <- tuning_predictor_stack # Use Biotic+PC1 stack
  slog("DEBUG", "Loading/cleaning/thinning occurrences using Biotic+PC1 tuning stack.")
  occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL)
  if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) {
    msg <- paste0(log_prefix, "Skipping: Insufficient occurrences (", occ_data_result$count %||% 0, ")."); slog("WARN", msg)
    return(list(status = "skipped_occurrences", species = species_name, occurrence_count = occ_data_result$count %||% 0, message = msg))
  }
  occs_coords <- occ_data_result$coords
  occurrence_count_after_thinning <- occ_data_result$count
  slog("INFO", "Occurrence count after clean/thin:", occurrence_count_after_thinning)
  
  # --- Define File Paths ---
  tuning_results_file <- file.path(group_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, ".rds"))
  final_model_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  # --- Check if Final Model Exists ---
  if (!config$force_rerun$run_biotic_sdms && file.exists(final_model_file)) { # Check biotic flag
    slog("INFO", "Skipping tuning/training: Final biotic+pc1 model exists.")
    load_existing_model <- TRUE
    final_model <- tryCatch(readRDS(final_model_file), error=function(e) { slog("ERROR", "Failed load existing model:", e$message); NULL})
    if(is.null(final_model) || !inherits(final_model, "SDMmodel")){
      msg <- paste0(log_prefix, "Existing final model invalid. Will re-run."); slog("WARN", msg); load_existing_model <- FALSE
    }
  } else { load_existing_model <- FALSE }
  
  # --- Run Tuning & Training if needed (Using Biotic+PC1 Stack) ---
  if (!load_existing_model) {
    slog("INFO", "Final biotic+pc1 model not found or rerun forced. Tuning/training...")
    background_points <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config, logger = NULL, seed = species_aphia_id)
    if (is.null(background_points)) {
      msg <- paste0(log_prefix, "Skipping: Failed background points."); slog("ERROR", msg)
      return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("INFO", "Starting hyperparameter tuning (using biotic+pc1 predictor).")
    tuning_output <- run_sdm_tuning_kfold(occs_coords, tuning_predictor_stack, background_points, config, logger = NULL, species_name)
    if (is.null(tuning_output) || is.null(tuning_output$best_hypers)) {
      # *** MORE SPECIFIC ERROR LOGGING FROM TUNING HELPER ***
      tuning_error_msg <- attributes(tuning_output)$error_message %||% "Unknown tuning error"
      msg <- paste0(log_prefix, "Skipping: Tuning failed. Reason: ", tuning_error_msg); slog("ERROR", msg)
      return(list(status = "error_tuning", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    best_hypers <- tuning_output$best_hypers
    slog("INFO", "Tuning complete. Best hypers:", paste(names(best_hypers), best_hypers[1,], collapse=", "))
    save_tuning_result <- tryCatch({ saveRDS(tuning_output, tuning_results_file); slog("DEBUG", "Tuning results saved."); NULL }, error = function(e){ slog("ERROR", "Failed save tuning results:", e$message); return(list(status = "error_saving_tuning_results", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0(log_prefix, "Failed save tuning: ", e$message))) })
    if (!is.null(save_tuning_result)) return(save_tuning_result)
    
    slog("INFO", "Starting final model training (using biotic+pc1 predictor).")
    final_model <- train_final_sdm(occs_coords, tuning_predictor_stack, background_points, best_hypers, config, logger = NULL, species_name)
    if(is.null(final_model)){
      msg <- paste0(log_prefix, "Skipping: Training failed."); slog("ERROR", msg)
      return(list(status = "error_training", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("INFO", "Final biotic+pc1 model training complete.")
    save_model_result <- tryCatch({ saveRDS(final_model, final_model_file); slog("DEBUG", "Final biotic+pc1 model saved."); NULL }, error = function(e){ slog("ERROR", "Failed save final model:", e$message); return(list(status = "error_saving_model", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0(log_prefix, "Failed save final model: ", e$message))) })
    if (!is.null(save_model_result)) return(save_model_result)
    
    rm(background_points, tuning_output); gc()
  } # End if !load_existing_model
  
  # --- Prediction Loop (Using Biotic+PC1 Stack ONLY) ---
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- names(env_predictor_paths_list) # Predict for all scenarios where we have geometry
  slog("INFO", paste("Starting predictions for", length(scenarios_to_predict), "scenarios using biotic+pc1 model."))
  
  if (!is.null(final_model)) {
    for (pred_scenario in scenarios_to_predict) {
      slog("DEBUG", paste("  Predicting scenario:", pred_scenario))
      pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", pred_scenario, predictor_type_suffix, ".tif"))
      
      if (!config$force_rerun$run_biotic_sdms && file.exists(pred_file)) { slog("DEBUG", "  Prediction exists. Skipping."); next }
      
      # 1. Load PC1 layer for prediction scenario
      pred_env_path <- env_predictor_paths_list[[pred_scenario]]
      if (is.null(pred_env_path) || !file.exists(pred_env_path)) { slog("WARN", paste("  Env stack missing (for PC1):", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      pc1_pred_layer <- tryCatch(terra::rast(pred_env_path)[[1]], error = function(e) { slog("WARN", paste("  Error load env PC1:", e$message)); NULL })
      if(is.null(pc1_pred_layer)) { prediction_errors <- prediction_errors + 1; next }
      reference_geom_pred <- pc1_pred_layer # Use scenario's PC1 for geometry
      
      # 2. Get host paths for prediction scenario
      host_paths_for_pred <- all_host_prediction_paths[[pred_scenario]][names(all_host_prediction_paths[[pred_scenario]]) %in% associated_hosts_sanitized]
      if (length(host_paths_for_pred) == 0) { slog("WARN", paste("  No host paths found:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      
      # 3. Load/Resample/Stack Hosts for prediction
      host_rasters_pred_list <- list()
      for(host_name in names(host_paths_for_pred)){
        host_path <- host_paths_for_pred[[host_name]]
        host_rast <- tryCatch(terra::rast(host_path), error=function(e){slog("WARN",paste("Failed load host",basename(host_path),":",e$message));NULL})
        if(!is.null(host_rast)){
          if (!terra::compareGeom(reference_geom_pred, host_rast, stopOnError=FALSE)) {
            slog("DEBUG", paste("Resampling host", basename(host_path)))
            host_rast <- tryCatch(terra::resample(host_rast, reference_geom_pred, method="bilinear"), error=function(e){slog("ERROR",paste("Failed resample",basename(host_path),":",e$message));NULL})
          }
          if(!is.null(host_rast)){ names(host_rast) <- host_name; host_rasters_pred_list[[host_name]] <- host_rast }
        }
      }
      if (length(host_rasters_pred_list) == 0) { slog("WARN", paste("  Failed load/resample hosts:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      host_stack_scenario <- tryCatch(terra::rast(host_rasters_pred_list), error = function(e) { slog("WARN", paste("  Error stacking hosts:", e$message)); NULL })
      if (is.null(host_stack_scenario)) { prediction_errors <- prediction_errors + 1; next }
      
      # 4. Calculate Max Host Suitability for prediction
      max_host_suitability_scenario <- tryCatch(terra::app(host_stack_scenario, fun = "max", na.rm = TRUE), error = function(e) { slog("WARN", paste("  Error calc max host:", e$message)); NULL })
      if(is.null(max_host_suitability_scenario)) { prediction_errors <- prediction_errors + 1; next }
      names(max_host_suitability_scenario) <- "host_suitability_max"
      
      # 5. Combine PC1 + Biotic for prediction stack
      prediction_predictor_stack <- tryCatch(c(pc1_pred_layer, max_host_suitability_scenario), error=function(e){slog("WARN", paste("  Error combining PC1+Biotic:", e$message)); NULL})
      if(is.null(prediction_predictor_stack)) { prediction_errors <- prediction_errors + 1; next }
      names(prediction_predictor_stack) <- c("PC1", "host_suitability_max") # Ensure names
      
      # 6. Predict & Save
      prediction_output <- predict_sdm_suitability(final_model, prediction_predictor_stack, config, logger = NULL)
      if (inherits(prediction_output, "SpatRaster") && !is.null(prediction_output) && terra::nlyr(prediction_output) > 0) {
        save_success <- tryCatch({ terra::writeRaster(prediction_output, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES")); slog("DEBUG", paste("  Prediction raster saved:", basename(pred_file))); TRUE }, error = function(e) { slog("ERROR", paste("  Failed save prediction:", e$message)); FALSE })
        if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
      } else { slog("WARN", paste("  Prediction failed. Output:", class(prediction_output))); prediction_errors <- prediction_errors + 1 }
      
      rm(pc1_pred_layer, host_stack_scenario, max_host_suitability_scenario, prediction_predictor_stack, prediction_output, host_rasters_pred_list); gc()
    } # End prediction scenario loop
  } else {
    prediction_errors <- length(scenarios_to_predict); msg <- paste0(log_prefix, "Skipping predictions: Final model unavailable."); slog("ERROR", msg)
    status <- if(load_existing_model) "error_loading_model" else "error_training"; return(list(status = status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
  }
  
  # --- Prepare return status ---
  final_status <- "success"
  status_message <- paste0(log_prefix, "Finished Biotic+PC1. Occs:", occurrence_count_after_thinning, ". Preds attempted:", length(scenarios_to_predict), ". Made:", predictions_made, ". Errors/Skipped:", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) { final_status <- "error_prediction_all"; slog("ERROR", status_message) }
  else if (prediction_errors > 0) { final_status <- "success_with_pred_errors"; slog("WARN", status_message) }
  else { slog("INFO", status_message) }
  
  rm(tuning_predictor_stack, occs_coords, host_stack_tuning, max_host_suitability_tuning, host_rasters_tuning_list, pc1_tuning_layer); gc()
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm_biotic_pc1 function


# --- 10. Setup Parallel Backend & Run ---
if (config$use_parallel && config$num_cores > 1) {
  log4r::info(logger, paste("Setting up parallel backend with", config$num_cores, "cores."))
  gc(full=TRUE)
  future::plan(future::multisession, workers = config$num_cores, gc = TRUE)
} else {
  log4r::info(logger, "Running sequentially.")
  future::plan(future::sequential)
}
progressr::handlers(global = TRUE); progressr::handlers("txtprogressbar")

log4r::info(logger, paste("Starting Biotic+PC1 SDM processing for", nrow(species_df), "anemonefish..."))

results_list <- progressr::with_progress({
  furrr::future_map(1:nrow(species_df), ~{
    process_species_sdm_biotic_pc1( # Call the Biotic+PC1 function
      species_row = species_df[.x, ],
      config = config,
      env_predictor_paths_list = env_predictor_paths_list, # Still needed for PC1
      all_host_prediction_paths = all_host_prediction_paths,
      associations_df = associations_df,
      group_pred_dir = group_pred_dir,
      group_results_dir = group_results_dir,
      group_models_dir = group_models_dir,
      species_log_dir = config$species_log_dir,
      predictor_type_suffix = predictor_type_suffix, # _biotic_pc1
      occurrence_dir = occurrence_dir
    )
  }, .options = furrr_options(seed = TRUE, scheduling = 1.0))
})

log4r::info(logger, "Parallel/sequential processing complete.")

# --- 11. Process Results ---
occurrence_counts <- list()
success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0

log4r::info(logger, "--- Processing Results Summary (Anemonefish Biotic+PC1 Models) ---")
for (res in results_list) {
  if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "Received NULL result."); next }
  log_level_func <- switch(res$status, success=log4r::info, success_with_pred_errors=log4r::warn, skipped_occurrences=log4r::warn, skipped_no_hosts=log4r::warn, log4r::error)
  log_level_func(logger, res$message)
  occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
  if (res$status == "success") success_count <- success_count + 1
  else if (res$status == "success_with_pred_errors") partial_success_count <- partial_success_count + 1
  else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
  else error_count <- error_count + 1
}

log4r::info(logger, paste("--- Overall Summary (Anemonefish Biotic+PC1 Models) ---"))
log4r::info(logger, paste("Total Species Targets:", length(results_list)))
log4r::info(logger, paste("Fully Successful Runs:", success_count))
log4r::info(logger, paste("Partially Successful Runs:", partial_success_count))
log4r::info(logger, paste("Skipped:", skipped_count))
log4r::info(logger, paste("Errors:", error_count))

if (length(occurrence_counts) > 0) {
  occ_count_df <- data.frame( Species = names(occurrence_counts), OccurrenceCountAfterThinning = unlist(occurrence_counts)) %>% dplyr::arrange(Species)
  occ_count_df_print <- occ_count_df; occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
  occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, ".csv"))
  tryCatch({ readr::write_csv(occ_count_df, occ_count_file); log4r::info(logger, paste("Occurrence counts saved:", occ_count_file))}, error = function(e) { log4r::error(logger, paste("Failed save counts:", e$message)) })
  cat("\n--- Final Occurrence Counts (Anemonefish Biotic+PC1, after cleaning/thinning) ---\n"); print(occ_count_df_print)
} else { log4r::warn(logger, "No occurrence counts recorded.") }

# --- Shutdown & Finish ---
future::plan(future::sequential)
gc(full=TRUE)

log4r::info(logger, "--- Script 06c (Biotic+PC1 v4) finished. ---")
#-------------------------------------------------------------------------------