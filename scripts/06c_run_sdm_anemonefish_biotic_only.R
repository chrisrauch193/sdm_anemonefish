# scripts/06c_run_sdm_anemonefish_biotic_only.R
#-------------------------------------------------------------------------------
# Run Biotic-Only SDMs for Anemonefish Species using Host Anemone
# Suitability Predictions (SDMtune Workflow) with Parallel, Logging, Progress
# Includes species-specific detailed log files.
#-------------------------------------------------------------------------------
cat("--- Running Script 06c: Run Anemonefish Biotic-Only SDMs (Parallel & Logged) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr, log4r, future, furrr, progressr)

# Source helpers
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R")) # Needed for load_selected_env_data if using VIF hosts
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))

# --- 2. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path,
                       log_level = config$log_level,
                       append = config$log_append,
                       log_to_console = config$log_to_console,
                       console_level = config$log_console_level)

log4r::info(logger, "--- Starting Script 06c: Run Anemonefish Biotic-Only SDMs ---")

# --- 3. Define Group Specifics & Predictor Type ---
group_name <- "anemonefish"
species_list_file <- config$anemonefish_species_list_file
occurrence_dir <- config$anemonefish_occurrence_dir
host_group_name <- "anemone" # Group name for host predictions
host_predictor_type_suffix <- ifelse(config$use_pca_predictors, "_pca", "_vif") # Suffix of host predictions
predictor_type_suffix <- "_hostonly" # Specific suffix for this run

log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors: Host Anemone Suitability Only", host_predictor_type_suffix ,"---")) # Clarify source

# --- 4. Load Processed Host Association Data ---
association_file_path <- file.path(config$data_dir, "processed_anemonefish_host_associations.csv") # Path to processed file
if(!file.exists(association_file_path)) {
  log4r::fatal(logger, "Processed association file missing. Run the processing snippet first.")
  stop("Processed association file not found at: ", association_file_path)
}
associations_df <- readr::read_csv(association_file_path, show_col_types = FALSE)
log4r::info(logger, paste("Loaded", nrow(associations_df), "fish-host associations."))

# --- 5. Load ALL Host Prediction Rasters (for Current Scenario) ---
log4r::info(logger, "Loading CURRENT prediction rasters for ALL potential host anemone species...")
host_predictions_dir <- file.path(config$predictions_dir, paste0(host_group_name, host_predictor_type_suffix))
if(!dir.exists(host_predictions_dir)) {
  log4r::fatal(logger, paste("Host predictions directory not found:", host_predictions_dir))
  stop("Host prediction directory missing. Run script 06a first.")
}

host_species_list_file <- config$anemone_species_list_file # Get host species list
if (!file.exists(host_species_list_file)) stop("Anemone list missing.")
host_species_df <- readr::read_csv(host_species_list_file, show_col_types = FALSE)

all_host_predictions <- list()
reference_raster_geom <- NULL # To store geometry of the first loaded raster

for (i in 1:nrow(host_species_df)) {
  host_sci_name <- host_species_df$scientificName[i]
  host_sci_name_sanitized <- gsub(" ", "_", host_sci_name)
  # Construct path for CURRENT prediction only
  pred_file <- file.path(host_predictions_dir, paste0("sdm_prediction_", host_sci_name_sanitized, "_current", host_predictor_type_suffix, ".tif"))
  
  if (file.exists(pred_file)) {
    pred_rast <- tryCatch(terra::rast(pred_file), error = function(e) {
      log4r::warn(logger, paste("Failed to load host prediction for", host_sci_name, ":", e$message))
      NULL
    })
    if (!is.null(pred_rast)) {
      # Check and store reference geometry from the first valid raster
      if (is.null(reference_raster_geom)) {
        reference_raster_geom <- pred_rast
        log4r::debug(logger, paste("Set reference geometry from:", basename(pred_file)))
      } else {
        # Check geometry and resample if needed
        if (!terra::compareGeom(reference_raster_geom, pred_rast, stopOnError=FALSE)) {
          log4r::warn(logger, paste("Resampling host predictor", basename(pred_file), "to match reference geometry."))
          pred_rast <- tryCatch(terra::resample(pred_rast, reference_raster_geom, method="bilinear"), error = function(e) {
            log4r::error(logger, paste("Failed to resample", basename(pred_file), ":", e$message))
            NULL
          })
        }
      }
      # Only add if raster is valid (and possibly resampled)
      if (!is.null(pred_rast)) {
        names(pred_rast) <- host_sci_name_sanitized # Use sanitized name for stacking later
        all_host_predictions[[host_sci_name_sanitized]] <- pred_rast
      }
    }
  } else {
    log4r::warn(logger, paste("Host prediction file missing for CURRENT scenario:", basename(pred_file)))
  }
}

if (length(all_host_predictions) == 0) {
  log4r::fatal(logger, "No host prediction rasters were loaded successfully.")
  stop("Failed to load any host predictions.")
}
log4r::info(logger, paste("Successfully loaded", length(all_host_predictions), "host prediction rasters."))

# --- 6. Create Output Directories (Specific to this run) ---
group_pred_dir <- file.path(config$predictions_dir, paste0(group_name, predictor_type_suffix))
group_results_dir <- file.path(config$results_dir, paste0(group_name, predictor_type_suffix))
group_models_dir <- file.path(config$models_dir, paste0(group_name, predictor_type_suffix))
species_log_dir <- config$species_log_dir # Use the main species log dir

dir.create(group_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_models_dir, recursive = TRUE, showWarnings = FALSE)
# Species log dir created in 06a

log4r::debug(logger, paste("Output directories created/checked for:", group_name, predictor_type_suffix))

# --- 7. Load Anemonefish Species List ---
tryCatch({
  species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) {
  log4r::fatal(logger, paste("Failed load species list:", e$message))
  stop("Anemonefish species list loading failed.")
})
log4r::info(logger, paste("Loaded", nrow(species_df), "anemonefish species from", basename(species_list_file)))

# --- 8. Define Function to Process Single Species (Biotic-Only) ---
process_species_sdm_biotic_only <- function(species_row, config, all_host_predictions, associations_df,
                                            group_pred_dir, group_results_dir, group_models_dir,
                                            species_log_dir, predictor_type_suffix, occurrence_dir) {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  log_prefix <- paste0("[", species_name, "] ")
  species_log_file <- file.path(species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), log_prefix, paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
  slog("INFO", "--- Starting Biotic-Only processing ---")
  
  # --- Get Associated Hosts ---
  associated_hosts <- associations_df %>%
    filter(AnemonefishScientificName == species_name) %>%
    pull(AssociatedAnemoneScientificName)
  
  associated_hosts_sanitized <- gsub(" ", "_", associated_hosts)
  
  if (length(associated_hosts_sanitized) == 0) {
    msg <- paste0(log_prefix, "Skipping: No associated hosts found in the processed association file.")
    slog("WARN", msg); return(list(status = "skipped_no_hosts", species = species_name, occurrence_count = NA, message = msg))
  }
  slog("DEBUG", paste("Associated hosts:", paste(associated_hosts_sanitized, collapse=", ")))
  
  # --- Create Host-Specific Predictor Stack ---
  host_preds_for_species <- all_host_predictions[names(all_host_predictions) %in% associated_hosts_sanitized]
  
  if (length(host_preds_for_species) == 0) {
    msg <- paste0(log_prefix, "Skipping: None of the associated hosts had prediction rasters loaded.")
    slog("ERROR", msg); return(list(status = "error_missing_host_preds", species = species_name, occurrence_count = NA, message = msg))
  }
  if (length(host_preds_for_species) < length(associated_hosts_sanitized)) {
    slog("WARN", "Missing prediction rasters for some associated hosts.")
  }
  
  # Stack the found host predictors
  predictor_stack <- tryCatch(terra::rast(host_preds_for_species), error = function(e) {
    slog("ERROR", paste("Failed to stack host predictors:", e$message)); NULL
  })
  if (is.null(predictor_stack) || terra::nlyr(predictor_stack) == 0) {
    msg <- paste0(log_prefix, "Skipping: Failed to create valid host predictor stack.")
    slog("ERROR", msg); return(list(status = "error_host_stack", species = species_name, occurrence_count = NA, message = msg))
  }
  slog("INFO", paste("Using", terra::nlyr(predictor_stack), "host suitability layers as predictors."))
  
  # Define file paths for this specific species and run type
  tuning_results_file <- file.path(group_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, ".rds"))
  final_model_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_current", predictor_type_suffix, ".tif")) # Only current prediction
  
  # --- Load, Clean, Thin Occurrences (using host stack) ---
  config_for_occ_load <- config
  config_for_occ_load$predictor_stack_for_thinning <- predictor_stack # Use host stack for thinning
  slog("DEBUG", "Loading/cleaning/thinning occurrences using host stack.")
  occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL)
  if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) {
    msg <- paste0(log_prefix, "Skipping: Insufficient valid/thinned occurrences (found ", occ_data_result$count %||% 0, ", need ", config$min_occurrences_sdm, ").")
    slog("WARN", msg); return(list(status = "skipped_occurrences", species = species_name, occurrence_count = occ_data_result$count %||% 0, message = msg))
  }
  occs_coords <- occ_data_result$coords
  occurrence_count_after_thinning <- occ_data_result$count
  slog("INFO", "Occurrence count after clean/thin:", occurrence_count_after_thinning)
  
  # --- Check if Final Model Exists ---
  if (!config$force_rerun$run_biotic_sdms && file.exists(final_model_file)) {
    slog("INFO", "Skipping tuning/training: Final model exists.")
    load_existing_model <- TRUE
    final_model <- tryCatch(readRDS(final_model_file), error=function(e) { slog("ERROR", "Failed to load existing model:", e$message); NULL})
    if(is.null(final_model) || !inherits(final_model, "SDMmodel")){
      msg <- paste0(log_prefix, "Existing final model file is invalid or failed to load. Will attempt re-run.")
      slog("WARN", msg); load_existing_model <- FALSE
    }
  } else { load_existing_model <- FALSE }
  
  # --- Run Tuning & Training if needed ---
  if (!load_existing_model) {
    slog("INFO", "Final model not found or rerun forced. Proceeding with tuning/training.")
    slog("DEBUG", "Generating background points based on host stack.")
    background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config, logger = NULL, seed = species_aphia_id)
    if (is.null(background_points)) {
      msg <- paste0(log_prefix, "Skipping: Failed background point generation."); slog("ERROR", msg)
      return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    
    slog("INFO", "Starting hyperparameter tuning (using host predictors).")
    tuning_output <- run_sdm_tuning_kfold(occs_coords, predictor_stack, background_points, config, logger = NULL, species_name)
    if (is.null(tuning_output) || is.null(tuning_output$best_hypers)) {
      msg <- paste0(log_prefix, "Skipping: Hyperparameter tuning failed."); slog("ERROR", msg)
      return(list(status = "error_tuning", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    best_hypers <- tuning_output$best_hypers
    slog("INFO", "Tuning complete. Best hypers:", paste(names(best_hypers), best_hypers[1,], collapse=", "))
    save_tuning_result <- tryCatch({ saveRDS(tuning_output, tuning_results_file); slog("DEBUG", "Tuning results saved."); NULL }, error = function(e){ slog("ERROR", "Failed save tuning results:", e$message); return(list(status = "error_saving_tuning_results", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0(log_prefix, "Failed save tuning results: ", e$message))) })
    if (!is.null(save_tuning_result)) return(save_tuning_result)
    
    slog("INFO", "Starting final model training.")
    final_model <- train_final_sdm(occs_coords, predictor_stack, background_points, best_hypers, config, logger = NULL, species_name)
    if(is.null(final_model)){
      msg <- paste0(log_prefix, "Skipping: Final model training failed."); slog("ERROR", msg)
      return(list(status = "error_training", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("INFO", "Final model training complete.")
    save_model_result <- tryCatch({ saveRDS(final_model, final_model_file); slog("DEBUG", "Final model saved."); NULL }, error = function(e){ slog("ERROR", "Failed save final model:", e$message); return(list(status = "error_saving_model", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0(log_prefix, "Failed save final model: ", e$message))) })
    if (!is.null(save_model_result)) return(save_model_result)
    
    rm(background_points, tuning_output); gc()
  }
  
  # --- Prediction (Current Only) ---
  predictions_made = 0; prediction_errors = 0
  if (!is.null(final_model)) {
    if (!config$force_rerun$run_biotic_sdms && file.exists(pred_file)) {
      slog("DEBUG", "  Prediction exists. Skipping.")
    } else {
      slog("INFO", "Predicting onto current host suitability predictors...")
      # Predictor stack is already the correct host stack
      prediction_output <- predict_sdm_suitability(final_model, predictor_stack, config, logger = NULL)
      if (inherits(prediction_output, "SpatRaster") && !is.null(prediction_output) && terra::nlyr(prediction_output) > 0) {
        save_success <- tryCatch({ terra::writeRaster(prediction_output, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES")); slog("DEBUG", paste("  Prediction raster saved:", basename(pred_file))); TRUE}, error = function(e) { slog("ERROR", paste("  Failed save prediction raster:", e$message)); FALSE })
        if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
      } else { slog("WARN", paste("  Prediction failed. Output type:", class(prediction_output))); prediction_errors <- prediction_errors + 1 }
      rm(prediction_output); gc()
    }
  } else {
    prediction_errors <- 1 # Only one prediction attempted
    msg <- paste0(log_prefix, "Skipping prediction: Final model was not available."); slog("ERROR", msg)
    status <- if(load_existing_model) "error_loading_model" else "error_training"
    return(list(status = status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
  }
  
  # --- Prepare return status ---
  final_status <- if(predictions_made == 1) "success" else "error_prediction_all"
  status_message <- paste0(log_prefix, "Finished Biotic-Only. Occurrences (thin/clean): ", occurrence_count_after_thinning, ". Prediction status: ", ifelse(predictions_made==1, "Success", "Failed/Skipped"), ".")
  if(final_status == "success") slog("INFO", status_message) else slog("ERROR", status_message)
  
  rm(predictor_stack, occs_coords); gc()
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm_biotic_only function

# --- 9. Setup Parallel Backend & Run ---
if (config$use_parallel && config$num_cores > 1) {
  log4r::info(logger, paste("Setting up parallel backend with", config$num_cores, "cores (multisession)."))
  future::plan(future::multisession, workers = config$num_cores)
} else {
  log4r::info(logger, "Running sequentially.")
  future::plan(future::sequential)
}

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

log4r::info(logger, paste("Starting Biotic-Only SDM processing for", nrow(species_df), "anemonefish..."))

results_list <- progressr::with_progress({
  furrr::future_map(1:nrow(species_df), ~{
    process_species_sdm_biotic_only(
      species_row = species_df[.x, ],
      config = config,
      all_host_predictions = all_host_predictions, # Pass pre-loaded host preds
      associations_df = associations_df,           # Pass associations
      group_pred_dir = group_pred_dir,
      group_results_dir = group_results_dir,
      group_models_dir = group_models_dir,
      species_log_dir = config$species_log_dir,
      predictor_type_suffix = predictor_type_suffix,
      occurrence_dir = occurrence_dir
    )
  }, .options = furrr_options(seed = TRUE))
})

log4r::info(logger, "Parallel/sequential processing complete.")

# --- 10. Process Results ---
# (Identical summary logic as in 06a/06b)
occurrence_counts <- list()
success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0

log4r::info(logger, "--- Processing Results Summary (Anemonefish Biotic-Only Models) ---")
for (res in results_list) {
  if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "Received NULL result from a species process."); next }
  log_level_func <- switch(res$status, success = log4r::info, skipped_occurrences = log4r::warn, log4r::error) # Simplified level mapping
  log_level_func(logger, res$message)
  occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
  if (res$status == "success") success_count <- success_count + 1
  else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
  else error_count <- error_count + 1
}

log4r::info(logger, paste("--- Overall Summary (Anemonefish Biotic-Only Models) ---"))
log4r::info(logger, paste("Total Species Targets:", length(results_list)))
log4r::info(logger, paste("Successfully Trained & Predicted:", success_count))
log4r::info(logger, paste("Skipped (Occurrences/No Hosts):", skipped_count))
log4r::info(logger, paste("Errors during processing:", error_count))

if (length(occurrence_counts) > 0) {
  occ_count_df <- data.frame( Species = names(occurrence_counts), OccurrenceCountAfterThinning = unlist(occurrence_counts)) %>% dplyr::arrange(Species)
  occ_count_df_print <- occ_count_df; occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
  occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, ".csv"))
  tryCatch({ readr::write_csv(occ_count_df, occ_count_file); log4r::info(logger, paste("Occurrence counts saved to:", occ_count_file))}, error = function(e) { log4r::error(logger, paste("Failed to save occurrence counts:", e$message)) })
  cat("\n--- Final Occurrence Counts (Anemonefish, after cleaning/thinning) ---\n"); print(occ_count_df_print)
} else { log4r::warn(logger, "No occurrence counts were recorded.") }

# Shut down parallel workers
future::plan(future::sequential)

log4r::info(logger, "--- Script 06c finished. ---")
#-------------------------------------------------------------------------------