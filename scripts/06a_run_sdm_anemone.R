# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using either PCA or VIF-selected Env
# Predictors (SDMtune Workflow) with Parallel Processing, Logging, and Progress
# Saves Predictions, CV Results, and VarImp to paper-like structure.
# Includes species-specific detailed log files.
# v2: Uses explicit arguments in future_map, simplified globals.
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (v2 Explicit Args) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr, log4r, future, furrr, progressr)

# Source helpers
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R")) # Ensure helpers are updated

# --- 2. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path,
                       log_level = config$log_level,
                       append = config$log_append,
                       log_to_console = config$log_to_console,
                       console_level = config$log_console_level)

log4r::info(logger, "--- Starting Script 06a: Run Standard Anemone SDMs ---")

# --- 3. Define Group Specifics & Predictor Type ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors
predictor_type_suffix <- ifelse(use_pca, "_pca", "_vif") # Suffix for intermediate dirs
# Suffix for paper CV/VI files (env-only models don't have a special suffix in the paper)
paper_file_suffix <- "" # Empty for anemone env-only

log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---"))

# --- 4. Load Predictor Information ---
predictor_paths_or_list <- NULL
if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (!file.exists(pca_paths_rds)) { log4r::fatal(logger, "PCA paths missing."); stop("PCA paths missing.") }
  predictor_paths_or_list <- tryCatch(readRDS(pca_paths_rds), error=function(e){log4r::fatal(logger,"Failed load PCA paths.");NULL})
  if (!is.list(predictor_paths_or_list) || length(predictor_paths_or_list) == 0) { log4r::fatal(logger, "PCA paths list invalid."); stop("PCA paths list invalid.") }
  log4r::info(logger, paste("Loaded PCA raster paths for scenarios:", paste(names(predictor_paths_or_list), collapse=", ")))
} else {
  predictor_paths_or_list <- config$final_vars_vif_anemone
  if(is.null(predictor_paths_or_list) || length(predictor_paths_or_list) < 2) { log4r::fatal(logger, "VIF vars missing."); stop("VIF vars missing.") }
  log4r::info(logger, paste("Using VIF-selected core variables:", paste(predictor_paths_or_list, collapse=", ")))
}
# Determine scenarios to process based on available predictor paths/config
scenarios_to_process <- if(use_pca) names(predictor_paths_or_list) else config$env_scenarios
if(!"current" %in% scenarios_to_process && config$force_rerun$run_standard_sdms) {
  log4r::warn(logger, "The 'current' scenario is needed for tuning but not found in PCA paths or env_scenarios. Tuning might fail if models don't exist.")
}


# --- 5. Define Output Paths (Both Structures) ---
# Intermediate (YOUR structure)
group_results_dir_intermediate <- file.path(config$results_dir, paste0(group_name, predictor_type_suffix))
group_models_dir_intermediate <- file.path(config$models_dir, paste0(group_name, predictor_type_suffix))
species_log_dir <- config$species_log_dir
# Paper Structure (from config)
paper_cv_dir_base <- config$paper_results_base_dir # Env-only CV -> base dir
paper_vi_dir_base <- config$paper_vi_base_dir      # Env-only VI -> base dir
paper_pred_base_dir <- config$paper_pred_base_dir
paper_pred_future_base_dir <- config$paper_pred_future_base_dir

# Create directories (only needs to be done once)
dir.create(group_results_dir_intermediate, recursive = TRUE, showWarnings = FALSE)
dir.create(group_models_dir_intermediate, recursive = TRUE, showWarnings = FALSE)
dir.create(species_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_cv_dir_base, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_vi_dir_base, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_pred_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_pred_future_base_dir, recursive = TRUE, showWarnings = FALSE)
for(ssp in config$original_future_scenarios) { # Ensure future subdirs exist
  dir.create(file.path(paper_pred_future_base_dir, ssp), showWarnings=FALSE, recursive=TRUE)
}

log4r::debug(logger, paste("Intermediate output dirs:", group_results_dir_intermediate, group_models_dir_intermediate))
log4r::debug(logger, paste("Paper output dirs:", paper_cv_dir_base, paper_vi_dir_base, paper_pred_base_dir))


# --- 6. Load Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE) }, error = function(e) { log4r::fatal(logger, "Failed load species list."); stop("Species list loading failed.") })
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))

# --- 7. Define Function to Process Single Species (Corrected Arguments) ---
process_species_sdm <- function(species_row, config, predictor_paths_or_list,
                                paper_cv_dir, paper_vi_dir, paper_pred_base_dir, paper_pred_future_base_dir, # Paper output dirs
                                intermediate_results_dir, intermediate_models_dir, # Intermediate dirs
                                species_log_dir,
                                predictor_type_suffix, paper_file_suffix, # Suffixes
                                use_pca, occurrence_dir, group_name, # Added group_name back
                                tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  log_prefix <- paste0("[", species_name, "] ")
  species_log_file <- file.path(species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), log_prefix, paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
  slog("INFO", "--- Starting processing (Standard Env Model v2) ---")
  
  # Define PAPER Structure Output Paths (using passed args)
  paper_cv_results_file <- file.path(paper_cv_dir, paste0("CV_Results_", species_name_sanitized, paper_file_suffix, ".csv")) # Suffix might be empty
  paper_vi_results_file <- file.path(paper_vi_dir, paste0("vi_", species_name_sanitized, paper_file_suffix, ".csv")) # Suffix might be empty
  
  # Intermediate Paths (using passed args)
  final_model_file_intermediate <- file.path(intermediate_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  tuning_results_file_intermediate <- file.path(intermediate_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  slog("DEBUG", "Paper CV Path:", paper_cv_results_file)
  slog("DEBUG", "Paper VI Path:", paper_vi_results_file)
  slog("DEBUG", "Intermediate Model Path:", final_model_file_intermediate)
  slog("DEBUG", "Intermediate Tuning Path:", tuning_results_file_intermediate)
  
  # --- Load Tuning Predictors ---
  slog("DEBUG", "Loading tuning predictors for scenario:", tuning_scenario)
  tuning_predictor_stack <- NULL
  if(use_pca){
    tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
    if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) { msg <- paste0(log_prefix, "Skipping: PCA stack missing for tuning."); slog("ERROR", msg); return(list(status="error", species=species_name, msg=msg)) }
    tuning_predictor_stack <- tryCatch(terra::rast(tuning_predictor_path), error=function(e) {slog("ERROR", "Failed load tuning PCA stack:", e$message); NULL})
  } else {
    current_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, tuning_scenario, config)
    if(length(current_vif_vars) < 2) { msg <- paste0(log_prefix, "Skipping: Failed generate VIF vars for tuning."); slog("ERROR", msg); return(list(status="error", species=species_name, msg=msg)) }
    tuning_predictor_stack <- load_selected_env_data(tuning_scenario, current_vif_vars, config)
  }
  if(is.null(tuning_predictor_stack)) { msg <- paste0(log_prefix, "Skipping: Failed load tuning predictor stack."); slog("ERROR", msg); return(list(status="error", species=species_name, msg=msg)) }
  slog("DEBUG", "Tuning predictor stack loaded.")
  
  # --- Load/Clean/Thin Occurrences ---
  config_for_occ_load <- config; config_for_occ_load$predictor_stack_for_thinning <- tuning_predictor_stack
  slog("DEBUG", "Loading/cleaning/thinning occurrences.")
  occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL)
  if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) { msg <- paste0(log_prefix, "Skipping: Insufficient occs (", occ_data_result$count %||% 0, ")."); slog("WARN", msg); return(list(status="skipped", species=species_name, msg=msg, count=occ_data_result$count %||% 0)) }
  occs_coords <- occ_data_result$coords; occurrence_count_after_thinning <- occ_data_result$count
  slog("INFO", "Occ count after clean/thin:", occurrence_count_after_thinning)
  
  # --- Check if Final Model Exists ---
  final_model <- NULL; tuning_results_object <- NULL
  if (!config$force_rerun$run_standard_sdms && file.exists(final_model_file_intermediate)) {
    slog("INFO", "Skipping tuning/training: Intermediate model exists.")
    load_existing_model <- TRUE
    final_model <- tryCatch(readRDS(final_model_file_intermediate), error=function(e) { slog("ERROR", "Failed load existing model:", e$message); NULL})
    if(is.null(final_model) || !inherits(final_model, "SDMmodel")){ msg <- paste0(log_prefix, "Invalid existing model. Re-running."); slog("WARN", msg); load_existing_model <- FALSE }
    else { if (file.exists(tuning_results_file_intermediate)) { tuning_results_object <- tryCatch(readRDS(tuning_results_file_intermediate), error = function(e) { slog("WARN", "Could not load tuning object for VI."); NULL }) } else { slog("WARN", "Tuning results file missing for VI.") } }
  } else { load_existing_model <- FALSE }
  
  # --- Run Tuning & Training if needed ---
  if (!load_existing_model) {
    slog("INFO", "Running tuning/training...")
    print("FUCK1")
    background_points <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config, logger = NULL, seed = species_aphia_id)
    print("FUCK2")
    if (is.null(background_points)) { msg <- paste0(log_prefix, "Skipping: Failed background."); slog("ERROR", msg); return(list(status="error", species=species_name, msg=msg, count=occurrence_count_after_thinning)) }
    print("FUCK3")
    tuning_output_list <- run_sdm_tuning_kfold(occs_coords, tuning_predictor_stack, background_points, config, logger = NULL, species_name)
    print("FUCK4")
    if (is.null(tuning_output_list)) { msg <- paste0(log_prefix, "Skipping: Tuning failed."); slog("ERROR", msg); return(list(status="error", species=species_name, msg=msg, count=occurrence_count_after_thinning)) }
    print("FUCK5")
    best_hypers <- tuning_output_list$best_hypers; tuning_results_object <- tuning_output_list$tuning_results_object
    print("FUCK6")
    slog("INFO", "Tuning complete. Best hypers:", paste(names(best_hypers), best_hypers[1,], collapse=", "))
    
    # Save CV results to PAPER structure
    if (!is.null(tuning_results_object) && inherits(tuning_results_object, "SDMtune") && !is.null(tuning_results_object@results) && nrow(tuning_results_object@results) > 0) {
      tryCatch({ dir.create(dirname(paper_cv_results_file), recursive = TRUE, showWarnings = FALSE); readr::write_csv(tuning_results_object@results, paper_cv_results_file); slog("DEBUG", "CV results saved (paper):", basename(paper_cv_results_file)) }, error = function(e) { slog("ERROR", "Failed save paper CV results:", e$message) })
    } else { slog("WARN", "Cannot save paper CV results.") }
    
    # Save INTERMEDIATE tuning object
    save_tuning_result <- tryCatch({ saveRDS(tuning_results_object, tuning_results_file_intermediate); slog("DEBUG", "Intermediate Tuning object saved."); NULL }, error = function(e){ slog("ERROR", "Failed save intermediate tuning obj:", e$message); return(list(status="error", species=species_name, msg=paste0(log_prefix, "Failed save tuning intermediate: ", e$message), count=occurrence_count_after_thinning)) })
    if (!is.null(save_tuning_result)) return(save_tuning_result)
    
    final_model <- train_final_sdm(occs_coords, tuning_predictor_stack, background_points, best_hypers, config, logger = NULL, species_name)
    if(is.null(final_model)) { msg <- paste0(log_prefix, "Skipping: Training failed."); slog("ERROR", msg); return(list(status="error", species=species_name, msg=msg, count=occurrence_count_after_thinning)) }
    slog("INFO", "Final model training complete.")
    
    # Save INTERMEDIATE model object
    save_model_result <- tryCatch({ saveRDS(final_model, final_model_file_intermediate); slog("DEBUG", "Intermediate final model saved."); NULL }, error = function(e){ slog("ERROR", "Failed save final model:", e$message); return(list(status="error", species=species_name, msg=paste0(log_prefix, "Failed save final model: ", e$message), count=occurrence_count_after_thinning)) })
    if (!is.null(save_model_result)) return(save_model_result)
    rm(background_points); gc()
  }
  
  # --- Variable Importance Calculation & Saving ---
  if (config$run_variable_importance && !is.null(final_model)) {
    test_swd_data <- NULL
    if (!is.null(tuning_results_object) && inherits(tuning_results_object@models[[1]], "SDMmodelCV")) { test_swd_data <- tuning_results_object@models[[1]]@data }
    else { slog("WARN", "Re-preparing SWD for VI."); if (!is.null(tuning_predictor_stack) && !is.null(occs_coords)) { background_points_vi <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config, logger = NULL, seed = species_aphia_id + 1); if (!is.null(background_points_vi)){ test_swd_data <- tryCatch(SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_points_vi, env = tuning_predictor_stack, verbose=FALSE), error=function(e) NULL); rm(background_points_vi); gc() }} }
    if (!is.null(test_swd_data)) { calculate_and_save_vi(final_model, test_swd_data, config$vi_permutations, paper_vi_results_file, logger = NULL) } else { slog("WARN", "Skipping VI: Couldn't get test SWD.") }
  } else if (config$run_variable_importance) { slog("WARN", "Skipping VI: Final model missing.") }
  rm(tuning_predictor_stack, tuning_results_object); gc()
  
  # --- Prediction Loop ---
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- if(use_pca) names(predictor_paths_or_list) else config$env_scenarios
  slog("INFO", paste("Starting predictions for", length(scenarios_to_predict), "scenarios."))
  if (!is.null(final_model)) {
    for (pred_scenario in scenarios_to_predict) {
      slog("DEBUG", paste("  Predicting scenario:", pred_scenario))
      # Construct PAPER Structure Output Path
      if (pred_scenario == "current") { target_pred_dir <- paper_pred_base_dir; target_pred_filename <- paste0("mean_pred_", species_name_sanitized, paper_file_suffix, ".tif") } else { ssp_match <- stringr::str_match(pred_scenario, "(ssp\\d{3})_(\\d{4})"); if(is.na(ssp_match[1,1])) { slog("WARN", "Cannot parse SSP/Year."); prediction_errors <- prediction_errors + 1; next }; ssp_code <- ssp_match[1, 2]; time_tag_clean <- ifelse(ssp_match[1, 3] == "2050", "dec50", "dec100"); target_pred_dir <- file.path(paper_pred_future_base_dir, ssp_code); target_pred_filename <- paste0("mean_pred_", species_name_sanitized, "_", ssp_code, "_", time_tag_clean, paper_file_suffix, ".tif") }
      pred_file_paper <- file.path(target_pred_dir, target_pred_filename)
      if (!config$force_rerun$run_standard_sdms && file.exists(pred_file_paper)) { slog("DEBUG", "  Prediction exists. Skipping."); next }
      # Load predictor stack
      pred_predictor_stack <- NULL
      if(use_pca) { pred_predictor_path <- predictor_paths_or_list[[pred_scenario]]; if (is.null(pred_predictor_path) || !file.exists(pred_predictor_path)) { slog("WARN", "PCA stack missing"); prediction_errors <- prediction_errors + 1; next }; pred_predictor_stack <- tryCatch(terra::rast(pred_predictor_path), error = function(e) {slog("WARN", "Error loading PCA stack"); NULL})
      } else { scenario_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, pred_scenario, config); if(length(scenario_vif_vars) < 1) { slog("WARN", "No VIF vars"); prediction_errors <- prediction_errors + 1; next }; pred_predictor_stack <- load_selected_env_data(pred_scenario, scenario_vif_vars, config) }
      if(is.null(pred_predictor_stack)) { slog("WARN", "Failed load predictor stack."); prediction_errors <- prediction_errors + 1; next }
      # Predict
      prediction_output <- predict_sdm_suitability(final_model, pred_predictor_stack, config, logger = NULL)
      # Save
      if (inherits(prediction_output, "SpatRaster")) { save_success <- tryCatch({ dir.create(dirname(pred_file_paper), recursive = TRUE, showWarnings = FALSE); terra::writeRaster(prediction_output, filename = pred_file_paper, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES")); slog("DEBUG", paste("  Prediction saved:", basename(pred_file_paper))); TRUE }, error = function(e) { slog("ERROR", paste("  Failed save prediction:", e$message)); FALSE }); if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1; rm(prediction_output); gc()
      } else { slog("WARN", paste("  Prediction failed:", prediction_output)); prediction_errors <- prediction_errors + 1 }
      rm(pred_predictor_stack); gc()
    }
  } else { prediction_errors <- length(scenarios_to_predict); msg <- paste0(log_prefix, "Skipping preds: Final model unavailable."); slog("ERROR", msg); status <- if(load_existing_model) "error_loading_model" else "error_training"; return(list(status=status, species=species_name, count=occurrence_count_after_thinning, message=msg)) }
  
  # --- Prepare return status ---
  final_status <- "success"; status_message <- paste0(log_prefix, "Finished. Occs:", occurrence_count_after_thinning, ". Preds attempted:", length(scenarios_to_predict), ". Made:", predictions_made, ". Errors/Skipped:", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) { final_status <- "error_prediction_all"; slog("ERROR", status_message) } else if (prediction_errors > 0) { final_status <- "success_with_pred_errors"; slog("WARN", status_message) } else { slog("INFO", status_message) }
  rm(occs_coords); gc()
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm

# --- 8. Setup Parallel Backend & Run ---
if (config$use_parallel && config$num_cores > 1) { log4r::info(logger, paste("Setting up parallel with", config$num_cores, "cores.")); gc(full=TRUE); future::plan(future::multisession, workers = config$num_cores, gc = TRUE) } else { log4r::info(logger, "Running sequentially."); future::plan(future::sequential) }
progressr::handlers(global = TRUE); progressr::handlers("txtprogressbar")
log4r::info(logger, paste("Starting Standard Anemone SDM processing for", nrow(species_df), "species..."))

# Updated future_map call with explicit globals
results_list <- progressr::with_progress({
  furrr::future_map(1:nrow(species_df), ~{
    process_species_sdm(
      species_row = species_df[.x, ],
      config = config,
      predictor_paths_or_list = predictor_paths_or_list,
      paper_cv_dir = config$paper_results_base_dir, # Env-only CV -> base dir
      paper_vi_dir = config$paper_vi_base_dir,      # Env-only VI -> base dir
      paper_pred_base_dir = config$paper_pred_base_dir,
      paper_pred_future_base_dir = config$paper_pred_future_base_dir,
      intermediate_results_dir = group_results_dir_intermediate, # Use intermediate path
      intermediate_models_dir = group_models_dir_intermediate,   # Use intermediate path
      species_log_dir = config$species_log_dir,
      predictor_type_suffix = predictor_type_suffix,
      paper_file_suffix = paper_file_suffix, # Pass empty suffix for env-only
      use_pca = use_pca,
      occurrence_dir = occurrence_dir,
      group_name = group_name # Pass group_name explicitly
    )
  },
  .options = furrr_options(
    seed = TRUE,
    globals = c("config",                     # Main config list
                "species_df",                 # Data frame being iterated
                "predictor_type_suffix",      # Suffix for intermediate files
                "paper_file_suffix",          # <<< ADDED paper_file_suffix HERE
                "group_name",                 # Name of the group (anemone/anemonefish)
                "intermediate_models_dir",    # Path for intermediate models
                "intermediate_results_dir",   # Path for intermediate tuning results
                # --- Helper Functions ---
                "process_species_sdm",        # The main function being called
                "load_clean_individual_occ_coords",
                "generate_sdm_background",
                "run_sdm_tuning_kfold",
                "train_final_sdm",
                "predict_sdm_suitability",
                "calculate_and_save_vi",
                "load_selected_env_data",      # Needed if use_pca = FALSE
                "generate_scenario_variable_list", # Needed if use_pca = FALSE
                "group_models_dir_intermediate",
                "group_results_dir_intermediate",
                "use_pca",
                "predictor_paths_or_list",
                "occurrence_dir"
    )
  )
  )
})

log4r::info(logger, "Parallel/sequential processing complete.")

# --- 9. Process Results ---
# (Identical summary logic as before)
occurrence_counts <- list(); success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0
log4r::info(logger, "--- Processing Results Summary (Standard Anemone Models) ---")
for (res in results_list) {
  if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "NULL result."); next }
  log_level_func <- switch(res$status, success=log4r::info, success_with_pred_errors=log4r::warn, skipped_occurrences=log4r::warn, log4r::error)
  log_level_func(logger, res$message)
  occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
  if (res$status == "success") success_count <- success_count + 1
  else if (res$status == "success_with_pred_errors") partial_success_count <- partial_success_count + 1
  else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
  else error_count <- error_count + 1
}
log4r::info(logger, paste("--- Overall Summary (Standard Anemone Models) ---"))
log4r::info(logger, paste("Total Species:", length(results_list))); log4r::info(logger, paste("Fully Successful:", success_count)); log4r::info(logger, paste("Partial Success:", partial_success_count)); log4r::info(logger, paste("Skipped:", skipped_count)); log4r::info(logger, paste("Errors:", error_count))

if (length(occurrence_counts) > 0) {
  occ_count_df <- data.frame( Species = names(occurrence_counts), OccurrenceCountAfterThinning = unlist(occurrence_counts)) %>% dplyr::arrange(Species)
  occ_count_df_print <- occ_count_df; occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
  occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, ".csv"))
  tryCatch({ readr::write_csv(occ_count_df, occ_count_file); log4r::info(logger, paste("Occurrence counts saved:", occ_count_file))}, error = function(e) { log4r::error(logger, paste("Failed save counts:", e$message)) })
  cat("\n--- Final Occurrence Counts (Anemone, after cleaning/thinning) ---\n"); print(occ_count_df_print)
} else { log4r::warn(logger, "No occurrence counts recorded.") }

future::plan(future::sequential); gc(full=TRUE)
log4r::info(logger, "--- Script 06a finished. ---")
#-------------------------------------------------------------------------------