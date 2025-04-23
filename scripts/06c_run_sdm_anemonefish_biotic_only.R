# scripts/06c_run_sdm_anemonefish_biotic_pc1.R
#-------------------------------------------------------------------------------
# Run SDMs for Anemonefish Species using PC1 + Max Host Suitability
# (SDMtune Workflow). Includes Parallel Processing, Logging, Progress,
# and Species-Specific Logs. Saves outputs to target structure.
# Assumes 06a produced host predictions (_pca suffix) for ALL scenarios.
# v6 - FINAL VI Call matching working 06a.
#-------------------------------------------------------------------------------
cat("--- Running Script 06c: Run Anemonefish Biotic+PC1 SDMs (v6 - Final VI Call) ---\n")

# --- 1. Setup: Load Config FIRST ---
if (file.exists("scripts/config.R")) {
  source("scripts/config.R")
  if (!exists("config") || !is.list(config)) {
    stop("FATAL: 'config' list object not found or invalid after sourcing scripts/config.R")
  }
} else {
  stop("FATAL: Configuration file 'scripts/config.R' not found.")
}

# --- 2. Load Required Packages ---
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr, log4r, future, furrr, progressr)

# --- 3. Source Helper Functions ---
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R")) # Needs v8 of calculate_and_save_vi

# --- 4. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path,
                       log_level = config$log_level,
                       append = config$log_append,
                       log_to_console = config$log_to_console,
                       console_level = config$log_console_level)

log4r::info(logger, "--- Starting Script 06c: Run Anemonefish Biotic+PC1 SDMs (v6) ---")

# --- 5. Define Group Specifics & Predictor Type ---
group_name <- "anemonefish"
species_list_file <- config$anemonefish_species_list_file
occurrence_dir <- config$anemonefish_occurrence_dir
host_group_name <- "anemone"

if(!config$use_pca_predictors) { log4r::fatal(logger, "Script 06c requires config$use_pca_predictors=TRUE."); stop("This script requires PCA env predictors.") }
predictor_type_suffix <- "_biotic_pc1"
host_predictor_type_suffix <- "_pca"

log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors: Biotic (Max Host Suitability", host_predictor_type_suffix, ") + PC1 ---"))

# --- 6. Load Environmental Predictor Paths (PCA) ---
pca_paths_rds <- config$pca_raster_paths_rds_path
if (is.null(pca_paths_rds) || !is.character(pca_paths_rds) || nchar(pca_paths_rds) == 0 || !file.exists(pca_paths_rds)) { log4r::fatal(logger, "PCA paths RDS missing/invalid."); stop("PCA paths RDS invalid.") }
env_predictor_paths_list <- tryCatch(readRDS(pca_paths_rds), error = function(e) { log4r::fatal(logger, paste("Failed load PCA paths RDS:", e$message)); NULL })
if (!is.list(env_predictor_paths_list) || length(env_predictor_paths_list) == 0) { log4r::fatal(logger, "PCA paths list empty/invalid."); stop("PCA paths list invalid.") }
scenarios_to_process <- names(env_predictor_paths_list)
log4r::info(logger, paste("Loaded Env PCA paths for scenarios:", paste(scenarios_to_process, collapse=", ")))

# --- 7. Load Processed Host Association Data ---
association_file_path <- config$anemone_fish_association_file
if(!file.exists(association_file_path)) { log4r::fatal(logger, "Association file missing."); stop("Association file missing.") }
associations_df <- tryCatch(readr::read_csv(association_file_path, show_col_types = FALSE), error = function(e){ log4r::fatal(logger, paste("Failed load associations:", e$message)); NULL })
if(is.null(associations_df)) stop("Failed to load associations file.")
log4r::info(logger, paste("Loaded", nrow(associations_df), "fish-host associations."))

# --- 8. Get File Paths for ALL Host Anemone Predictions ---
log4r::info(logger, "Gathering file paths for ALL potential host anemone predictions...")
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
    host_pred_file_path <- construct_prediction_filename(host_sci_name_sanitized, scenario, host_predictor_type_suffix, config) # Use helper
    if(!is.null(host_pred_file_path) && file.exists(host_pred_file_path)) {
      all_host_prediction_paths[[scenario]][[host_sci_name_sanitized]] <- host_pred_file_path
    } else {
      log4r::warn(logger, paste("Host prediction file missing or path construction failed:", host_pred_file_path %||% paste(host_sci_name_sanitized, scenario, host_predictor_type_suffix)))
      missing_files_scenario = missing_files_scenario + 1
    }
  }
  if(length(all_host_prediction_paths[[scenario]]) == 0 && missing_files_scenario > 0){ log4r::error(logger, paste("FATAL: No host predictions found for ANY species in scenario:", scenario)) }
  else if (missing_files_scenario > 0) { log4r::warn(logger, paste("Found", length(all_host_prediction_paths[[scenario]]), "host files, but", missing_files_scenario, "missing for:", scenario)) }
  else { log4r::info(logger, paste("Found all", length(all_host_prediction_paths[[scenario]]), "host files for:", scenario)) }
}

# --- 9. Create Intermediate Output Dirs ---
base_intermediate_model_path <- config$models_dir_intermediate
base_intermediate_results_path <- config$results_dir_intermediate
base_species_log_path <- config$species_log_dir
if (any(sapply(list(base_intermediate_model_path, base_intermediate_results_path, base_species_log_path, group_name, predictor_type_suffix), function(x) is.null(x) || !is.character(x) || nchar(x) == 0))) { stop("FATAL: Base path/group/suffix invalid.") }
intermediate_models_dir <- file.path(base_intermediate_model_path, paste0(group_name, predictor_type_suffix))
intermediate_results_dir <- file.path(base_intermediate_results_path, paste0(group_name, predictor_type_suffix))
tryCatch({
  dir.create(intermediate_models_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(intermediate_results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_species_log_path, recursive = TRUE, showWarnings = FALSE)
  log4r::debug(logger, paste("Intermediate dirs created/checked:", group_name, predictor_type_suffix))
}, error = function(e) { log4r::fatal(logger, paste("Failed create intermediate dirs:", e$message)); stop("Dir creation failed.") })

# --- 10. Load Anemonefish Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE) }, error = function(e) { log4r::fatal(logger, paste("Failed load species list:", e$message)); stop("Species list failed.") })
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))

# --- 11. Define Function to Process Single Species (Biotic + PC1) ---
# Uses the v8 VI call logic (final_model + training_swd)
process_species_sdm_biotic_pc1 <- function(species_row, config, env_predictor_paths_list, all_host_prediction_paths, associations_df,
                                           group_name, predictor_type_suffix, occurrence_dir,
                                           species_log_dir, tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  species_log_file <- file.path(species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), paste0("[", species_name, "]"), paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
  slog("INFO", "--- Starting Biotic+PC1 processing (v6) ---")
  
  intermediate_models_dir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
  intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix))
  final_model_file <- file.path(intermediate_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  tuning_rds_file <- file.path(intermediate_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, "_object.rds"))
  
  associated_hosts <- associations_df %>% filter(AnemonefishScientificName == species_name) %>% pull(AssociatedAnemoneScientificName)
  associated_hosts_sanitized <- gsub(" ", "_", associated_hosts)
  if (length(associated_hosts_sanitized) == 0) { msg <- paste0("Skipping: No associated hosts found."); slog("WARN", msg); return(list(status = "skipped_no_hosts", species = species_name, occurrence_count = NA, message = msg)) }
  slog("DEBUG", paste("Associated hosts:", paste(associated_hosts_sanitized, collapse=", ")))
  
  slog("DEBUG", "Preparing Biotic+PC1 predictor stack for tuning scenario:", tuning_scenario)
  tuning_env_path <- env_predictor_paths_list[[tuning_scenario]]
  if (is.null(tuning_env_path) || !file.exists(tuning_env_path)) { msg <- paste0("Skipping: Env PCA stack path missing for tuning."); slog("ERROR", msg); return(list(status = "error_tuning_env_pca_path", species = species_name, occurrence_count = NA, message = msg)) }
  pc1_tuning_layer <- tryCatch(terra::rast(tuning_env_path)[[1]], error = function(e) { slog("ERROR", "Failed load tuning PC1:", e$message); NULL })
  if (is.null(pc1_tuning_layer)) { msg <- paste0("Skipping: Failed load tuning PC1 layer."); slog("ERROR", msg); return(list(status = "error_tuning_env_pc1_load", species = species_name, occurrence_count = NA, message = msg)) }
  reference_geom_tuning <- pc1_tuning_layer
  
  host_paths_for_tuning <- all_host_prediction_paths[[tuning_scenario]][names(all_host_prediction_paths[[tuning_scenario]]) %in% associated_hosts_sanitized]
  if (length(host_paths_for_tuning) == 0) { msg <- paste0("Skipping: No host prediction paths found for tuning."); slog("ERROR", msg); return(list(status = "error_tuning_biotic_paths", species = species_name, occurrence_count = NA, message = msg)) }
  if (length(host_paths_for_tuning) < length(associated_hosts_sanitized)) { slog("WARN", "Missing paths for some hosts in tuning.") }
  
  host_rasters_tuning_list <- list()
  for(host_name in names(host_paths_for_tuning)){
    host_path <- host_paths_for_tuning[[host_name]]
    host_rast <- tryCatch(terra::rast(host_path), error=function(e){slog("WARN", paste("Failed load host", basename(host_path),":", e$message)); NULL})
    if(!is.null(host_rast)){
      if (!terra::compareGeom(reference_geom_tuning, host_rast, stopOnError=FALSE, res=TRUE)) { slog("DEBUG", paste("Resampling host", basename(host_path))); host_rast <- tryCatch(terra::resample(host_rast, reference_geom_tuning, method="bilinear"), error = function(e) { slog("ERROR", paste("Failed resample host", basename(host_path),":", e$message)); NULL }) }
      if(!is.null(host_rast)){ names(host_rast) <- host_name; host_rasters_tuning_list[[host_name]] <- host_rast }
    }
  }
  if (length(host_rasters_tuning_list) == 0) { msg <- paste0("Skipping: Failed load/resample ANY hosts for tuning."); slog("ERROR", msg); return(list(status = "error_tuning_biotic_load", species = species_name, occurrence_count = NA, message = msg)) }
  host_stack_tuning <- tryCatch(terra::rast(host_rasters_tuning_list), error = function(e) { slog("ERROR", paste("Failed stack host predictors:", e$message)); NULL })
  if (is.null(host_stack_tuning) || terra::nlyr(host_stack_tuning) == 0) { msg <- paste0("Skipping: Failed create valid host stack."); slog("ERROR", msg); return(list(status = "error_tuning_biotic_stack", species = species_name, occurrence_count = NA, message = msg)) }
  max_host_suitability_tuning <- tryCatch(terra::app(host_stack_tuning, fun = "max", na.rm = TRUE), error = function(e) { slog("ERROR", paste("Failed calc max host:", e$message)); NULL })
  if(is.null(max_host_suitability_tuning)) { msg <- paste0("Skipping: Failed calc max host."); slog("ERROR", msg); return(list(status = "error_tuning_biotic_max", species = species_name, occurrence_count = NA, message = msg)) }
  names(max_host_suitability_tuning) <- "host_suitability_max"
  slog("DEBUG", "Max host suitability layer created for tuning.")
  
  tuning_predictor_stack <- tryCatch(c(pc1_tuning_layer, max_host_suitability_tuning), error = function(e) { slog("ERROR", paste("Failed combine PC1+Biotic:", e$message)); NULL })
  if(is.null(tuning_predictor_stack)) { msg <- paste0("Skipping: Failed combine stacks."); slog("ERROR", msg); return(list(status = "error_tuning_combine_stack", species = species_name, occurrence_count = NA, message = msg)) }
  names(tuning_predictor_stack) <- c("PC1", "host_suitability_max")
  slog("INFO", paste("Tuning predictor stack created:", paste(names(tuning_predictor_stack), collapse=", ")))
  
  raster_crs_terra <- terra::crs(tuning_predictor_stack)
  
  config_for_occ_load <- config; config_for_occ_load$predictor_stack_for_thinning <- tuning_predictor_stack
  slog("DEBUG", "Loading/cleaning/thinning occurrences.")
  occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL, species_log_file=species_log_file)
  if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) { msg <- paste0("Skipping: Insufficient occurrences (", occ_data_result$count %||% 0, ")."); slog("WARN", msg); return(list(status = "skipped_occurrences", species = species_name, occurrence_count = occ_data_result$count %||% 0, message = msg)) }
  occs_coords <- occ_data_result$coords; occurrence_count_after_thinning <- occ_data_result$count
  occs_sf_clean <- sf::st_as_sf(as.data.frame(occs_coords), coords=c("decimalLongitude","decimalLatitude"), crs=config$occurrence_crs) # sf for spatial ops
  if(sf::st_crs(occs_sf_clean) != sf::st_crs(raster_crs_terra)){ occs_sf_clean <- sf::st_transform(occs_sf_clean, crs=sf::st_crs(raster_crs_terra))}
  slog("INFO", "Occurrence count after clean/thin:", occurrence_count_after_thinning)
  
  final_model <- NULL; tuning_output <- NULL; load_existing_model <- FALSE; full_swd_data <- NULL
  
  if (!config$force_rerun$run_biotic_sdms && file.exists(final_model_file)) { # Check biotic flag
    slog("INFO", "Final model file exists. Attempting to load...")
    final_model <- tryCatch(readRDS(final_model_file), error=function(e) { slog("ERROR", "Failed load existing model:", e$message); NULL})
    if(!is.null(final_model) && inherits(final_model, "SDMmodel")){
      slog("INFO", "Successfully loaded existing final model. Skipping tuning/training.")
      load_existing_model <- TRUE
      if(file.exists(tuning_rds_file)) { tuning_output <- tryCatch(readRDS(tuning_rds_file), error=function(e){slog("WARN","Could not load tuning RDS for VI."); NULL}) }
      if(is.null(tuning_output)) { slog("WARN", "Tuning RDS file not found/loaded, VI will use final model object."); tuning_output <- final_model }
      background_points_for_vi <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack, config, logger=NULL, species_log_file=species_log_file, seed_offset = species_aphia_id)
      if(!is.null(background_points_for_vi)){ full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_points_for_vi, env = tuning_predictor_stack, verbose = FALSE) }, error = function(e) { slog("WARN", "Failed prepareSWD for VI:", e$message); NULL }) } else { slog("WARN", "Failed background gen for VI.") }
      rm(background_points_for_vi); gc()
    } else { msg <- paste0("Existing final model file invalid/failed load. Will re-run."); slog("WARN", msg); load_existing_model <- FALSE }
  }
  
  if (!load_existing_model) {
    slog("INFO", "Final model not found/invalid or rerun forced. Proceeding with tuning/training.")
    slog("DEBUG", "Generating background points.")
    background_points <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack, config, logger=NULL, species_log_file=species_log_file, seed_offset = species_aphia_id)
    if (is.null(background_points)) { msg <- paste0("Skipping: Failed background generation."); slog("ERROR", msg); return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg)) }
    slog("DEBUG", "Background points generated.")
    
    full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_points, env = tuning_predictor_stack, verbose = FALSE) }, error = function(e) { slog("ERROR", paste("Failed prepareSWD:", e$message)); NULL })
    if (is.null(full_swd_data)) { rm(background_points); gc(); return(list(status = "error_swd_preparation", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("log_prefix", "SWD prep failed."))) }
    slog("DEBUG", "SWD object prepared.")
    
    slog("INFO", "Starting hyperparameter tuning.")
    tuning_output <- run_sdm_tuning_scv(occs_coords, tuning_predictor_stack, background_points, config, logger=NULL, species_name, species_log_file=species_log_file)
    if (is.null(tuning_output) || is.null(attr(tuning_output, "best_hypers"))) { msg <- paste0("Skipping: Tuning failed."); slog("ERROR", msg); rm(background_points, full_swd_data); gc(); return(list(status = "error_tuning", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg)) }
    best_hypers <- attr(tuning_output, "best_hypers")
    if(!save_tuning_results(tuning_output, species_name_sanitized, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)) { rm(background_points, full_swd_data); gc(); return(list(status = "error_saving_tuning_results", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("Failed save tuning results."))) }
    
    slog("INFO", "Starting final model training.")
    final_model <- train_final_sdm(occs_coords, tuning_predictor_stack, background_points, best_hypers, config, logger=NULL, species_name, species_log_file=species_log_file)
    if(is.null(final_model)){ msg <- paste0("Skipping: Training failed."); slog("ERROR", msg); rm(background_points, full_swd_data); gc(); return(list(status = "error_training", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg)) }
    slog("INFO", "Final model training complete.")
    if(!save_final_model(final_model, species_name_sanitized, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)) { rm(background_points, full_swd_data); gc(); return(list(status = "error_saving_model", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("Failed save final model."))) }
    
    rm(background_points); gc()
  } # End if !load_existing_model
  
  # --- Variable Importance (Using the same logic as final 06a) ---
  if (!is.null(final_model)) {
    if(!is.null(full_swd_data)){ # Check SWD data
      slog("INFO", "Calculating and saving variable importance...")
      vi_success <- calculate_and_save_vi(
        final_model = final_model,             # Pass the final SDMmodel
        training_swd = full_swd_data,        # Pass the full training SWD data
        species_name_sanitized = species_name_sanitized,
        group_name = group_name,
        predictor_type_suffix = predictor_type_suffix,
        config = config,
        logger = NULL,
        species_log_file = species_log_file)
      if(!vi_success) slog("WARN", "Variable importance calculation/saving failed.")
    } else {
      slog("WARN", "Skipping variable importance: required training SWD data is missing.")
    }
  } else {
    slog("WARN", "Skipping variable importance: final model object is unavailable.")
  }
  
  # --- Prediction Loop ---
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- names(env_predictor_paths_list) # Only scenarios where PCA paths exist
  slog("INFO", paste("Starting predictions for", length(scenarios_to_predict), "scenarios."))
  
  if (!is.null(final_model)) {
    for (pred_scenario in scenarios_to_predict) {
      slog("DEBUG", paste("  Predicting scenario:", pred_scenario))
      target_pred_file_path <- construct_prediction_filename(species_name_sanitized, pred_scenario, predictor_type_suffix, config)
      if (is.null(target_pred_file_path)) { slog("WARN", "  Could not construct prediction filename. Skipping."); prediction_errors <- prediction_errors + 1; next }
      if (!config$force_rerun$run_biotic_sdms && file.exists(target_pred_file_path)) { slog("DEBUG", "  Prediction exists. Skipping."); next } # Check biotic flag
      
      # 1. Load PC1 layer
      pred_env_path <- env_predictor_paths_list[[pred_scenario]]
      if (is.null(pred_env_path) || !file.exists(pred_env_path)) { slog("WARN", paste("  Env PCA stack missing:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      pc1_pred_layer <- tryCatch(terra::rast(pred_env_path)[[1]], error = function(e) { slog("WARN", paste("  Error load env PC1:", e$message)); NULL })
      if(is.null(pc1_pred_layer)) { prediction_errors <- prediction_errors + 1; next }
      reference_geom_pred <- pc1_pred_layer
      
      # 2. Get host paths
      host_paths_for_pred <- all_host_prediction_paths[[pred_scenario]][names(all_host_prediction_paths[[pred_scenario]]) %in% associated_hosts_sanitized]
      if (length(host_paths_for_pred) == 0) { slog("WARN", paste("  No host paths found:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      
      # 3. Load/Resample/Stack Hosts
      host_rasters_pred_list <- list()
      for(host_name in names(host_paths_for_pred)){
        host_path <- host_paths_for_pred[[host_name]]
        host_rast <- tryCatch(terra::rast(host_path), error=function(e){slog("WARN",paste("Failed load host",basename(host_path),":",e$message));NULL})
        if(!is.null(host_rast)){
          if (!terra::compareGeom(reference_geom_pred, host_rast, stopOnError=FALSE, res=TRUE)) { slog("DEBUG", paste("Resampling host", basename(host_path))); host_rast <- tryCatch(terra::resample(host_rast, reference_geom_pred, method="bilinear"), error=function(e){slog("ERROR",paste("Failed resample",basename(host_path),":",e$message));NULL}) }
          if(!is.null(host_rast)){ names(host_rast) <- host_name; host_rasters_pred_list[[host_name]] <- host_rast }
        }
      }
      if (length(host_rasters_pred_list) == 0) { slog("WARN", paste("  Failed load/resample hosts:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      host_stack_scenario <- tryCatch(terra::rast(host_rasters_pred_list), error = function(e) { slog("WARN", paste("  Error stacking hosts:", e$message)); NULL })
      if (is.null(host_stack_scenario)) { prediction_errors <- prediction_errors + 1; next }
      
      # 4. Calculate Max Host Suitability
      max_host_suitability_scenario <- tryCatch(terra::app(host_stack_scenario, fun = "max", na.rm = TRUE), error = function(e) { slog("WARN", paste("  Error calc max host:", e$message)); NULL })
      if(is.null(max_host_suitability_scenario)) { prediction_errors <- prediction_errors + 1; next }
      names(max_host_suitability_scenario) <- "host_suitability_max"
      
      # 5. Combine PC1 + Biotic for prediction
      prediction_predictor_stack <- tryCatch(c(pc1_pred_layer, max_host_suitability_scenario), error=function(e){slog("WARN", paste("  Error combining PC1+Biotic:", e$message)); NULL})
      if(is.null(prediction_predictor_stack)) { prediction_errors <- prediction_errors + 1; next }
      names(prediction_predictor_stack) <- c("PC1", "host_suitability_max")
      
      # 6. Predict & Save
      prediction_output <- predict_sdm_suitability(final_model, prediction_predictor_stack, config, logger=NULL, species_log_file=species_log_file)
      if (inherits(prediction_output, "SpatRaster")) {
        save_success <- save_sdm_prediction(prediction_output, species_name_sanitized, pred_scenario, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)
        if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
        rm(prediction_output); gc()
      } else { slog("WARN", paste("  Prediction failed:", prediction_output)); prediction_errors <- prediction_errors + 1 }
      rm(pc1_pred_layer, host_stack_scenario, max_host_suitability_scenario, prediction_predictor_stack, host_rasters_pred_list); gc()
    } # End prediction scenario loop
  } else {
    prediction_errors <- length(scenarios_to_predict); msg <- paste0("Skipping predictions: Final model unavailable."); slog("ERROR", msg)
    status <- if(load_existing_model) "error_loading_model" else "error_training"; return(list(status = status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
  }
  
  # --- Prepare return status ---
  final_status <- "success"; status_message <- paste0("Finished Biotic+PC1 SDM. Occs:", occurrence_count_after_thinning, ". Preds attempted:", length(scenarios_to_predict), ". Made:", predictions_made, ". Errors/Skipped:", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) { final_status <- "error_prediction_all"; slog("ERROR", status_message) } else if (prediction_errors > 0) { final_status <- "success_with_pred_errors"; slog("WARN", status_message) } else { slog("INFO", status_message) }
  
  # --- Clean up ---
  if (!is.null(tuning_predictor_stack)) rm(tuning_predictor_stack)
  if (!is.null(occs_coords)) rm(occs_coords)
  if (!is.null(full_swd_data)) rm(full_swd_data)
  if (!load_existing_model && !is.null(tuning_output) && !identical(tuning_output, final_model)) rm(tuning_output)
  gc()
  
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm_biotic_pc1 function


# --- 12. Setup Parallel Backend & Run ---
if (config$use_parallel && config$num_cores > 1) { log4r::info(logger, paste("Setting up parallel backend with", config$num_cores, "cores.")); gc(full=TRUE); future::plan(future::multisession, workers = config$num_cores, gc = TRUE) } else { log4r::info(logger, "Running sequentially."); future::plan(future::sequential) }
progressr::handlers(global = TRUE); progressr::handlers("txtprogressbar")

log4r::info(logger, paste("Starting Biotic+PC1 SDM processing for", nrow(species_df), group_name, "species..."))

results_list <- progressr::with_progress({
  furrr::future_map(1:nrow(species_df), ~{
    process_species_sdm_biotic_pc1( # Call the correct function
      species_row = species_df[.x, ], config = config, env_predictor_paths_list = env_predictor_paths_list,
      all_host_prediction_paths = all_host_prediction_paths, associations_df = associations_df,
      group_name = group_name, predictor_type_suffix = predictor_type_suffix, occurrence_dir = occurrence_dir,
      species_log_dir = config$species_log_dir, tuning_scenario = "current"
    )
  }, .options = furrr_options(seed = TRUE, scheduling = 1.0))
})

log4r::info(logger, "Parallel/sequential processing complete.")

# --- 13. Process Results ---
# (Identical summary logic as 06a/06b)
occurrence_counts <- list(); success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0
log4r::info(logger, paste("--- Processing Results Summary (", group_name, predictor_type_suffix, "Models) ---"))
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
log4r::info(logger, paste("--- Overall Summary (", group_name, predictor_type_suffix, "Models) ---"))
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
  cat("\n--- Final Occurrence Counts (", group_name, predictor_type_suffix, ", after cleaning/thinning) ---\n"); print(occ_count_df_print)
} else { log4r::warn(logger, "No occurrence counts were recorded.") }

future::plan(future::sequential); gc(full=TRUE)
log4r::info(logger, "--- Script 06c (v6 - Final VI Call Fix) finished. ---") # Updated version marker
#-------------------------------------------------------------------------------