# scripts/sdm_runs/experiments/tweaking_base_run.R

cat("--- Running Script tweaking_base_run.R ---\n")

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
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))

# --- 4. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path,
                       log_level = config$log_level,
                       append = config$log_append,
                       log_to_console = config$log_to_console,
                       console_level = config$log_console_level)

log4r::info(logger, "--- Starting Script tweaking_base_run.R ---")

# --- 5. Define Group Specifics & Predictor Type ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors
predictor_type_suffix <- ifelse(use_pca, "_pca", "_vif")

# # --- 5. Define Group Specifics & Predictor Type ---
# group_name <- "anemonefish"
# species_list_file <- config$anemonefish_species_list_file
# occurrence_dir <- config$anemonefish_occurrence_dir

log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---"))

# --- 6. Load Predictor Information ---
predictor_paths_or_list <- NULL
if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (is.null(pca_paths_rds) || !is.character(pca_paths_rds) || nchar(pca_paths_rds) == 0 || !file.exists(pca_paths_rds)) {
    log4r::fatal(logger, paste("PCA raster paths RDS file path missing/invalid in config OR file not found at:", pca_paths_rds %||% "NULL"))
    stop("PCA paths RDS missing or invalid.")
  }
  predictor_paths_or_list <- tryCatch(readRDS(pca_paths_rds), error = function(e) {
    log4r::fatal(logger, paste("Failed load PCA paths RDS:", e$message)); NULL
  })
  if (!is.list(predictor_paths_or_list) || length(predictor_paths_or_list) == 0) {
    log4r::fatal(logger, "PCA raster paths list empty/invalid after loading RDS.")
    stop("PCA paths list invalid.")
  }
  log4r::info(logger, paste("Loaded PCA raster paths for scenarios:", paste(names(predictor_paths_or_list), collapse=", ")))
} else {
  predictor_paths_or_list <- config$final_vars_vif_anemone
  if(is.null(predictor_paths_or_list) || !is.character(predictor_paths_or_list) || length(predictor_paths_or_list) < 2) {
    log4r::fatal(logger, "Final VIF var list `final_vars_vif_anemone` missing, invalid, or too short.")
    stop("VIF vars missing.")
  }
  log4r::info(logger, paste("Using VIF-selected core variables:", paste(predictor_paths_or_list, collapse=", ")))
}

# --- 7. Create *Intermediate* Output Dirs NOW (with checks on paths from config) ---
log4r::debug(logger, "Attempting to create intermediate output directories...")

# Check base paths BEFORE constructing full paths
base_intermediate_model_path <- config$models_dir_intermediate # Use the explicitly included intermediate path
base_intermediate_results_path <- config$results_dir_intermediate # Use the explicitly included intermediate path
base_species_log_path <- config$species_log_dir

if (is.null(base_intermediate_model_path) || !is.character(base_intermediate_model_path) || length(base_intermediate_model_path) != 1 || nchar(base_intermediate_model_path) == 0) stop("FATAL: config$models_dir_intermediate is missing or invalid.")
if (is.null(base_intermediate_results_path) || !is.character(base_intermediate_results_path) || length(base_intermediate_results_path) != 1 || nchar(base_intermediate_results_path) == 0) stop("FATAL: config$results_dir_intermediate is missing or invalid.")
if (is.null(base_species_log_path) || !is.character(base_species_log_path) || length(base_species_log_path) != 1 || nchar(base_species_log_path) == 0) stop("FATAL: config$species_log_dir is missing or invalid.")
if (is.null(group_name) || !is.character(group_name) || nchar(group_name) == 0) stop("FATAL: group_name is invalid.")
if (is.null(predictor_type_suffix) || !is.character(predictor_type_suffix) ) stop("FATAL: predictor_type_suffix is invalid.")

# Construct intermediate paths *after* checks
intermediate_models_dir <- file.path(base_intermediate_model_path, paste0(group_name, predictor_type_suffix))
intermediate_results_dir <- file.path(base_intermediate_results_path, paste0(group_name, predictor_type_suffix))

log4r::debug(logger, paste("Constructed intermediate models dir path:", intermediate_models_dir))
log4r::debug(logger, paste("Constructed intermediate results dir path:", intermediate_results_dir))

# Attempt directory creation with error handling
tryCatch({
  dir.create(intermediate_models_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(intermediate_results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_species_log_path, recursive = TRUE, showWarnings = FALSE) # Ensure base species log dir exists
  log4r::debug(logger, paste("Intermediate output directories created/checked for:", group_name, predictor_type_suffix))
}, error = function(e) {
  log4r::fatal(logger, paste("Failed to create intermediate directories. Path Model:", intermediate_models_dir, "Path Results:", intermediate_results_dir, "Path Species Log:", base_species_log_path, "Error:", e$message))
  stop("Failed to create intermediate directories. Check paths and permissions.")
})
# --- End Section 7 ---


# --- 8. Load Species List ---
tryCatch({
  species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) {
  log4r::fatal(logger, paste("Failed load species list:", e$message))
  stop("Species list loading failed.")
})
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))



# --- 9. Define Function to Process Single Species (v5 - Tweaking Integrated) ---
process_species_sdm <- function(species_row, config, predictor_paths_or_list, group_name, predictor_type_suffix, use_pca, occurrence_dir,
                                tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  species_log_file <- file.path(config$species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), paste0("[", species_name, "]"), paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
  cat("INFO --- Starting processing (Script 06a v5 - Tweaking Integrated) ---")
  
  # --- Define Intermediate Paths ---
  intermediate_models_dir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
  intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix))
  final_model_file <- file.path(intermediate_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  tuning_rds_file <- file.path(intermediate_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, "_object.rds"))
  
  # --- Load GLOBAL Predictor Stack for Tuning Scenario (Needed for SAC & Background Gen) ---
  cat("DEBUG Loading GLOBAL tuning predictors for scenario:", tuning_scenario)
  tuning_predictor_stack_global <- NULL
  if(use_pca){
    tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
    if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) { msg <- paste0("Skipping: GLOBAL PCA stack path for tuning scenario '", tuning_scenario, "' not found."); cat("ERROR", msg); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
    tuning_predictor_stack_global <- tryCatch(terra::rast(tuning_predictor_path), error = function(e) {cat("ERROR Failed to load GLOBAL tuning PCA stack:", e$message); NULL})
  } else {
    # This part is currently inactive based on config$use_pca_predictors = TRUE
    current_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, tuning_scenario, config)
    tuning_predictor_stack_global <- load_selected_env_data(tuning_scenario, current_vif_vars, config)
  }
  if(is.null(tuning_predictor_stack_global)) { msg <- paste0("Skipping: Failed load GLOBAL predictor stack for tuning scenario '", tuning_scenario, "'."); cat("ERROR", msg); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
  cat("DEBUG GLOBAL Tuning predictor stack loaded.")
  raster_crs_terra <- terra::crs(tuning_predictor_stack_global)
  
  # --- Load and Clean Initial Occurrences (NO CELL THINNING HERE) ---
  config_for_occ_load <- config # Copy config
  # Ensure cell thinning is disabled in the loading function call
  config_for_occ_load$thinning_method <- NULL
  cat("DEBUG Loading/cleaning initial occurrences (before SAC thinning).")
  occ_data_result_raw <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger=NULL, species_log_file=species_log_file)
  
  if (is.null(occ_data_result_raw) || is.null(occ_data_result_raw$coords) || occ_data_result_raw$count < config$min_occurrences_sdm) {
    msg <- paste0("Skipping: Insufficient occurrences before SAC thinning (", occ_data_result_raw$count %||% 0, ").");
    cat("WARN", msg);
    return(list(status = "skipped_occurrences_raw", species = species_name, occurrence_count = occ_data_result_raw$count %||% 0, message = msg))
  }
  occs_coords_raw <- occ_data_result_raw$coords
  occurrence_count_raw <- occ_data_result_raw$count
  cat("INFO", paste("Occurrence count after basic cleaning:", occurrence_count_raw))
  
  # --- <<< Spatial Autocorrelation (SAC) Thinning >>> ---
  occs_coords_thinned <- occs_coords_raw # Default to raw if thinning skipped/fails
  occurrence_count_after_thinning <- occurrence_count_raw # Use raw count initially
  thinning_dist_applied <- NA
  
  sac_thin_result <- thin_occurrences_by_sac(
    occs_coords = occs_coords_raw,
    predictor_stack = tuning_predictor_stack_global, # Use the GLOBAL stack for env data extraction during SAC check
    config = config,
    logger = NULL, # Pass NULL logger to helper if needed
    species_log_file = species_log_file
  )
  
  if (!is.null(sac_thin_result)) {
    occs_coords_thinned <- sac_thin_result$coords_thinned
    occurrence_count_after_thinning <- sac_thin_result$n_thinned # Update count
    thinning_dist_applied <- sac_thin_result$thinning_distance_km
    cat("INFO", paste("Occurrence count after SAC thinning:", occurrence_count_after_thinning, "(Distance:", thinning_dist_applied %||% "NA", "km)"))
    # Check if enough points remain *after* SAC thinning
    if(occurrence_count_after_thinning < config$min_occurrences_sdm) {
      msg <- paste0("Skipping: Insufficient occurrences after SAC thinning (", occurrence_count_after_thinning, ").");
      cat("WARN", msg);
      return(list(status = "skipped_occurrences_sac", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
  } else {
    cat("WARN SAC thinning function returned NULL or failed. Proceeding with unthinned data.")
    # occs_coords_thinned and occurrence_count_after_thinning retain original values
  }
  # --- <<< END SPATIAL AUTOCORRELATION THINNING >>> ---
  
  # --- Use thinned coordinates for subsequent steps ---
  occs_coords <- occs_coords_thinned # Use potentially thinned coords from now on
  # Convert to sf object for background generation
  occs_sf_clean <- sf::st_as_sf(as.data.frame(occs_coords), coords=c("longitude","latitude"), crs=config$occurrence_crs)
  # Add AphiaID column needed by generate_sdm_background_obis
  occs_sf_clean$AphiaID <- species_aphia_id
  # Ensure CRS matches raster
  if(sf::st_crs(occs_sf_clean) != sf::st_crs(raster_crs_terra)){
    occs_sf_clean <- sf::st_transform(occs_sf_clean, crs=sf::st_crs(raster_crs_terra))
  }
  
  # --- Check Existing Model (Optional) ---
  final_model <- NULL; tuning_output <- NULL; load_existing_model <- FALSE
  full_swd_data <- NULL # Initialize SWD object
  species_specific_stack <- NULL # Initialize species specific stack
  
  if (!config$force_rerun$run_standard_sdms && file.exists(final_model_file)) {
    cat("INFO Final model file exists. Attempting to load...")
    final_model <- tryCatch(readRDS(final_model_file), error=function(e) { cat("ERROR Failed to load existing model:", e$message); NULL})
    if(!is.null(final_model) && inherits(final_model, "SDMmodel")){
      cat("INFO Successfully loaded existing final model. Skipping tuning/training.")
      load_existing_model <- TRUE
      
      # Regenerate background points and species stack for VI/Eval using OBIS method
      cat("DEBUG Regenerating background points and species stack for loaded model VI/Eval...")
      background_return_eval <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack_global, config, logger=NULL, species_log_file=species_log_file, seed_offset = species_aphia_id)
      
      if(!is.null(background_return_eval) && !is.null(background_return_eval$background_points) && !is.null(background_return_eval$species_specific_stack)) {
        background_points_eval <- background_return_eval$background_points
        species_specific_stack <- background_return_eval$species_specific_stack # Assign to outer scope variable for later use
        full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_points_eval, env = species_specific_stack, verbose = FALSE) }, error = function(e) { cat("WARN Failed prepareSWD for VI/Eval on loaded model:", e$message); NULL })
        rm(background_points_eval); gc() # Clean up local eval points
      } else { cat("WARN Failed background/stack generation for VI/Eval on loaded model.") }
      
      # Load tuning object if it exists, needed for test metrics
      if(file.exists(tuning_rds_file)) {
        tuning_output <- tryCatch(readRDS(tuning_rds_file), error=function(e){cat("WARN  Could not load tuning RDS."); NULL})
      } else { cat("WARN Tuning RDS file not found, cannot log test metrics for loaded model.") }
      if(is.null(tuning_output)) { tuning_output <- NULL } # Ensure it's NULL if loading failed
      
    } else {
      msg <- paste0("Existing final model file invalid/failed load. Will re-run."); cat("WARN", msg);
      final_model <- NULL # Ensure final_model is NULL if loading failed
      load_existing_model <- FALSE
    }
  } # End check for existing model
  
  # --- Tuning & Training (if not loading existing) ---
  if (!load_existing_model) {
    cat("INFO Final model not found/invalid or rerun forced. Proceeding with tuning/training.")
    
    cat("DEBUG Generating background points and species stack for training...")
    background_return <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack_global, config, logger=NULL, species_log_file=species_log_file, seed_offset = species_aphia_id)
    if (is.null(background_return) || is.null(background_return$background_points) || is.null(background_return$species_specific_stack)) {
      msg <- paste0("Skipping: Failed background point/stack generation for training."); cat("ERROR", msg);
      return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    background_points <- background_return$background_points
    species_specific_stack <- background_return$species_specific_stack # Use this masked stack
    cat("DEBUG Background points generated and stack masked for training.")
    
    # Prepare SWD data using the species_specific_stack
    full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_points, env = species_specific_stack, verbose = FALSE) }, error = function(e) { cat("ERROR", paste("Failed prepareSWD for tuning/training:", e$message)); NULL })
    if (is.null(full_swd_data)) { rm(background_points, species_specific_stack); gc(); return(list(status = "error_swd_preparation", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("SWD prep failed."))) }
    cat("DEBUG SWD object prepared for tuning/training.")
    
    # Create Spatial Folds (using species_specific_stack)
    cat("DEBUG Creating spatial folds...")
    spatial_folds <- create_spatial_cv_folds_simplified(full_swd_data, species_specific_stack, config, logger=NULL, species_log_file)
    if(is.null(spatial_folds)) {
      msg <- "Failed to create spatial folds."; cat("ERROR", msg);
      rm(background_points, species_specific_stack, full_swd_data); gc();
      return(list(status = "error_cv_folds", species=species_name, occurrence_count = occurrence_count_after_thinning, message=msg))
    }
    
    # Tune Hyperparameters (using species_specific_stack via SWD data & folds)
    cat("INFO Starting hyperparameter tuning.")
    tuning_output <- run_sdm_tuning_scv(occs_coords, species_specific_stack, background_points, config, logger=NULL, species_name, species_log_file=species_log_file)
    
    if (is.null(tuning_output) || is.null(attr(tuning_output, "best_hypers"))) {
      msg <- paste0("Skipping: Tuning failed."); cat("ERROR", msg);
      rm(background_points, species_specific_stack, full_swd_data, spatial_folds); gc();
      return(list(status = "error_tuning", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    best_hypers <- attr(tuning_output, "best_hypers")
    
    # # Save Tuning Results
    # if(!save_tuning_results(tuning_output, species_name_sanitized, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)) {
    #   rm(background_points, species_specific_stack, full_swd_data, spatial_folds, tuning_output); gc();
    #   return(list(status = "error_saving_tuning_results", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("Failed save tuning results.")))
    # }
    
    # Train Final Model (using species_specific_stack)
    cat("INFO Starting final model training.")
    final_model <- train_final_sdm(occs_coords, species_specific_stack, background_points, best_hypers, config, logger=NULL, species_name, species_log_file=species_log_file)
    
    if(is.null(final_model)){
      msg <- paste0("Skipping: Training failed."); cat("ERROR", msg);
      rm(background_points, species_specific_stack, full_swd_data, spatial_folds, tuning_output); gc();
      return(list(status = "error_training", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    cat("INFO Final model training complete.")
    
    # Save Final Model (to intermediate location)
    if(!save_final_model(final_model, species_name_sanitized, predictor_type_suffix, group_name, config, logger=NULL, species_log_file=species_log_file)) {
      rm(background_points, species_specific_stack, full_swd_data, spatial_folds, tuning_output, final_model); gc();
      return(list(status = "error_saving_model", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("Failed save final model.")))
    }
    
    rm(background_points, spatial_folds); gc() # Keep species_specific_stack & full_swd_data for logging/VI
    
  } # End if !load_existing_model
  
  # --- Log Metrics (Ensure this happens AFTER model training or loading) ---
  if (!is.null(final_model) && !is.null(species_specific_stack)) { # Need stack for AICc
    log_final_model_metrics(
      final_model = final_model,
      full_swd_data = full_swd_data, # Pass SWD data (might be NULL if loaded model and SWD regen failed)
      tuning_predictor_stack = species_specific_stack, # Pass stack used for training/AICc
      tuning_output = tuning_output, # Pass tuning results (might be NULL if loaded model and RDS missing)
      species_name_sanitized = species_name_sanitized,
      group_name = group_name,
      predictor_type_suffix = predictor_type_suffix,
      config = config,
      logger = NULL, # Pass NULL logger to helper if needed
      species_log_file = species_log_file
    )
  } else {
    cat("WARN Final model or species-specific stack object is NULL, cannot log metrics.")
  }
  
  # --- Variable Importance ---
  if (!is.null(final_model)) {
    if(!is.null(full_swd_data)){ # Check if SWD data is available
      cat("INFO Calculating and saving variable importance...")
      vi_success <- calculate_and_save_vi(
        final_model = final_model,
        training_swd = full_swd_data,
        species_name_sanitized = species_name_sanitized,
        group_name = group_name,
        predictor_type_suffix = predictor_type_suffix,
        config = config,
        logger = NULL,
        species_log_file = species_log_file)
      if(!vi_success) cat("WARN Variable importance calculation/saving failed.")
    } else {
      cat("WARN Skipping variable importance: training SWD data is missing.")
    }
  } else {
    cat("WARN Skipping variable importance: final model object is unavailable.")
  }
  
  # --- Prediction Loop ---
  # Requires the GLOBAL predictor stacks, not species_specific_stack
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- if(use_pca) names(predictor_paths_or_list) else config$env_scenarios
  cat("INFO", paste("Starting predictions for", length(scenarios_to_predict), "scenarios."))
  
  if (!is.null(final_model)) {
    for (pred_scenario in scenarios_to_predict) {
      cat("DEBUG", paste("  Predicting scenario:", pred_scenario))
      target_pred_file_path <- construct_prediction_filename(species_name_sanitized, pred_scenario, predictor_type_suffix, config)
      if (is.null(target_pred_file_path)) { cat("WARN   Could not construct prediction filename. Skipping."); prediction_errors <- prediction_errors + 1; next }
      
      # Check force_rerun flag
      rerun_flag <- config$force_rerun$run_standard_sdms
      if (!rerun_flag && file.exists(target_pred_file_path)) {
        cat("DEBUG   Prediction exists in target location. Skipping."); next
      }
      
      # --- Load GLOBAL predictor stack for this prediction scenario ---
      pred_predictor_stack_global <- NULL
      if(use_pca) {
        pred_predictor_path <- predictor_paths_or_list[[pred_scenario]]
        if (is.null(pred_predictor_path) || !file.exists(pred_predictor_path)) { cat("WARN", paste("  GLOBAL PCA stack path missing:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
        pred_predictor_stack_global <- tryCatch(terra::rast(pred_predictor_path), error = function(e) {cat("WARN", paste("  Error loading GLOBAL PCA stack:", e$message)); NULL})
      } else {
        # VIF logic (currently inactive)
        scenario_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, pred_scenario, config)
        if(length(scenario_vif_vars) < 1) { cat("WARN", paste("  No VIF vars for scenario:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
        pred_predictor_stack_global <- load_selected_env_data(pred_scenario, scenario_vif_vars, config) # Assumes this loads global extent
      }
      if(is.null(pred_predictor_stack_global)) { cat("WARN", paste("  Failed load GLOBAL predictor stack:", pred_scenario)); prediction_errors <- prediction_errors + 1; next }
      # --- END Load GLOBAL predictor stack ---
      
      # Predict using the final model and the GLOBAL stack
      prediction_output <- predict_sdm_suitability(final_model, pred_predictor_stack_global, config, logger=NULL, species_log_file=species_log_file)
      
      if (inherits(prediction_output, "SpatRaster")) {
        save_success <- save_sdm_prediction(prediction_output, species_name_sanitized, pred_scenario, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)
        if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
        rm(prediction_output); gc()
      } else { cat("WARN", paste("  Prediction failed:", prediction_output)); prediction_errors <- prediction_errors + 1 }
      rm(pred_predictor_stack_global); gc() # Clean up the global stack for this scenario
    } # End prediction scenario loop
  } else {
    prediction_errors <- length(scenarios_to_predict); msg <- paste0("Skipping all predictions: Final model was not available."); cat("ERROR", msg)
    status <- if(load_existing_model) "error_loading_model" else "error_training"; return(list(status = status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
  }
  
  # --- Prepare return status ---
  final_status <- "success"
  status_message <- paste0("Finished Std Env SDM. Occs:", occurrence_count_after_thinning, ". Preds attempted:", length(scenarios_to_predict), ". Made:", predictions_made, ". Errors/Skipped:", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) { final_status <- "error_prediction_all"; cat("ERROR", status_message) }
  else if (prediction_errors > 0) { final_status <- "success_with_pred_errors"; cat("WARN", status_message) }
  else { cat("INFO", status_message) }
  
  # --- Clean up ---
  if (!is.null(species_specific_stack)) rm(species_specific_stack)
  if (!is.null(tuning_predictor_stack_global)) rm(tuning_predictor_stack_global)
  if (!is.null(occs_coords)) rm(occs_coords)
  if (!is.null(full_swd_data)) rm(full_swd_data)
  if (!load_existing_model && !is.null(tuning_output) && inherits(tuning_output, "SDMtune")) rm(tuning_output)
  if (exists("final_model", inherits=FALSE)) rm(final_model)
  gc()
  
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm function



# --- 10. Setup Run ---
log4r::info(logger, paste("Starting SDM processing for", nrow(species_df), group_name, "species..."))
res <- process_species_sdm(
  species_row = species_df[1, ],
  config = config,
  predictor_paths_or_list = predictor_paths_or_list,
  group_name = group_name,
  predictor_type_suffix = predictor_type_suffix,
  use_pca = use_pca,
  occurrence_dir = occurrence_dir,
  tuning_scenario = "current"
)


# --- 11. Process Results ---
# (Result processing and summary logic remains the same as v4)
occurrence_counts <- list()
success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0
log4r::info(logger, paste("--- Processing Results Summary (", group_name, predictor_type_suffix, "Models) ---"))
if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "Received NULL result."); next }
log_level_func <- switch(res$status, success=log4r::info, success_with_pred_errors=log4r::warn, skipped_occurrences=log4r::warn, log4r::error)
log_level_func(logger, res$message)
occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
if (res$status == "success") success_count <- success_count + 1

log4r::info(logger, paste("--- Overall Summary (", group_name, predictor_type_suffix, "Models) ---"))
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

gc(full=TRUE)
log4r::info(logger, "--- Script tweaking_base_run.R finished. ---")
#-------------------------------------------------------------------------------











# 
# # --- 10. Setup Parallel Backend & Run ---
# # (Parallel setup and furrr::future_map call remain the same as v4)
# if (config$use_parallel && config$num_cores > 1) {
#   log4r::info(logger, paste("Setting up parallel backend with", config$num_cores, "cores (multisession)."))
#   gc(full=TRUE); future::plan(future::multisession, workers = config$num_cores, gc = TRUE)
# } else {
#   log4r::info(logger, "Running sequentially."); future::plan(future::sequential)
# }
# progressr::handlers(global = TRUE); progressr::handlers("txtprogressbar")
# 
# log4r::info(logger, paste("Starting SDM processing for", nrow(species_df), group_name, "species..."))
# 
# results_list <- progressr::with_progress({
#   furrr::future_map(1:nrow(species_df), ~{
#     process_species_sdm(
#       species_row = species_df[.x, ],
#       config = config,
#       predictor_paths_or_list = predictor_paths_or_list,
#       group_name = group_name,
#       predictor_type_suffix = predictor_type_suffix,
#       use_pca = use_pca,
#       occurrence_dir = occurrence_dir,
#       tuning_scenario = "current"
#     )
#   }, .options = furrr_options(seed = TRUE, scheduling = 1.0))
# })
# 
# log4r::info(logger, "Parallel/sequential processing complete.")
# 
# 
# # --- 11. Process Results ---
# # (Result processing and summary logic remains the same as v4)
# occurrence_counts <- list()
# success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0
# log4r::info(logger, paste("--- Processing Results Summary (", group_name, predictor_type_suffix, "Models) ---"))
# for (res in results_list) {
#   if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "Received NULL result."); next }
#   log_level_func <- switch(res$status, success=log4r::info, success_with_pred_errors=log4r::warn, skipped_occurrences=log4r::warn, log4r::error)
#   log_level_func(logger, res$message)
#   occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
#   if (res$status == "success") success_count <- success_count + 1
#   else if (res$status == "success_with_pred_errors") partial_success_count <- partial_success_count + 1
#   else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
#   else error_count <- error_count + 1
# }
# log4r::info(logger, paste("--- Overall Summary (", group_name, predictor_type_suffix, "Models) ---"))
# log4r::info(logger, paste("Total Species Targets:", length(results_list)))
# log4r::info(logger, paste("Fully Successful Runs:", success_count))
# log4r::info(logger, paste("Partially Successful Runs:", partial_success_count))
# log4r::info(logger, paste("Skipped:", skipped_count))
# log4r::info(logger, paste("Errors:", error_count))
# if (length(occurrence_counts) > 0) {
#   occ_count_df <- data.frame( Species = names(occurrence_counts), OccurrenceCountAfterThinning = unlist(occurrence_counts)) %>% dplyr::arrange(Species)
#   occ_count_df_print <- occ_count_df; occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
#   occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, ".csv"))
#   tryCatch({ readr::write_csv(occ_count_df, occ_count_file); log4r::info(logger, paste("Occurrence counts saved:", occ_count_file))}, error = function(e) { log4r::error(logger, paste("Failed save counts:", e$message)) })
#   cat("\n--- Final Occurrence Counts (", group_name, predictor_type_suffix, ", after cleaning/thinning) ---\n"); print(occ_count_df_print)
# } else { log4r::warn(logger, "No occurrence counts were recorded.") }
# 
# future::plan(future::sequential); gc(full=TRUE)
# log4r::info(logger, "--- Script 06a (v7 - Corrected VI Call) finished. ---")
# #-------------------------------------------------------------------------------