# scripts/sdm_runs/experiments/tweaking_base_run.R

cat("--- Running Script tweaking_base_run.R ---\n")


path_to_maxent_jar_local <- "/home/bi-server-kyoto/R/x86_64-pc-linux-gnu-library/4.4/dismo/java/maxent.jar" 
# If you installed dismo and it includes maxent.jar, it might be:
# path_to_maxent_jar_local <- file.path(system.file(package="dismo"), "java", "maxent.jar")
# Check if it exists:
if (!file.exists(path_to_maxent_jar_local)) {
  cat("WARNING: maxent.jar not found at specified path for MAXENT.Phillips:", path_to_maxent_jar_local, "\n")
  # You might want to stop or only run MAXENT.Phillips.2 if jar is missing
}


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


# --- 9. Define Function to Process Single Species (v5 - Tweaking Integrated, with BIOMOD2 step) ---
process_species_sdm <- function(species_row, config, predictor_paths_or_list, group_name, predictor_type_suffix, use_pca, occurrence_dir,
                                tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  # Use cat for logging within this function for the tweaking script
  log_prefix <- paste(Sys.time(), paste0("[", species_name, "]"))
  cat(log_prefix, "INFO --- Starting processing (tweaking_base_run.R with BIOMOD2 step) ---\n")
  
  # --- Define Intermediate Paths (SDMtune) ---
  # (Keep this section as is)
  intermediate_models_dir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
  intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix))
  final_model_file <- file.path(intermediate_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  tuning_rds_file <- file.path(intermediate_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, "_object.rds"))
  
  # --- Load GLOBAL Predictor Stack for Tuning Scenario ---
  # (Keep this section as is)
  cat(log_prefix, "DEBUG Loading GLOBAL tuning predictors for scenario:", tuning_scenario, "\n")
  tuning_predictor_stack_global <- NULL
  if(use_pca){
    tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
    if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) { msg <- paste0("Skipping: GLOBAL PCA stack path for tuning scenario '", tuning_scenario, "' not found."); cat(log_prefix, "ERROR", msg, "\n"); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
    tuning_predictor_stack_global <- tryCatch(terra::rast(tuning_predictor_path), error = function(e) {cat(log_prefix, "ERROR Failed to load GLOBAL tuning PCA stack:", e$message, "\n"); NULL})
  } else {
    current_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, tuning_scenario, config)
    tuning_predictor_stack_global <- load_selected_env_data(tuning_scenario, current_vif_vars, config)
  }
  if(is.null(tuning_predictor_stack_global)) { msg <- paste0("Skipping: Failed load GLOBAL predictor stack for tuning scenario '", tuning_scenario, "'."); cat(log_prefix, "ERROR", msg, "\n"); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
  cat(log_prefix, "DEBUG GLOBAL Tuning predictor stack loaded.\n")
  raster_crs_terra <- terra::crs(tuning_predictor_stack_global)
  
  # --- Load and Clean Initial Occurrences (NO CELL THINNING HERE) ---
  # (Keep this section as is)
  config_for_occ_load <- config
  config_for_occ_load$thinning_method <- NULL
  cat(log_prefix, "DEBUG Loading/cleaning initial occurrences (before SAC thinning).\n")
  occ_data_result_raw <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger=NULL, species_log_file=NULL) # Pass NULL for species_log_file to use cat
  
  if (is.null(occ_data_result_raw) || is.null(occ_data_result_raw$coords) || occ_data_result_raw$count < config$min_occurrences_sdm) {
    msg <- paste0("Skipping: Insufficient occurrences before SAC thinning (", occ_data_result_raw$count %||% 0, ").");
    cat(log_prefix, "WARN", msg, "\n");
    return(list(status = "skipped_occurrences_raw", species = species_name, occurrence_count = occ_data_result_raw$count %||% 0, message = msg))
  }
  occs_coords_raw <- occ_data_result_raw$coords
  occurrence_count_raw <- occ_data_result_raw$count
  cat(log_prefix, "INFO Occurrence count after basic cleaning:", occurrence_count_raw, "\n")
  
  # --- Spatial Autocorrelation (SAC) Thinning ---
  # (Keep this section as is)
  occs_coords_thinned <- occs_coords_raw
  occurrence_count_after_thinning <- occurrence_count_raw
  thinning_dist_applied <- NA
  sac_thin_result <- thin_occurrences_by_sac(occs_coords_raw, tuning_predictor_stack_global, config, logger=NULL, species_log_file=NULL) # Pass NULL for species_log_file
  
  if (!is.null(sac_thin_result)) {
    occs_coords_thinned <- sac_thin_result$coords_thinned
    occurrence_count_after_thinning <- sac_thin_result$n_thinned
    thinning_dist_applied <- sac_thin_result$thinning_distance_km
    cat(log_prefix, "INFO Occurrence count after SAC thinning:", occurrence_count_after_thinning, "(Distance:", thinning_dist_applied %||% "NA", "km)\n")
    if(occurrence_count_after_thinning < config$min_occurrences_sdm) {
      msg <- paste0("Skipping: Insufficient occurrences after SAC thinning (", occurrence_count_after_thinning, ").");
      cat(log_prefix, "WARN", msg, "\n");
      return(list(status = "skipped_occurrences_sac", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
  } else {
    cat(log_prefix, "WARN SAC thinning function returned NULL or failed. Proceeding with unthinned data.\n")
  }
  
  # --- Use thinned coordinates for subsequent steps ---
  # (Keep this section as is)
  occs_coords <- occs_coords_thinned
  occs_sf_clean <- sf::st_as_sf(as.data.frame(occs_coords), coords=c("longitude","latitude"), crs=config$occurrence_crs)
  occs_sf_clean$AphiaID <- species_aphia_id
  if(sf::st_crs(occs_sf_clean) != sf::st_crs(raster_crs_terra)){
    occs_sf_clean <- sf::st_transform(occs_sf_clean, crs=sf::st_crs(raster_crs_terra))
  }
  
  # ===========================================================================
  # --- SDMtune Workflow (Existing Logic) ---
  # ===========================================================================
  cat(log_prefix, "INFO --- Starting SDMtune Workflow ---\n")
  sdmtune_final_model <- NULL; sdmtune_tuning_output <- NULL; sdmtune_load_existing_model <- FALSE
  sdmtune_full_swd_data <- NULL; sdmtune_species_specific_stack <- NULL; sdmtune_background_points <- NULL # For BIOMOD2 too
  sdmtune_best_hypers <- NULL # To store best hyperparameters from SDMtune for BIOMOD2
  
  if (!config$force_rerun$run_standard_sdms && file.exists(final_model_file)) {
    cat(log_prefix, "INFO [SDMtune] Final model file exists. Attempting to load...\n")
    sdmtune_final_model <- tryCatch(readRDS(final_model_file), error=function(e) { cat(log_prefix, "ERROR [SDMtune] Failed to load existing model:", e$message, "\n"); NULL})
    if(!is.null(sdmtune_final_model) && inherits(sdmtune_final_model, "SDMmodel")){
      cat(log_prefix, "INFO [SDMtune] Successfully loaded existing final model.\n")
      sdmtune_load_existing_model <- TRUE
      
      # Regenerate background points and species stack - CRUCIAL for BIOMOD2 later
      cat(log_prefix, "DEBUG [SDMtune] Regenerating background points and species stack for loaded model...\n")
      background_return_eval <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack_global, config, logger=NULL, species_log_file=NULL, seed_offset = species_aphia_id)
      if(!is.null(background_return_eval) && !is.null(background_return_eval$background_points) && !is.null(background_return_eval$species_specific_stack)) {
        sdmtune_background_points <- background_return_eval$background_points # Store for BIOMOD2
        sdmtune_species_specific_stack <- background_return_eval$species_specific_stack # Store for BIOMOD2 and VI/Metrics
        sdmtune_full_swd_data <- tryCatch({ SDMtune::prepareSWD(species=species_name, p=occs_coords, a=sdmtune_background_points, env=sdmtune_species_specific_stack, verbose=FALSE) }, error = function(e) { cat(log_prefix, "WARN [SDMtune] Failed prepareSWD for VI/Eval:", e$message, "\n"); NULL })
      } else { cat(log_prefix, "WARN [SDMtune] Failed background/stack generation for loaded model.\n") }
      
      if(file.exists(tuning_rds_file)) {
        sdmtune_tuning_output <- tryCatch(readRDS(tuning_rds_file), error=function(e){cat(log_prefix, "WARN [SDMtune] Could not load tuning RDS.\n"); NULL})
        if(!is.null(sdmtune_tuning_output) && !is.null(attr(sdmtune_tuning_output, "best_hypers"))) {
          sdmtune_best_hypers <- attr(sdmtune_tuning_output, "best_hypers") # Get hypers for BIOMOD2
        }
      } else { cat(log_prefix, "WARN [SDMtune] Tuning RDS file not found.\n") }
    } else { cat(log_prefix, "WARN [SDMtune] Existing final model file invalid. Will re-run.\n"); sdmtune_final_model <- NULL; sdmtune_load_existing_model <- FALSE }
  }
  
  if (!sdmtune_load_existing_model) {
    cat(log_prefix, "INFO [SDMtune] Proceeding with tuning/training.\n")
    background_return <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack_global, config, logger=NULL, species_log_file=NULL, seed_offset = species_aphia_id)
    if (is.null(background_return) || is.null(background_return$background_points) || is.null(background_return$species_specific_stack)) { msg <- "[SDMtune] Failed background/stack generation."; cat(log_prefix, "ERROR", msg, "\n"); return(list(status="error_background_sdmtune", species=species_name, occurrence_count=occurrence_count_after_thinning, message=msg)) }
    sdmtune_background_points <- background_return$background_points # Store for BIOMOD2
    sdmtune_species_specific_stack <- background_return$species_specific_stack # Store for BIOMOD2 and SDMtune
    
    sdmtune_full_swd_data <- tryCatch({ SDMtune::prepareSWD(species=species_name, p=occs_coords, a=sdmtune_background_points, env=sdmtune_species_specific_stack, verbose=FALSE) }, error = function(e) { cat(log_prefix, "ERROR [SDMtune] Failed prepareSWD:", e$message, "\n"); NULL })
    if (is.null(sdmtune_full_swd_data)) { return(list(status="error_swd_sdmtune", species=species_name, occurrence_count=occurrence_count_after_thinning, message="[SDMtune] SWD prep failed.")) }
    
    spatial_folds <- create_spatial_cv_folds_simplified(sdmtune_full_swd_data, sdmtune_species_specific_stack, config, logger=NULL, species_log_file=NULL)
    if(is.null(spatial_folds)) { msg <- "[SDMtune] Failed spatial folds."; cat(log_prefix, "ERROR", msg, "\n"); return(list(status="error_cv_folds_sdmtune", species=species_name, occurrence_count=occurrence_count_after_thinning, message=msg)) }
    
    sdmtune_tuning_output <- run_sdm_tuning_scv(occs_coords, sdmtune_species_specific_stack, sdmtune_background_points, config, logger=NULL, species_name, species_log_file=NULL)
    if (is.null(sdmtune_tuning_output) || is.null(attr(sdmtune_tuning_output, "best_hypers"))) { msg <- "[SDMtune] Tuning failed."; cat(log_prefix, "ERROR", msg, "\n"); return(list(status="error_tuning_sdmtune", species=species_name, occurrence_count=occurrence_count_after_thinning, message=msg)) }
    sdmtune_best_hypers <- attr(sdmtune_tuning_output, "best_hypers") # Get hypers for BIOMOD2
    save_tuning_results(sdmtune_tuning_output, species_name_sanitized, predictor_type_suffix, config, logger=NULL, species_log_file=NULL)
    
    sdmtune_final_model <- train_final_sdm(occs_coords, sdmtune_species_specific_stack, sdmtune_background_points, sdmtune_best_hypers, config, logger=NULL, species_name, species_log_file=NULL)
    if(is.null(sdmtune_final_model)){ msg <- "[SDMtune] Training failed."; cat(log_prefix, "ERROR", msg, "\n"); return(list(status="error_training_sdmtune", species=species_name, occurrence_count=occurrence_count_after_thinning, message=msg)) }
    save_final_model(sdmtune_final_model, species_name_sanitized, predictor_type_suffix, group_name, config, logger=NULL, species_log_file=NULL)
  }
  
  if (!is.null(sdmtune_final_model) && !is.null(sdmtune_species_specific_stack)) {
    log_final_model_metrics(sdmtune_final_model, sdmtune_full_swd_data, sdmtune_species_specific_stack, sdmtune_tuning_output, species_name_sanitized, group_name, predictor_type_suffix, config, logger=NULL, species_log_file=NULL)
    if(!is.null(sdmtune_full_swd_data)){
      calculate_and_save_vi(sdmtune_final_model, sdmtune_full_swd_data, species_name_sanitized, group_name, predictor_type_suffix, config, logger=NULL, species_log_file=NULL)
    } else { cat(log_prefix, "WARN [SDMtune] Skipping VI: SWD data unavailable.\n")}
  } else { cat(log_prefix, "WARN [SDMtune] Model or stack unavailable, skipping metrics/VI.\n") }
  cat(log_prefix, "INFO --- Finished SDMtune Workflow ---\n")
  
  # ===========================================================================
  # --- BIOMOD2 Workflow (First Step: MAXNET with SDMtune Hypers) ---
  # ===========================================================================
  cat(log_prefix, "INFO --- Starting BIOMOD2 Workflow (Default MAXNET with blockCV) ---\n")
  
  if (is.null(occs_coords) || is.null(sdmtune_background_points) || is.null(sdmtune_species_specific_stack)) {
    cat(log_prefix, "ERROR [BIOMOD2] Prerequisite data from SDMtune for BIOMOD2 is missing. Skipping BIOMOD2.\n")
  } else {
    myBiomodData <- format_data_for_biomod2(
      species_name = species_name_sanitized,
      pres_coords = occs_coords,
      bg_coords = sdmtune_background_points,
      env_stack = sdmtune_species_specific_stack, 
      species_log_file = NULL 
    )
    
    if (!is.null(myBiomodData)) {
      cat(log_prefix, "INFO [BIOMOD2] Data formatted successfully.\n")
      
      # Create blockCV table for BIOMOD2
      # Note: sdmtune_species_specific_stack is the one masked to species extent, good for blockCV context
      myBiomodCVTable <- create_biomod2_block_cv_table(
        biomod_formatted_data = myBiomodData,
        predictor_stack = sdmtune_species_specific_stack, # Important: use the stack BIOMOD_FormatingData saw
        config = config, # Pass full config for blockCV settings
        species_log_file = NULL
      )
      
      if (!is.null(myBiomodCVTable)) {
        cat(log_prefix, "INFO [BIOMOD2] blockCV table created successfully.\n")
        
        # Add group_name to config temporarily if helper needs it and it's not passed
        config$group_name_for_biomod_output <- group_name 
        
        myBiomodModelOut <- run_biomod2_default_maxnet_with_blockcv(
          biomod_formatted_data = myBiomodData,
          biomod_cv_table = myBiomodCVTable,
          species_name_sanitized = species_name_sanitized,
          predictor_type_suffix = predictor_type_suffix, 
          config = config,
          species_log_file = NULL
        )
        
        if (!is.null(myBiomodModelOut) && length(myBiomodModelOut@models.computed) > 0) {
          cat(log_prefix, "INFO [BIOMOD2] MAXNET model run successfully with blockCV.\n")
          
          algo_to_project <- 'MAXNET' # We only ran MAXNET
          cat(log_prefix, "INFO [BIOMOD2] Projecting", algo_to_project, "model (current scenario)...\n")
          biomod2_current_projection <- project_biomod2_models_current(
            biomod_model_out = myBiomodModelOut,
            env_stack_current = env_predictor_stack_current, 
            species_name_for_saving = species_name_sanitized_for_files, # The underscore version
            predictor_type_suffix = predictor_type_suffix,   
            selected_algo = algo_short_name, # e.g., 'MAXNET'
            config = cfg, 
            species_log_file = NULL
          )
          if(!is.null(biomod2_current_projection)){
            cat(log_prefix, "INFO [BIOMOD2] Current projection for", algo_to_project, "created.\n")
            save_biomod2_projection(biomod2_current_projection, species_name_sanitized, "current", 
                                    paste0(predictor_type_suffix, "_biomod2_", algo_to_project, "_default"), # Add default suffix
                                    config, logger=NULL, species_log_file=NULL)
          } else { cat(log_prefix, "ERROR [BIOMOD2]", algo_to_project, "projection failed.\n") }
        } else {
          cat(log_prefix, "ERROR [BIOMOD2] MAXNET model run failed or no models computed.\n")
        }
      } else {
        cat(log_prefix, "ERROR [BIOMOD2] Failed to create blockCV table.\n")
      }
    } else {
      cat(log_prefix, "ERROR [BIOMOD2] Data formatting for BIOMOD2 failed.\n")
    }
  }
  cat(log_prefix, "INFO --- Finished BIOMOD2 Workflow (Default MAXNET with blockCV) ---\n")
  
  
  
  
    
  # --- SDMtune Prediction Loop (Keep as is for now) ---
  # (Your existing SDMtune prediction loop)
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- if(use_pca) names(predictor_paths_or_list) else config$env_scenarios
  cat(log_prefix, "INFO [SDMtune] Starting predictions for", length(scenarios_to_predict), "scenarios.\n")
  # ... (rest of your SDMtune prediction loop) ...
  if (!is.null(sdmtune_final_model)) {
    for (pred_scenario in scenarios_to_predict) {
      cat(log_prefix, "DEBUG [SDMtune] Predicting scenario:", pred_scenario, "\n")
      target_pred_file_path <- construct_prediction_filename(species_name_sanitized, pred_scenario, predictor_type_suffix, config)
      if (is.null(target_pred_file_path)) { cat(log_prefix, "WARN [SDMtune] Could not construct prediction filename. Skipping.\n"); prediction_errors <- prediction_errors + 1; next }
      rerun_flag <- config$force_rerun$run_standard_sdms
      if (!rerun_flag && file.exists(target_pred_file_path)) { cat(log_prefix, "DEBUG [SDMtune] Prediction exists. Skipping.\n"); next }
      
      pred_predictor_stack_global <- NULL
      if(use_pca) { pred_predictor_path <- predictor_paths_or_list[[pred_scenario]]; if (is.null(pred_predictor_path) || !file.exists(pred_predictor_path)) { cat(log_prefix, "WARN [SDMtune] PCA stack missing for scenario:", pred_scenario, "\n"); prediction_errors <- prediction_errors + 1; next }; pred_predictor_stack_global <- tryCatch(terra::rast(pred_predictor_path), error=function(e)NULL) }
      else { scenario_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, pred_scenario, config); if(length(scenario_vif_vars)<1) { cat(log_prefix, "WARN [SDMtune] No VIF vars for scenario:", pred_scenario, "\n"); prediction_errors <- prediction_errors + 1; next }; pred_predictor_stack_global <- load_selected_env_data(pred_scenario, scenario_vif_vars, config)}
      if(is.null(pred_predictor_stack_global)) { cat(log_prefix, "WARN [SDMtune] Failed load GLOBAL stack for scenario:", pred_scenario, "\n"); prediction_errors <- prediction_errors + 1; next }
      
      prediction_output <- predict_sdm_suitability(sdmtune_final_model, pred_predictor_stack_global, config, logger=NULL, species_log_file=NULL)
      if (inherits(prediction_output, "SpatRaster")) {
        if(save_sdm_prediction(prediction_output, species_name_sanitized, pred_scenario, predictor_type_suffix, config, logger=NULL, species_log_file=NULL)) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
      } else { cat(log_prefix, "WARN [SDMtune] Prediction failed:", prediction_output, "\n"); prediction_errors <- prediction_errors + 1 }
      rm(pred_predictor_stack_global, prediction_output); gc()
    }
  } else { prediction_errors <- length(scenarios_to_predict); msg <- "[SDMtune] Skipping predictions: Final model unavailable."; cat(log_prefix, "ERROR", msg, "\n"); return(list(status = "error_training_sdmtune", species=species_name, occurrence_count=occurrence_count_after_thinning, message=msg)) }
  
  # --- Prepare return status ---
  # (Keep this section as is, it refers to SDMtune predictions)
  final_status <- "success"
  status_message <- paste0("Finished Std Env SDM (SDMtune). Occs:", occurrence_count_after_thinning, ". Preds attempted:", length(scenarios_to_predict), ". Made:", predictions_made, ". Errors/Skipped:", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) { final_status <- "error_prediction_all"; cat(log_prefix, "ERROR", status_message, "\n") }
  else if (prediction_errors > 0) { final_status <- "success_with_pred_errors"; cat(log_prefix, "WARN", status_message, "\n") }
  else { cat(log_prefix, "INFO", status_message, "\n") }
  
  # --- Clean up (SDMtune objects) ---
  # (Keep this section as is)
  if (!is.null(sdmtune_species_specific_stack)) rm(sdmtune_species_specific_stack)
  if (!is.null(tuning_predictor_stack_global)) rm(tuning_predictor_stack_global)
  if (!is.null(occs_coords)) rm(occs_coords)
  if (!is.null(sdmtune_full_swd_data)) rm(sdmtune_full_swd_data)
  if (!sdmtune_load_existing_model && !is.null(sdmtune_tuning_output) && inherits(sdmtune_tuning_output, "SDMtune")) rm(sdmtune_tuning_output)
  if (exists("sdmtune_final_model", inherits=FALSE) && !is.null(sdmtune_final_model)) rm(sdmtune_final_model)
  # Add BIOMOD2 objects to cleanup if they were created
  if (exists("myBiomodData", inherits=FALSE) && !is.null(myBiomodData)) rm(myBiomodData)
  # myBiomodModelOut and biomod2_current_projection are mostly paths/metadata, actual models are on disk.
  gc()
  
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm function


# --- 10. Setup Run (Tweaking Script: Single Species, No Parallel) ---
# (Keep this section as is - it will run the modified process_species_sdm for one species)
config$use_parallel <- FALSE
config$num_cores <- 1
future::plan(future::sequential)
cat("--- LOGGING: Tweaking run: Starting SDM processing for 1 species:", species_df$scientificName[1], "---\n")
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

# --- 11. Process Results (Tweaking Script: Single Result) ---
# (Keep this section as is)
occurrence_counts <- list(); success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0
cat("--- LOGGING: Tweaking run: Processing Results Summary (", group_name, predictor_type_suffix, "Models) ---\n")
if (is.null(res)) { 
  error_count <- error_count + 1; 
  cat("--- LOGGING: Tweaking run: ERROR - Received NULL result for the single species run.---\n")
  # If `next` was intended for a loop, it's not applicable here.
  # For a single run, we just note the error.
} else {
  # Determine log level based on status (simplified for cat)
  log_level_prefix <- switch(res$status, 
                             success="INFO", 
                             success_with_pred_errors="WARN", 
                             skipped_occurrences_raw="WARN", 
                             skipped_occurrences_sac="WARN", 
                             "ERROR") # Default to ERROR
  cat("--- LOGGING: Tweaking run:", log_level_prefix, "-", res$message, "---\n")
  
  occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
  if (res$status == "success") success_count <- success_count + 1
  else if (res$status == "success_with_pred_errors") partial_success_count <- partial_success_count + 1
  else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
  else error_count <- error_count + 1
}

cat("--- LOGGING: Tweaking run: Overall Summary (", group_name, predictor_type_suffix, "Models) ---\n")
cat("Total Species Targets: 1 (tweaking run)\n")
cat("Fully Successful Runs:", success_count, "\n")
cat("Partially Successful Runs:", partial_success_count, "\n")
cat("Skipped:", skipped_count, "\n")
cat("Errors:", error_count, "\n")

if (length(occurrence_counts) > 0) {
  occ_count_df <- data.frame( Species = names(occurrence_counts), OccurrenceCountAfterThinning = unlist(occurrence_counts)) %>% dplyr::arrange(Species)
  occ_count_df_print <- occ_count_df; occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
  occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, "_tweaking.csv"))
  tryCatch({ readr::write_csv(occ_count_df, occ_count_file); cat("--- LOGGING: Tweaking run: Occurrence counts saved:", occ_count_file, "---\n")}, error = function(e) { cat("--- LOGGING: Tweaking run: ERROR - Failed save counts:", e$message, "---\n") })
  cat("\n--- Final Occurrence Counts (", group_name, predictor_type_suffix, ", after final thinning) ---\n"); print(occ_count_df_print)
} else { cat("--- LOGGING: Tweaking run: WARN - No occurrence counts were recorded.---\n") }

gc(full=TRUE)
cat("--- Script tweaking_base_run.R (with BIOMOD2 step) finished. ---\n")
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