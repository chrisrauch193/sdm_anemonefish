# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using either PCA or VIF-selected Env
# Predictors (SDMtune Workflow)
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr)
# Source helpers
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R"); source(env_helper_path)
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R"); source(sdm_helper_path)

# --- 2. Define Group Specifics & Predictor Type ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors
predictor_type_suffix <- ifelse(use_pca, "_pca", "_vif")
cat("--- Processing Group:", group_name, "---\n")
cat("--- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---\n")

# --- 3. Load Predictor Information ---
predictor_paths_or_list <- NULL
if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path # Get path from config
  if (!file.exists(pca_paths_rds)) {
    stop("PCA raster paths file not found. Run script 05b first: ", pca_paths_rds)
  }
  predictor_paths_or_list <- readRDS(pca_paths_rds)
  if (!is.list(predictor_paths_or_list) || length(predictor_paths_or_list) == 0) {
    stop("PCA raster paths list is empty or invalid.")
  }
  cat("  Loaded PCA raster paths for scenarios:", paste(names(predictor_paths_or_list), collapse=", "), "\n")
} else {
  predictor_paths_or_list <- config$final_vars_vif_anemone # Get VIF list
  if(is.null(predictor_paths_or_list) || length(predictor_paths_or_list) < 2) {
    stop("Final VIF variable list for '", group_name, "' (`final_vars_vif_anemone` in config) is missing or too short.")
  }
  cat("  Using VIF-selected core variables:", paste(predictor_paths_or_list, collapse=", "), "\n")
}

# --- 4. Create Output Directories (Include predictor type in name) ---
group_pred_dir <- file.path(config$predictions_dir, paste0(group_name, predictor_type_suffix))
group_results_dir <- file.path(config$results_dir, paste0(group_name, predictor_type_suffix))
group_models_dir <- file.path(config$models_dir, paste0(group_name, predictor_type_suffix))
dir.create(group_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 5. Load Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) { stop("Failed load species list: ", e$message) })

# --- 6. Main Loop: Species -> Tuning -> Final Model -> Prediction ---
total_sdms_run <- 0; total_sdms_skipped <- 0; tuning_scenario <- "current" # Always tune using current predictors

for (i in 1:nrow(species_df)) {
  species_name <- species_df$scientificName[i]
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_df$AphiaID[i]
  cat("\n---------------- Processing Species:", species_name, "(", species_aphia_id, ") ----------------\n")
  
  # Define species-specific file paths (include predictor type suffix)
  tuning_results_file <- file.path(group_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, ".rds"))
  final_model_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  # Check if final model already exists to potentially skip tuning/training
  if (!config$force_rerun$run_standard_sdms && file.exists(final_model_file)) {
    cat("  Final model exists. Skipping tuning and training for this species.\n")
    # Flag to load the existing model later for predictions
    load_existing_model <- TRUE
  } else {
    load_existing_model <- FALSE
    
    # --- Data Loading & Prep for Tuning ---
    occs_coords <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config)
    if(is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm) {
      warning("Skipping species: Not enough occurrences for AphiaID ", species_aphia_id, ".", call. = FALSE)
      total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
    }
    
    # Load predictors for TUNING (always 'current' scenario)
    tuning_predictor_stack <- NULL
    if(use_pca){
      tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
      if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) {
        warning("PCA stack for tuning scenario '", tuning_scenario, "' not found. Skipping species.", call.=FALSE)
        total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
      }
      tuning_predictor_stack <- tryCatch(terra::rast(tuning_predictor_path), error = function(e) NULL)
    } else {
      # Generate the specific list for the 'current' scenario from the VIF core list
      current_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, tuning_scenario, config)
      tuning_predictor_stack <- load_selected_env_data(tuning_scenario, current_vif_vars, config)
    }
    
    if(is.null(tuning_predictor_stack)) {
      warning("Failed load predictor stack for tuning scenario '", tuning_scenario, "'. Skipping species.", call.=FALSE)
      total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
    }
    
    # Background points
    background_points <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config, seed = i*100) # Vary seed slightly per species
    if (is.null(background_points)) {
      warning("Failed background gen for tuning. Skipping species.", call. = FALSE); total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
    }
    
    # --- Hyperparameter Tuning ---
    cat("  Tuning hyperparameters using scenario:", tuning_scenario, "\n")
    tuning_output <- run_sdm_tuning_kfold(occs_coords, tuning_predictor_stack, background_points, config, species_name)
    
    if (is.null(tuning_output) || is.null(tuning_output$best_hypers)) {
      warning("Hyperparameter tuning failed. Skipping species.", call. = FALSE); total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
    }
    best_hypers <- tuning_output$best_hypers
    tryCatch(saveRDS(tuning_output, tuning_results_file), error=function(e){warning("Failed save tuning results: ", e$message)})
    
    # --- Train Final Model ---
    cat("  Training final model using scenario:", tuning_scenario, "predictors and all data...\n")
    final_model <- train_final_sdm(occs_coords, tuning_predictor_stack, background_points, best_hypers, config, species_name)
    
    if(is.null(final_model)){
      warning("Final model training failed. Skipping species.", call.=FALSE); total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
    }
    tryCatch(saveRDS(final_model, final_model_file), error=function(e){warning("Failed save final model: ", e$message)})
    cat("  Final model saved:", basename(final_model_file), "\n")
    
    rm(tuning_predictor_stack, background_points, tuning_output, occs_coords); gc() # Clean up tuning data
    
  } # End of check for existing final model
  
  # --- Prediction Loop ---
  # Load the final model if it wasn't just trained
  if (load_existing_model) {
    if(file.exists(final_model_file)) {
      cat("  Loading existing final model for prediction...\n")
      final_model <- readRDS(final_model_file)
      if(!inherits(final_model, "SDMmodel")){
        warning("Loaded object is not a valid SDMmodel. Cannot predict.", call.=FALSE)
        final_model <- NULL
      }
    } else {
      warning("Final model file not found and training was skipped. Cannot predict.", call.=FALSE)
      final_model <- NULL
    }
  }
  
  # Proceed with predictions if the final model is available
  if (!is.null(final_model)) {
    # Predict for all scenarios for which we have predictor data
    scenarios_to_predict <- if(use_pca) names(predictor_paths_or_list) else config$env_scenarios
    
    for (pred_scenario in scenarios_to_predict) {
      cat("    -- Predicting Scenario:", pred_scenario, "--\n")
      pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", pred_scenario, predictor_type_suffix, ".tif"))
      
      if (!config$force_rerun$run_standard_sdms && file.exists(pred_file)) {
        cat("      Prediction exists. Skipping.\n"); total_sdms_skipped <- total_sdms_skipped + 1; next
      }
      
      # Load predictors for the prediction scenario
      pred_predictor_stack <- NULL
      if(use_pca) {
        pred_predictor_path <- predictor_paths_or_list[[pred_scenario]]
        if (is.null(pred_predictor_path) || !file.exists(pred_predictor_path)) {
          warning("PCA stack for prediction scenario '", pred_scenario, "' not found. Skipping prediction.", call.=FALSE)
          total_sdms_skipped <- total_sdms_skipped + 1; next
        }
        pred_predictor_stack <- tryCatch(terra::rast(pred_predictor_path), error = function(e) NULL)
      } else {
        # Generate the specific VIF var list for this scenario
        scenario_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, pred_scenario, config)
        if(length(scenario_vif_vars) < 1) { warning("No VIF vars generated for scenario '", pred_scenario, "'. Skipping pred.", call.=FALSE); total_sdms_skipped <- total_sdms_skipped + 1; next}
        pred_predictor_stack <- load_selected_env_data(pred_scenario, scenario_vif_vars, config)
      }
      
      if(is.null(pred_predictor_stack)) {
        warning("Failed load predictor stack for prediction scenario '", pred_scenario, "'. Skipping prediction.", call.=FALSE)
        total_sdms_skipped <- total_sdms_skipped + 1; next
      }
      
      # Predict suitability
      prediction_raster <- predict_sdm_suitability(final_model, pred_predictor_stack, config)
      
      if (!is.null(prediction_raster)) {
        tryCatch({
          terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
          cat("      Prediction raster saved:", basename(pred_file), "\n")
          total_sdms_run <- total_sdms_run + 1
        }, error=function(e){warning("Failed save prediction raster: ", e$message); total_sdms_skipped <- total_sdms_skipped + 1})
      } else {
        warning("Prediction failed for scenario '", pred_scenario, "'.", call. = FALSE); total_sdms_skipped <- total_sdms_skipped + 1
      }
      rm(pred_predictor_stack, prediction_raster); gc()
    } # End prediction scenario loop
  } else {
    cat("  Skipping prediction loop as final model is not available.\n")
    total_sdms_skipped <- total_sdms_skipped + length(scenarios_to_predict) # Skipped all scenario predictions
  }
  
  # Clean up model object before next species
  if (exists("final_model", inherits = FALSE)) rm(final_model)
  gc()
} # End species loop

cat("\n--- Script 06a finished. ---"); cat("Total Anemone SDMs run/predicted:", total_sdms_run, "\nSkipped:", total_sdms_skipped, "\n")
#-------------------------------------------------------------------------------