# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using PCA Predictors (SDMtune Workflow)
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (PCA + SDMtune) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr) # Ensure SDMtune is loaded
# Source helpers
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R"); source(env_helper_path)
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R"); source(sdm_helper_path)

# --- 2. Define Group Specifics ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
cat("--- Processing Group:", group_name, "---\n")

# --- 3. Load PCA Raster Paths ---
pca_paths_rds <- config$pca_raster_paths_rds_path # Get path from config
if (!file.exists(pca_paths_rds)) {
  stop("PCA raster paths file not found. Run script 05 (preprocessing) first: ", pca_paths_rds)
}
pca_raster_paths <- readRDS(pca_paths_rds)
if (!is.list(pca_raster_paths) || length(pca_raster_paths) == 0) {
  stop("PCA raster paths list is empty or invalid.")
}
cat("  Loaded PCA raster paths for scenarios:", paste(names(pca_raster_paths), collapse=", "), "\n")

# --- 4. Create Output Directories ---
group_pred_dir <- file.path(config$predictions_dir, group_name)
group_results_dir <- file.path(config$results_dir, group_name)
group_models_dir <- file.path(config$models_dir, group_name)
dir.create(group_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 5. Load Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) { stop("Failed load species list: ", e$message) })

# --- 6. Main Loop: Species -> Tuning -> Final Model -> Prediction ---
total_sdms_run <- 0; total_sdms_skipped <- 0; tuning_scenario <- "current" # Tune using current predictors

for (i in 1:nrow(species_df)) {
  species_name <- species_df$scientificName[i]
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_df$AphiaID[i]
  cat("\n---------------- Processing Species:", species_name, "(", species_aphia_id, ") ----------------\n")
  
  # Define species-specific file paths
  tuning_results_file <- file.path(group_results_dir, paste0("sdm_tuning_", species_name_sanitized, ".rds"))
  final_model_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, ".rds"))
  
  # Check if final model already exists to skip the whole species
  if (!config$force_rerun$run_standard_sdms && file.exists(final_model_file)) {
    cat("  Final model exists. Skipping tuning and training for this species.\n")
    # Still need to check predictions for all scenarios
  } else {
    
    # --- Data Loading & Prep for Tuning ---
    # Load occurrences (just coordinates needed)
    occs_coords <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config)
    if(is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm) {
      warning("Skipping: Not enough occurrences for AphiaID ", species_aphia_id, ".", call. = FALSE)
      total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next # Skip all scenarios for this species
    }
    
    # Load CURRENT PCA predictor stack for tuning
    tuning_pca_stack_path <- pca_raster_paths[[tuning_scenario]]
    if (is.null(tuning_pca_stack_path) || !file.exists(tuning_pca_stack_path)) {
      warning("PCA stack for tuning scenario '", tuning_scenario, "' not found. Skipping species.", call.=FALSE)
      total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
    }
    tuning_predictor_stack <- tryCatch(terra::rast(tuning_pca_stack_path), error = function(e) NULL)
    if(is.null(tuning_predictor_stack)) {
      warning("Failed load PCA stack for tuning scenario '", tuning_scenario, "'. Skipping species.", call.=FALSE)
      total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
    }
    
    # Generate background points using the tuning stack
    background_points <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config)
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
    rm(tuning_predictor_stack, background_points, tuning_output); gc() # Clean up tuning data
    
  } # End of check for existing final model
  
  # --- Prediction Loop (Run even if final model was loaded) ---
  if (!exists("final_model", inherits = FALSE)) { # Load final model if skipped training
    if(file.exists(final_model_file)) {
      cat("  Loading existing final model for prediction...\n")
      final_model <- readRDS(final_model_file)
      if(!inherits(final_model, "SDMmodel")){
        warning("Loaded object is not a valid SDMmodel. Cannot predict.", call.=FALSE)
        final_model <- NULL # Ensure it's null if loading failed
      }
    } else {
      warning("Final model file not found and training was skipped. Cannot predict.", call.=FALSE)
      final_model <- NULL
    }
  }
  
  if (!is.null(final_model)) {
    for (pred_scenario in names(pca_raster_paths)) { # Predict for all scenarios with PCA paths
      cat("    -- Predicting Scenario:", pred_scenario, "--\n")
      pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", pred_scenario, ".tif"))
      
      if (!config$force_rerun$run_standard_sdms && file.exists(pred_file)) {
        cat("      Prediction exists. Skipping.\n"); total_sdms_skipped <- total_sdms_skipped + 1; next
      }
      
      # Load PCA stack for the prediction scenario
      pred_pca_stack_path <- pca_raster_paths[[pred_scenario]]
      if (is.null(pred_pca_stack_path) || !file.exists(pred_pca_stack_path)) {
        warning("PCA stack for prediction scenario '", pred_scenario, "' not found. Skipping prediction.", call.=FALSE)
        total_sdms_skipped <- total_sdms_skipped + 1; next
      }
      pred_predictor_stack <- tryCatch(terra::rast(pred_pca_stack_path), error = function(e) NULL)
      if(is.null(pred_predictor_stack)) {
        warning("Failed load PCA stack for prediction scenario '", pred_scenario, "'. Skipping prediction.", call.=FALSE)
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
    total_sdms_skipped <- total_sdms_skipped + length(names(pca_raster_paths)) # Skipped all scenario predictions
  }
  
  # Clean up model object before next species
  if (exists("final_model", inherits = FALSE)) rm(final_model)
  if (exists("occs_coords", inherits = FALSE)) rm(occs_coords)
  gc()
} # End species loop

cat("\n--- Script 06a finished. ---"); cat("Total Anemone SDMs run/predicted:", total_sdms_run, "\nSkipped:", total_sdms_skipped, "\n")
#-------------------------------------------------------------------------------