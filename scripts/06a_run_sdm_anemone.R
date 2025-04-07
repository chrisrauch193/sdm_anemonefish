# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using selected Env Predictors (sdmtune)
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (sdmtune) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, sdmtune, tools)
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R"); source(sdm_helper_path)

# --- 2. Define Group Specifics ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
# *** IMPORTANT: Load the FINAL selected variable list (core names) for this group ***
# This list should be determined *after* iterative VIF analysis in 05a.
final_core_vars_group <- config$final_vars_anemone # Assumes this list exists in config or is loaded
if(is.null(final_core_vars_group)) stop("Please define 'config$final_vars_anemone' with the final core variable list.")

# Create output directories if they don't exist
group_pred_dir <- file.path(config$predictions_dir, group_name)
group_results_dir <- file.path(config$results_dir, group_name)
group_models_dir <- file.path(config$models_dir, group_name)
dir.create(group_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3. Load Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) { stop("Failed load species list: ", e$message) })

# --- 4. Main Loop: Species -> Scenario ---
total_sdms_run <- 0; total_sdms_skipped <- 0
for (i in 1:nrow(species_df)) {
  species_name_sanitized <- gsub(" ", "_", species_df$scientificName[i])
  species_aphia_id <- species_df$AphiaID[i]
  cat("\n---------------- Processing Species:", species_df$scientificName[i], "----------------\n")
  
  # Load occurrences once per species
  occ_sf_clean <- load_clean_individual_occ(species_aphia_id, occurrence_dir, config)
  if(is.null(occ_sf_clean) || nrow(occ_sf_clean) < config$min_occurrences_sdm) {
    warning("Skipping: Not enough occurrences."); total_sdms_skipped <- total_sdms_skipped + length(config$env_scenarios); next
  }
  
  for (scenario in config$env_scenarios) {
    cat("  --- Scenario:", scenario, "---\n")
    
    pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", scenario, ".tif"))
    results_file <- file.path(group_results_dir, paste0("sdm_results_", species_name_sanitized, "_", scenario, ".rds"))
    eval_file <- file.path(group_results_dir, paste0("sdm_eval_", species_name_sanitized, "_", scenario, ".csv"))
    model_obj_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, "_", scenario, ".rds"))
    
    if (!config$force_rerun$run_standard_sdms && file.exists(pred_file)) {
      cat("    Prediction exists. Skipping.\n"); total_sdms_skipped <- total_sdms_skipped + 1; next
    }
    
    # Generate variable names for this scenario from the FINAL core list
    scenario_vars <- generate_scenario_variable_list(final_core_vars_group, scenario, config)
    if(length(scenario_vars) < 1) { warning("No variables generated for scenario. Skipping."); total_sdms_skipped <- total_sdms_skipped + 1; next }
    
    # Load only the selected environmental predictors for this scenario
    predictor_stack <- load_selected_env_data(scenario, scenario_vars, config); if(is.null(predictor_stack)) { total_sdms_skipped <- total_sdms_skipped + 1; next }
    
    # Thin occurrences
    occ_sf_thinned <- thin_individual_occ(occ_sf_clean, predictor_stack, config)
    if(is.null(occ_sf_thinned) || nrow(occ_sf_thinned) < config$min_occurrences_sdm) { warning("Skipping: Not enough occs after thinning."); total_sdms_skipped <- total_sdms_skipped + 1; next }
    
    # Background Points
    background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config); if (is.null(background_points)) { warning("Failed background gen."); total_sdms_skipped <- total_sdms_skipped + 1; next}
    
    # Run SDM Tuning (using sdmtune helper)
    sdmtune_results <- run_sdm_sdmtune_grid(occ_sf_thinned, predictor_stack, background_points, config)
    if (is.null(sdmtune_results)) { warning("sdmtune::gridSearch failed."); total_sdms_skipped <- total_sdms_skipped + 1; next }
    
    # Save Results
    tryCatch(saveRDS(sdmtune_results, file = results_file), error=function(e){warning("Failed save sdmtune results.")})
    tryCatch(readr::write_csv(sdmtune::results(sdmtune_results), file = eval_file), error=function(e){warning("Failed save Eval table.")})
    cat("    sdmtune results saved.\n")
    
    # Predict Best Model
    prediction_raster <- predict_sdm_sdmtune(sdmtune_results, predictor_stack, config)
    if (is.null(prediction_raster)) { warning("Prediction failed."); total_sdms_skipped <- total_sdms_skipped + 1; next }
    
    # Save Prediction & Model
    tryCatch({ terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES")); cat("    Prediction raster saved.\n"); total_sdms_run <- total_sdms_run + 1 }, error=function(e){warning("Failed save prediction."); total_sdms_skipped <- total_sdms_skipped + 1})
    # Save the best model object (sdmtune stores models differently, save the main object or extract best)
    # For simplicity, let's save the whole sdmtune result object again, or extract best model if needed later
    # tryCatch({ best_model_obj <- sdmtune::get_best_model(sdmtune_results); saveRDS(best_model_obj, file = model_obj_file); cat("    Best model object saved.\n") }, error=function(e){warning("Failed save best model.", e$message)})
    
    rm(predictor_stack, occ_sf_thinned, background_points, sdmtune_results, prediction_raster); gc()
  } # End scenario loop
  rm(occ_sf_clean); gc()
} # End species loop
} # End group check (only runs for 'anemone')

cat("\n--- Script 06a finished. ---"); cat("Total SDMs run:", total_sdms_run, "\nSkipped:", total_sdms_skipped, "\n")
#-------------------------------------------------------------------------------