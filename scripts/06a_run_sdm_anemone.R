# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using selected Env Predictors (SDMtune)
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (SDMtune) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
# Ensure SDMtune is loaded
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr)
# Source helpers
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R"); source(env_helper_path)
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R"); source(sdm_helper_path)

# --- 2. Define Group Specifics ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir

# --- 3. *** IMPORTANT: Define/Load FINAL Variable List *** ---
# This list (core names) MUST be defined after VIF analysis in 05a.
# Option 1: Define it directly here
final_core_vars_group <- c(
  "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_range",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  "chl_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
# Option 2: Load from config if defined there (preferred for consistency)
# if (!exists("final_vars_anemone", where = config)) {
#    stop("Object 'final_vars_anemone' not found in config. Define it after VIF analysis.")
# }
# final_core_vars_group <- config$final_vars_anemone
if(is.null(final_core_vars_group) || length(final_core_vars_group) < 2) {
  stop("Final core variable list for '", group_name, "' is missing or too short.")
}
cat("--- Using FINAL core variables for SDMs:", paste(final_core_vars_group, collapse=", "), "---\n")

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

# --- 6. Main Loop: Species -> Scenario ---
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
    
    # *** Load only the selected environmental predictors for this scenario ***
    predictor_stack <- load_selected_env_data(scenario, scenario_vars, config)
    if(is.null(predictor_stack)) {
      warning("Failed to load selected env data for scenario '", scenario, "'. Skipping.", call.=FALSE)
      total_sdms_skipped <- total_sdms_skipped + 1; next
    }
    cat("    Loaded predictor stack with layers:", paste(names(predictor_stack), collapse=", "), "\n")
    
    # Thin occurrences
    occ_sf_thinned <- thin_individual_occ(occ_sf_clean, predictor_stack, config)
    if(is.null(occ_sf_thinned) || nrow(occ_sf_thinned) < config$min_occurrences_sdm) {
      warning("Skipping: Not enough occs after thinning."); total_sdms_skipped <- total_sdms_skipped + 1; rm(predictor_stack); gc(); next
    }
    cat("    Thinned occurrences:", nrow(occ_sf_thinned), "points.\n")
    
    # Background Points
    background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config)
    if (is.null(background_points)) {
      warning("Failed background gen."); total_sdms_skipped <- total_sdms_skipped + 1; rm(predictor_stack, occ_sf_thinned); gc(); next
    }
    
    # Run SDM Tuning (using SDMtune helper)
    SDMtune_results <- run_sdm_SDMtune_grid(occ_sf_thinned, predictor_stack, background_points, config)
    if (is.null(SDMtune_results)) {
      warning("SDMtune::gridSearch failed."); total_sdms_skipped <- total_sdms_skipped + 1; rm(predictor_stack, occ_sf_thinned, background_points); gc(); next
    }
    
    # Save Results
    tryCatch(saveRDS(SDMtune_results, file = results_file), error=function(e){warning("Failed save SDMtune results object: ", e$message)})
    tryCatch(readr::write_csv(SDMtune::results(SDMtune_results), file = eval_file), error=function(e){warning("Failed save Eval table: ", e$message)})
    cat("    SDMtune results saved.\n")
    
    # Predict Best Model
    prediction_raster <- predict_sdm_SDMtune(SDMtune_results, predictor_stack, config)
    if (is.null(prediction_raster)) {
      warning("Prediction failed."); total_sdms_skipped <- total_sdms_skipped + 1; rm(predictor_stack, occ_sf_thinned, background_points, SDMtune_results); gc(); next
    }
    
    # Save Prediction & Model Object (Saving the whole SDMtune result object which contains the best model)
    tryCatch({
      terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
      cat("    Prediction raster saved to:", basename(pred_file), "\n")
      total_sdms_run <- total_sdms_run + 1
    }, error=function(e){warning("Failed save prediction raster: ", e$message); total_sdms_skipped <- total_sdms_skipped + 1})
    
    # Save the full SDMtune object - it contains the best model
    tryCatch({ saveRDS(SDMtune_results, file = model_obj_file); cat("    SDMtune object (containing best model) saved to:", basename(model_obj_file), "\n")
    }, error=function(e){warning("Failed save SDMtune object: ", e$message)})
    
    
    rm(predictor_stack, occ_sf_thinned, background_points, SDMtune_results, prediction_raster); gc()
  } # End scenario loop
  rm(occ_sf_clean); gc()
} # End species loop

cat("\n--- Script 06a finished. ---"); cat("Total Anemone SDMs run:", total_sdms_run, "\nSkipped:", total_sdms_skipped, "\n")
#-------------------------------------------------------------------------------