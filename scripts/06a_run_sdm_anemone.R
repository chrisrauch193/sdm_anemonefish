# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using selected Env Predictors (SDMtune)
# PARALLELIZED VERSION over species. Saves thinned occurrence count.
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (SDMtune - Parallelized) ---\n")

# --- 1. Setup ---
if (!exists("config")) {
  source("scripts/config.R")
  if (!exists("config")) stop("Failed load config.")
}
# Load necessary packages
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr,
               foreach, doParallel) # Added foreach, doParallel

# Source helpers (ensure all needed functions are present)
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R")
if (!file.exists(env_helper_path)) stop("Env Helper missing.")
if (!file.exists(sdm_helper_path)) stop("SDM Helper missing.")
source(env_helper_path)
source(sdm_helper_path)

# --- 2. Define Group Specifics ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir

# --- 3. Load FINAL Selected Variable List ---
# Load from config (ensure it was defined there after VIF analysis)
if (!exists("final_vars_anemone", where = config)) {
  stop("Object 'final_vars_anemone' not found in config. Define it after VIF analysis in 05a.")
}
final_core_vars_group <- config$final_vars_anemone
if(is.null(final_core_vars_group) || length(final_core_vars_group) < 1) { # Need at least 1 predictor
  stop("Final core variable list for '", group_name, "' is missing or empty.")
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

# --- 6. Setup Parallel Backend ---
n_cores_to_use <- 1 # Default to sequential
if (config$use_parallel) {
  n_cores_detected <- parallel::detectCores()
  n_cores_to_use <- min(config$num_cores, n_cores_detected - 1, nrow(species_df))
  if (n_cores_to_use < 1) n_cores_to_use <- 1
  
  if (n_cores_to_use > 1) {
    cat("Setting up parallel backend with", n_cores_to_use, "cores...\n")
    cl <- parallel::makeCluster(n_cores_to_use)
    doParallel::registerDoParallel(cl)
    # Export necessary variables/functions to the workers
    # Include all functions loaded from helpers AND the config object
    parallel::clusterExport(cl, varlist = c("config", "load_clean_individual_occ",
                                            "thin_individual_occ", "generate_sdm_background",
                                            "run_sdm_sdmtune_grid", "predict_sdm_sdmtune",
                                            "load_selected_env_data", "generate_scenario_variable_list",
                                            "preprocess_env_rasters", "load_stack_env_data")) # Added load_stack
    # Load packages on workers
    parallel::clusterEvalQ(cl, { library(terra); library(sf); library(dplyr); library(SDMtune); library(readr); library(tools); library(stringr) })
  } else {
    cat("Parallel disabled or only 1 core available/needed. Running sequentially.\n")
    foreach::registerDoSEQ() # Explicitly register sequential backend
  }
} else {
  cat("Running sequentially.\n")
  foreach::registerDoSEQ()
}


# --- 7. Parallel Loop over Species ---
# Use foreach loop
results_list <- foreach::foreach(
  i = 1:nrow(species_df),
  .errorhandling = 'pass' # Continue if one species fails
) %dopar% { # Use %dopar% for parallel, %do% for sequential
  
  # Inside the loop, everything relates to a single species (index i)
  # Load required packages *inside* the dopar loop for safety
  require(terra)
  require(sf)
  require(dplyr)
  require(readr)
  require(SDMtune)
  require(tools)
  require(stringr)
  
  species_name_sanitized <- gsub(" ", "_", species_df$scientificName[i])
  species_aphia_id <- species_df$AphiaID[i]
  # Include process ID in log messages for easier debugging
  log_prefix <- paste0("[Worker ", Sys.getpid(), " | Species ", i, "/", nrow(species_df), " (", species_df$scientificName[i], ")] ")
  cat(log_prefix, "Starting processing...\n")
  
  # Initialize counts for this species
  species_sdms_run <- 0
  species_sdms_skipped <- 0
  species_results <- list() # To store results per scenario for this species
  
  # Load occurrences once per species
  occ_sf_clean <- load_clean_individual_occ(species_aphia_id, occurrence_dir, config)
  if(is.null(occ_sf_clean) || nrow(occ_sf_clean) < config$min_occurrences_sdm) {
    cat(log_prefix, "Skipping: Not enough occurrences (", NROW(occ_sf_clean), ").\n")
    return(list(species=species_df$scientificName[i], status="skipped_no_occ", run=0, skipped=length(config$env_scenarios)))
  }
  cat(log_prefix, "Loaded", nrow(occ_sf_clean), "initial occurrences.\n")
  
  # Loop through scenarios for THIS species
  for (scenario in config$env_scenarios) {
    cat(log_prefix, "--- Scenario:", scenario, "---\n")
    
    scenario_status <- "skipped_unknown" # Default status for this scenario run
    n_thinned_for_scenario <- NA # Default value
    
    # --- Define Output Paths ---
    pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", scenario, ".tif"))
    results_file <- file.path(group_results_dir, paste0("sdm_results_", species_name_sanitized, "_", scenario, ".rds"))
    eval_file <- file.path(group_results_dir, paste0("sdm_eval_", species_name_sanitized, "_", scenario, ".csv"))
    model_obj_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, "_", scenario, ".rds"))
    
    # --- Check Existence ---
    if (!config$force_rerun$run_standard_sdms && file.exists(pred_file)) {
      cat(log_prefix, "  Prediction exists. Skipping.\n"); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "skipped_exists"; next
    }
    
    # --- Prepare Predictors ---
    scenario_vars <- generate_scenario_variable_list(final_core_vars_group, scenario, config)
    if(length(scenario_vars) < 1) { warning(log_prefix, "No variables generated for scenario. Skipping."); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "skipped_no_vars"; next }
    
    predictor_stack <- load_selected_env_data(scenario, scenario_vars, config)
    if(is.null(predictor_stack)) { warning(log_prefix, "Failed to load selected env data. Skipping."); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "skipped_load_fail"; next }
    cat(log_prefix, "  Loaded predictor stack with layers:", paste(names(predictor_stack), collapse=", "), "\n")
    
    # --- Thin Occurrences ---
    occ_sf_thinned <- thin_individual_occ(occ_sf_clean, predictor_stack, config)
    if(is.null(occ_sf_thinned) || nrow(occ_sf_thinned) < config$min_occurrences_sdm) {
      warning(log_prefix, "Skipping: Not enough occs after thinning (", NROW(occ_sf_thinned), ")."); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "skipped_thinning"; rm(predictor_stack); gc(); next
    }
    n_thinned_for_scenario <- nrow(occ_sf_thinned) # STORE THE COUNT
    cat(log_prefix, "  Thinned occurrences:", n_thinned_for_scenario, "points.\n")
    
    # --- Background Points ---
    background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config)
    if (is.null(background_points)) {
      warning(log_prefix, "Failed background gen."); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "skipped_background"; rm(predictor_stack, occ_sf_thinned); gc(); next
    }
    
    # --- Run SDM Tuning ---
    sdmtune_results <- run_sdm_sdmtune_grid(occ_sf_thinned, predictor_stack, background_points, config)
    if (is.null(sdmtune_results)) {
      warning(log_prefix, "SDMtune::gridSearch failed."); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "failed_tuning"; rm(predictor_stack, occ_sf_thinned, background_points); gc(); next
    }
    
    # --- Save Results ---
    tryCatch(saveRDS(sdmtune_results, file = results_file), error=function(e){warning(log_prefix, "Failed save sdmtune results object: ", e$message)})
    tryCatch({
      results_df <- SDMtune::results(sdmtune_results)
      results_df$n_thinned_occurrences <- n_thinned_for_scenario # Add count here
      readr::write_csv(results_df, file = eval_file)
    }, error=function(e){warning(log_prefix, "Failed save Eval table: ", e$message)})
    cat(log_prefix, "  SDMtune results saved.\n")
    
    # --- Predict Best Model ---
    prediction_raster <- predict_sdm_sdmtune(sdmtune_results, predictor_stack, config)
    if (is.null(prediction_raster)) {
      warning(log_prefix, "Prediction failed."); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "failed_prediction"; rm(predictor_stack, occ_sf_thinned, background_points, sdmtune_results); gc(); next
    }
    
    # --- Save Prediction & Model ---
    save_success <- FALSE
    tryCatch({
      terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
      cat(log_prefix, "  Prediction raster saved to:", basename(pred_file), "\n")
      species_sdms_run <- species_sdms_run + 1
      save_success <- TRUE
      scenario_status <- "completed"
    }, error=function(e){ warning(log_prefix, "Failed save prediction raster: ", e$message); species_sdms_skipped <- species_sdms_skipped + 1; scenario_status <- "failed_save_pred" })
    
    # Save the full SDMtune object (contains best model info/object)
    if(save_success) {
      tryCatch({ saveRDS(sdmtune_results, file = model_obj_file); cat(log_prefix, "  SDMtune object saved to:", basename(model_obj_file), "\n")
      }, error=function(e){warning(log_prefix, "Failed save SDMtune object: ", e$message)})
    }
    
    # Store scenario result
    species_results[[scenario]] <- list(status=scenario_status, n_thinned=n_thinned_for_scenario)
    
    rm(predictor_stack, occ_sf_thinned, background_points, sdmtune_results, prediction_raster); gc()
    
  } # End scenario loop
  
  rm(occ_sf_clean); gc()
  
  # Return status for this species from the foreach loop
  list(species=species_df$scientificName[i], status="finished_species_loop", run=species_sdms_run, skipped=species_sdms_skipped, scenario_details = species_results)
  
} # End foreach species loop

# --- 7. Stop Parallel Backend ---
if (config$use_parallel && exists("cl")) {
  cat("\nStopping parallel backend...\n")
  tryCatch(parallel::stopCluster(cl), error = function(e) {warning("Error stopping cluster: ", e$message)})
}

# --- 8. Summarize Results ---
cat("\n=========================================================\n")
cat("Standard Anemone SDM runs finished.\n")
# Aggregate results from the list returned by foreach
total_sdms_run <- sum(sapply(results_list, function(x) if(is.list(x) && !is.null(x$run)) x$run else 0), na.rm = TRUE)
total_sdms_skipped <- sum(sapply(results_list, function(x) if(is.list(x) && !is.null(x$skipped)) x$skipped else 0), na.rm = TRUE)
# Add more detailed summary if needed by parsing results_list further
cat("Total SDMs successfully run and saved:", total_sdms_run, "\n")
cat("Total SDMs skipped (due to errors, missing data, or existing files):", total_sdms_skipped, "\n")
# You could print details about which species/scenarios failed here by checking results_list
cat("=========================================================\n")

cat("\n--- Script 06a finished. ---\n")
#-------------------------------------------------------------------------------