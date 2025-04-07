# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using selected Env Predictors (SDMtune)
# PARALLELIZED VERSION over species with doSNOW progress bar.
# Saves thinned occurrence count.
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (SDMtune - Parallelized) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
# Load necessary packages, including parallel ones
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr,
               foreach, snow, doSNOW) # Use doSNOW for progress

# Source helpers
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R"); source(env_helper_path)
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R"); source(sdm_helper_path)

# --- 2. Define Group Specifics ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
cat("Analyzing Group:", group_name, "\n")

# --- 3. Load FINAL Selected Variable List ---
if (!exists("final_vars_anemone", where = config)) {
  stop("Object 'final_vars_anemone' not found in config. Define it after VIF analysis in 05a.")
}
final_core_vars_group <- config$final_vars_anemone
if(is.null(final_core_vars_group) || length(final_core_vars_group) < 1) {
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

# --- 6. Setup Parallel Backend (doSNOW) ---
n_cores_to_use <- 1
pb <- NULL # Progress bar object
progress <- NULL # Progress reporting function
cl <- NULL # Cluster object

if (config$use_parallel && nrow(species_df) > 1) { # Only parallelize if more than 1 species
  n_cores_detected <- parallel::detectCores()
  n_cores_to_use <- min(config$num_cores, n_cores_detected - 1, nrow(species_df))
  if (n_cores_to_use < 1) n_cores_to_use <- 1
  
  if (n_cores_to_use > 1) {
    cat("Setting up parallel backend with", n_cores_to_use, "cores (doSNOW)...\n")
    cl <- parallel::makeCluster(n_cores_to_use)
    doSNOW::registerDoSNOW(cl)
    
    # Setup progress bar
    pb <- txtProgressBar(max = nrow(species_df), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Export necessary variables/functions to the workers
    # Ensure all functions from helpers and the config object are exported
    parallel::clusterExport(cl, varlist = ls(all.names = TRUE), envir = .GlobalEnv)
    # Load packages on workers
    parallel::clusterEvalQ(cl, { library(terra); library(sf); library(dplyr); library(SDMtune); library(readr); library(tools); library(stringr) })
    
  } else {
    cat("Parallel disabled or only 1 core available/needed. Running sequentially.\n")
    opts <- NULL
    foreach::registerDoSEQ()
  }
} else {
  cat("Running sequentially (parallel disabled or only 1 species).\n")
  opts <- NULL
  foreach::registerDoSEQ()
}


# --- 7. Parallel Loop over Species (with Progress Options) ---
results_list <- foreach::foreach(
  i = 1:nrow(species_df),
  .packages = c("terra", "sf", "dplyr", "readr", "SDMtune", "tools", "stringr"), # Ensure needed packages are listed
  .errorhandling = 'pass', # Continue if one species fails, result will be error object
  .options.snow = opts      # Pass progress options to doSNOW
) %dopar% { # Use %dopar% for parallel, %do% for sequential
  
  # Inside the loop, everything relates to a single species (index i)
  # NOTE: Standard cat/print output from workers might be jumbled or lost in parallel.
  # Consider using proper logging packages (e.g., futile.logger, logger) for complex parallel logs.
  
  species_name_sanitized <- gsub(" ", "_", species_df$scientificName[i])
  species_aphia_id <- species_df$AphiaID[i]
  species_log_prefix <- paste0("[Species ", i, "/", nrow(species_df), " (", species_df$scientificName[i], ")] ") # Simplified prefix
  
  # --- Initialize results for this species ---
  species_run_count = 0
  species_skip_count = 0
  scenario_details_list = list()
  
  # --- Load Occurrences ---
  occ_sf_clean <- load_clean_individual_occ(species_aphia_id, occurrence_dir, config)
  if(is.null(occ_sf_clean) || nrow(occ_sf_clean) < config$min_occurrences_sdm) {
    warning(species_log_prefix, "Skipping: Not enough occurrences (", NROW(occ_sf_clean), ").", call.=FALSE)
    return(list(species=species_df$scientificName[i], status="skipped_no_occ", run=0, skipped=length(config$env_scenarios)))
  }
  
  # --- Loop through Scenarios for this species ---
  for (scenario in config$env_scenarios) {
    # Scenario-specific setup
    scenario_status <- "skipped_unknown"; n_thinned_for_scenario <- NA
    pred_file <- file.path(group_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", scenario, ".tif"))
    results_file <- file.path(group_results_dir, paste0("sdm_results_", species_name_sanitized, "_", scenario, ".rds"))
    eval_file <- file.path(group_results_dir, paste0("sdm_eval_", species_name_sanitized, "_", scenario, ".csv"))
    model_obj_file <- file.path(group_models_dir, paste0("sdm_model_", species_name_sanitized, "_", scenario, ".rds"))
    
    # --- Skip if exists ---
    if (!config$force_rerun$run_standard_sdms && file.exists(pred_file)) {
      cat(species_log_prefix, "Scenario ", scenario, ": Prediction exists. Skipping.\n"); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_exists"; next
    }
    
    # --- Prepare Predictors ---
    scenario_vars <- generate_scenario_variable_list(final_core_vars_group, scenario, config)
    if(length(scenario_vars) < 1) { warning(species_log_prefix, scenario, ": No variables generated. Skipping.", call.=FALSE); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_no_vars"; next }
    predictor_stack <- load_selected_env_data(scenario, scenario_vars, config)
    if(is.null(predictor_stack)) { warning(species_log_prefix, scenario, ": Failed load env data. Skipping.", call.=FALSE); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_load_fail"; next }
    
    # --- Thin Occurrences ---
    occ_sf_thinned <- thin_individual_occ(occ_sf_clean, predictor_stack, config)
    if(is.null(occ_sf_thinned) || nrow(occ_sf_thinned) < config$min_occurrences_sdm) {
      warning(species_log_prefix, scenario, ": Not enough occs after thinning (", NROW(occ_sf_thinned), "). Skipping.", call.=FALSE); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_thinning"; rm(predictor_stack); gc(); next
    }
    n_thinned_for_scenario <- nrow(occ_sf_thinned)
    
    # --- Background Points ---
    background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config)
    if (is.null(background_points)) {
      warning(species_log_prefix, scenario, ": Failed background gen. Skipping.", call.=FALSE); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_background"; rm(predictor_stack, occ_sf_thinned); gc(); next
    }
    
    # --- Run SDM Tuning ---
    sdmtune_results <- run_sdm_sdmtune_grid(occ_sf_thinned, predictor_stack, background_points, config)
    if (is.null(sdmtune_results)) {
      warning(species_log_prefix, scenario, ": SDMtune::gridSearch failed. Skipping.", call.=FALSE); species_skip_count <- species_skip_count + 1; scenario_status <- "failed_tuning"; rm(predictor_stack, occ_sf_thinned, background_points); gc(); next
    }
    
    # --- Save Results ---
    tryCatch(saveRDS(sdmtune_results, file = results_file), error=function(e){warning(species_log_prefix, scenario, ": Failed save sdmtune results object: ", e$message)})
    tryCatch({ results_df <- SDMtune::results(sdmtune_results); results_df$n_thinned_occurrences <- n_thinned_for_scenario; readr::write_csv(results_df, file = eval_file) }, error=function(e){warning(species_log_prefix, scenario, ": Failed save Eval table: ", e$message)})
    cat(species_log_prefix, scenario, ": SDMtune results saved.\n")
    
    # --- Predict Best Model ---
    prediction_raster <- predict_sdm_sdmtune(sdmtune_results, predictor_stack, config)
    if (is.null(prediction_raster)) {
      warning(species_log_prefix, scenario, ": Prediction failed. Skipping save.", call.=FALSE); species_skip_count <- species_skip_count + 1; scenario_status <- "failed_prediction"; rm(predictor_stack, occ_sf_thinned, background_points, sdmtune_results); gc(); next
    }
    
    # --- Save Prediction & Model ---
    save_success <- FALSE
    tryCatch({
      terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
      cat(species_log_prefix, scenario, ": Prediction raster saved.\n")
      species_sdms_run <- species_sdms_run + 1
      save_success <- TRUE
      scenario_status <- "completed"
    }, error=function(e){ warning(species_log_prefix, scenario, ": Failed save prediction: ", e$message); species_skip_count <- species_skip_count + 1; scenario_status <- "failed_save_pred" })
    
    if(save_success) {
      tryCatch({ saveRDS(sdmtune_results, file = model_obj_file); cat(species_log_prefix, scenario, ": SDMtune object saved.\n")
      }, error=function(e){warning(species_log_prefix, scenario, ": Failed save SDMtune object: ", e$message)})
    }
    
    # Store scenario result
    scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=n_thinned_for_scenario)
    
    rm(predictor_stack, occ_sf_thinned, background_points, sdmtune_results, prediction_raster); gc()
    
  } # End scenario loop for this species
  
  rm(occ_sf_clean); gc()
  
  # Return status for this species from the foreach loop
  list(species=species_df$scientificName[i], status="finished_species_loop", run=species_sdms_run, skipped=species_skip_count, scenario_details = scenario_details_list)
  
} # End foreach species loop

# --- 8. Close progress bar if used ---
if(!is.null(pb)) close(pb)

# --- 9. Stop Parallel Backend ---
if (config$use_parallel && exists("cl") && inherits(cl, "cluster")) {
  cat("\nStopping parallel backend...\n")
  tryCatch(parallel::stopCluster(cl), error = function(e) {warning("Error stopping cluster: ", e$message)})
}

# --- 10. Summarize Results ---
cat("\n=========================================================\n")
cat("Standard Anemone SDM runs finished.\n")
# Aggregate results from the list
total_sdms_run <- sum(sapply(results_list, function(x) if(is.list(x) && !is.null(x$run)) x$run else 0), na.rm = TRUE)
total_sdms_skipped <- sum(sapply(results_list, function(x) if(is.list(x) && !is.null(x$skipped)) x$skipped else 0), na.rm = TRUE)
cat("Total SDMs successfully run and prediction saved:", total_sdms_run, "\n")
cat("Total SDM runs skipped (already exist, error, or few occs):", total_sdms_skipped, "\n")
# Optional: Print more details on errors/skips by examining results_list
cat("=========================================================\n")

cat("\n--- Script 06a finished. ---\n")
#-------------------------------------------------------------------------------