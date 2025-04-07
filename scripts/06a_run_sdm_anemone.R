# scripts/06a_run_sdm_anemone.R
#-------------------------------------------------------------------------------
# Run Standard SDMs for Anemone Species using selected Env Predictors (SDMtune)
# PARALLELIZED VERSION over species with future/progressr.
# Logs detailed output for EACH SPECIES to a separate file.
# Saves outputs to species-specific directories.
# Saves thinned occurrence count.
#-------------------------------------------------------------------------------
cat("--- Running Script 06a: Run Standard Anemone SDMs (SDMtune - Parallelized) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
# Load necessary packages
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr,
               future, future.apply, progressr, logger) # Added logger

# Source helpers
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R"); source(env_helper_path)
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R"); source(sdm_helper_path)
# No global logging setup needed here

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

# --- 4. Define BASE Output Directories (Species subdirs created in loop) ---
base_pred_dir <- file.path(config$predictions_dir, group_name)
base_results_dir <- file.path(config$results_dir, group_name)
base_models_dir <- file.path(config$models_dir, group_name)
base_log_dir <- file.path(config$log_dir, paste0("species_logs_", group_name)) # Specific dir for these logs
dir.create(base_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(base_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(base_models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(base_log_dir, recursive = TRUE, showWarnings = FALSE)
cat("Base output directories ensured. Species-specific logs will be in:", base_log_dir, "\n")

# --- 5. Load Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
}, error = function(e) { stop("Failed load species list: ", e$message) })
cat("Loaded", nrow(species_df), "species for group '", group_name, "'.\n")

# --- 6. Setup Parallel Backend (Future) & Progress Reporting ---
n_cores_to_use <- 1
if (config$use_parallel && nrow(species_df) > 1) {
  n_cores_detected <- future::availableCores()
  n_cores_to_use <- min(config$num_cores, n_cores_detected - 1, nrow(species_df))
  if (n_cores_to_use < 1) n_cores_to_use <- 1
  if (n_cores_to_use > 1) {
    cat("Setting up parallel backend with", n_cores_to_use, "cores (future::multisession)...\n")
    future::plan(future::multisession, workers = n_cores_to_use)
  } else {
    cat("Parallel disabled or only 1 core available/needed. Running sequentially.\n")
    future::plan(future::sequential)
  }
} else {
  cat("Running sequentially (parallel disabled or only 1 species).\n")
  future::plan(future::sequential)
}

# Setup progressr handlers
progressr::handlers(global = TRUE)
if (interactive() && requireNamespace("progress", quietly = TRUE)) {
  progressr::handlers("progress")
} else {
  progressr::handlers("txtprogressbar")
}

# --- 7. Parallel Loop over Species (using future_lapply with progressr) ---
cat("Starting SDM processing loop...\n")
results_list <- progressr::with_progress({
  
  p <- progressr::progressor(steps = nrow(species_df)) # Initialize progressor
  
  future.apply::future_lapply(1:nrow(species_df), function(i) {
    # --- Worker Setup ---
    suppressPackageStartupMessages({
      library(terra); library(sf); library(dplyr); library(readr);
      library(SDMtune); library(tools); library(stringr); library(logger)
    })
    # Define config access within future (might inherit, but safer to be explicit if needed)
    # current_config <- config
    
    # --- Species Info ---
    species_name_sci <- species_df$scientificName[i]
    species_name_sanitized <- gsub(" ", "_", species_name_sci)
    species_aphia_id <- species_df$AphiaID[i]
    
    # --- Species-Specific Directories & Logging ---
    species_pred_dir <- file.path(base_pred_dir, species_name_sanitized)
    species_results_dir <- file.path(base_results_dir, species_name_sanitized)
    species_models_dir <- file.path(base_models_dir, species_name_sanitized)
    species_log_file <- file.path(base_log_dir, paste0("sdm_log_", species_name_sanitized, ".log"))
    dir.create(species_pred_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(species_results_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(species_models_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Configure logger for THIS species run (appends to file)
    log_threshold(INFO)
    log_layout(layout_glue('{time} [{level}] {msg}')) # Simple layout for file
    log_appender(appender_file(species_log_file))
    
    log_info("--- Starting Species: {species_name_sci} (ID: {species_aphia_id}) ---")
    
    # --- Initialize results for this species ---
    species_run_count = 0; species_skip_count = 0; scenario_details_list = list()
    
    # --- Load Occurrences ---
    occ_sf_clean <- load_clean_individual_occ(species_aphia_id, occurrence_dir, config)
    if(is.null(occ_sf_clean) || nrow(occ_sf_clean) < config$min_occurrences_sdm) {
      log_warn("Skipping: Not enough occurrences ({NROW(occ_sf_clean)}).")
      p(amount = 1) # Increment progress bar even if skipped early
      return(list(species=species_name_sci, status="skipped_no_occ", run=0, skipped=length(config$env_scenarios)))
    }
    log_info("Loaded {nrow(occ_sf_clean)} initial occurrences.")
    
    # --- Loop through Scenarios ---
    for (scenario in config$env_scenarios) {
      log_info("--- Scenario: {scenario} ---")
      scenario_status <- "skipped_unknown"; n_thinned_for_scenario <- NA
      
      # --- Define Output Paths ---
      pred_file <- file.path(species_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", scenario, ".tif"))
      results_file <- file.path(species_results_dir, paste0("sdm_results_", species_name_sanitized, "_", scenario, ".rds"))
      eval_file <- file.path(species_results_dir, paste0("sdm_eval_", species_name_sanitized, "_", scenario, ".csv"))
      model_obj_file <- file.path(species_models_dir, paste0("sdm_model_", species_name_sanitized, "_", scenario, ".rds"))
      
      # --- Skip Check ---
      if (!config$force_rerun$run_standard_sdms && file.exists(pred_file)) {
        log_info("Prediction exists. Skipping."); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_exists"; scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=NA); next
      }
      
      # --- Prepare Predictors ---
      scenario_vars <- generate_scenario_variable_list(final_core_vars_group, scenario, config)
      if(length(scenario_vars) < 1) { log_warn("No variables generated. Skipping."); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_no_vars"; scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=NA); next }
      predictor_stack <- load_selected_env_data(scenario, scenario_vars, config)
      if(is.null(predictor_stack)) { log_warn("Failed load env data. Skipping."); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_load_fail"; scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=NA); next }
      log_info("Loaded predictor stack with {terra::nlyr(predictor_stack)} layers.")
      
      # --- Thin Occurrences ---
      occ_sf_thinned <- thin_individual_occ(occ_sf_clean, predictor_stack, config)
      if(is.null(occ_sf_thinned) || nrow(occ_sf_thinned) < config$min_occurrences_sdm) {
        log_warn("Not enough occs after thinning ({NROW(occ_sf_thinned)}). Skipping."); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_thinning"; rm(predictor_stack); gc(); scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=NROW(occ_sf_thinned)); next
      }
      n_thinned_for_scenario <- nrow(occ_sf_thinned)
      log_info("Thinned occurrences: {n_thinned_for_scenario} points.")
      
      # --- Background Points ---
      background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config) # Uses internal log/warn
      if (is.null(background_points)) {
        log_warn("Failed background gen. Skipping."); species_skip_count <- species_skip_count + 1; scenario_status <- "skipped_background"; rm(predictor_stack, occ_sf_thinned); gc(); scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=n_thinned_for_scenario); next
      }
      
      # --- Run SDM Tuning ---
      sdmtune_results <- run_sdm_sdmtune_grid(occ_sf_thinned, predictor_stack, background_points, config) # Uses internal log/warn
      if (is.null(sdmtune_results)) {
        log_warn("SDMtune::gridSearch failed. Skipping."); species_skip_count <- species_skip_count + 1; scenario_status <- "failed_tuning"; rm(predictor_stack, occ_sf_thinned, background_points); gc(); scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=n_thinned_for_scenario); next
      }
      
      # --- Save Results ---
      tryCatch(saveRDS(sdmtune_results, file = results_file), error=function(e){log_warn("Failed save sdmtune results object: {e$message}")})
      tryCatch({ results_df <- SDMtune::results(sdmtune_results); results_df$n_thinned_occurrences <- n_thinned_for_scenario; readr::write_csv(results_df, file = eval_file) }, error=function(e){log_warn("Failed save Eval table: {e$message}")})
      log_info("SDMtune results and evaluation table saved.")
      
      # --- Predict Best Model ---
      prediction_raster <- predict_sdm_sdmtune(sdmtune_results, predictor_stack, config) # Uses internal log/warn
      if (is.null(prediction_raster)) {
        log_warn("Prediction failed. Skipping save."); species_skip_count <- species_skip_count + 1; scenario_status <- "failed_prediction"; rm(predictor_stack, occ_sf_thinned, background_points, sdmtune_results); gc(); scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=n_thinned_for_scenario); next
      }
      
      # --- Save Prediction & Model ---
      save_success <- FALSE
      tryCatch({
        terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
        log_info("Prediction raster saved.")
        species_run_count <- species_run_count + 1
        save_success <- TRUE
        scenario_status <- "completed"
      }, error=function(e){ log_error("Failed save prediction: {e$message}"); species_skip_count <- species_skip_count + 1; scenario_status <- "failed_save_pred" })
      
      if(save_success) {
        tryCatch({ saveRDS(sdmtune_results, file = model_obj_file); log_info("SDMtune object saved.")
        }, error=function(e){log_warn("Failed save SDMtune object: {e$message}")})
      }
      
      scenario_details_list[[scenario]] <- list(status=scenario_status, n_thinned=n_thinned_for_scenario)
      rm(predictor_stack, occ_sf_thinned, background_points, sdmtune_results, prediction_raster); gc()
      
    } # End scenario loop
    
    rm(occ_sf_clean); gc()
    
    # Increment overall progress bar when a species finishes all its scenarios
    p()
    
    # Return results for this species
    list(species=species_name_sci, status="finished_species_loop", run=species_run_count, skipped=species_skip_count, scenario_details = scenario_details_list)
    
  }, future.seed = TRUE) # End future_lapply
  
}) # End with_progress

# --- 8. Stop Parallel Backend ---
current_plan <- class(future::plan())[1]
if (!current_plan %in% c("sequential", "uniprocess")) {
  log_info("Shutting down parallel workers ({current_plan})...")
  future::plan(future::sequential)
}

# --- 9. Summarize Results ---
cat("\n=========================================================\n")
log_info("Standard Anemone SDM runs finished.")
total_sdms_run <- sum(sapply(results_list, function(x) if(is.list(x) && !is.null(x$run)) x$run else 0), na.rm = TRUE)
total_sdms_skipped <- sum(sapply(results_list, function(x) if(is.list(x) && !is.null(x$skipped)) x$skipped else 0), na.rm = TRUE)
log_info("Total SDMs successfully run and prediction saved: {total_sdms_run}")
log_info("Total SDM runs skipped (already exist, error, or few occs): {total_sdms_skipped}")
errors <- results_list[sapply(results_list, inherits, "error")]
if (length(errors) > 0) {
  log_error("{length(errors)} species tasks failed entirely in the parallel loop.")
  for(err_idx in 1:min(5, length(errors))) { # Log first 5 errors
    log_error("Error for task {which(sapply(results_list, inherits, 'error'))[err_idx]}: {conditionMessage(errors[[err_idx]])}")
  }
}
cat("=========================================================\n")

cat("\n--- Script 06a finished. Check species log files in: ", base_log_dir, " ---\n")
#-------------------------------------------------------------------------------