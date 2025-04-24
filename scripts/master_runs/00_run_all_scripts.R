# scripts/00_run_all_scripts.R
#-------------------------------------------------------------------------------
# Master script to run the Anemone/Anemonefish SDM workflow (Steps 1-7)
# 1. Setup: Load config, install/load packages.
# 2. Data Acquisition: Get species lists, download occurrences, download env data.
# 3. Preprocessing: VIF/PCA analysis (requires manual variable selection input).
# 4. Standard SDM Runs: Execute SDMs for each species/scenario using PCA predictors.
# 5. Biotic SDM Runs: Execute Biotic-Only and Combined models for anemonefish.
#-------------------------------------------------------------------------------

# --- Clean Workspace ---
rm(list = ls())
gc() # Garbage collection

# --- 1. Load Configuration & Requirements ---
cat("\n--- Step 1: Loading Configuration and Requirements ---\n")
if (file.exists("scripts/config.R")) {
  source("scripts/config.R")
} else {
  stop("Configuration file 'scripts/config.R' not found.")
}
if (file.exists(file.path(config$scripts_dir, "01_install_requirements.R"))) {
  source(file.path(config$scripts_dir, "01_install_requirements.R"))
} else {
  stop("Requirements file 'scripts/01_install_requirements.R' not found.")
}
cat("--- Step 1: Finished ---\n")

# --- Helper Function for Output Check ---
check_output_exists <- function(files) {
  if (length(files) == 0) return(FALSE)
  all(file.exists(files))
}
check_any_output_exists <- function(directory, pattern) {
  if(!dir.exists(directory)) return(FALSE)
  length(list.files(directory, pattern = pattern, recursive = TRUE)) > 0
}


# --- 2a. Get Species Metadata ---
cat("\n--- Step 2a: Get Species Metadata ---\n")
output_files_02 <- c(config$anemone_species_list_file, config$anemonefish_species_list_file)
script_path_02 <- file.path(config$scripts_dir, "02_get_species_metadata.R")
if (!file.exists(script_path_02)) stop("Script not found: ", script_path_02)
if (config$force_rerun$get_species || !check_output_exists(output_files_02)) {
  cat("Running 02_get_species_metadata.R...\n")
  tryCatch({ source(script_path_02, echo = TRUE); if (!check_output_exists(output_files_02)) stop("Script 02 failed."); cat("Finished 02.\n")}, error = function(e) {stop("Error in 02: ", e$message)})
} else {cat("Skipping 02.\n")}
cat("--- Step 2a: Finished ---\n")


# --- 2b. Download Occurrence Data ---
cat("\n--- Step 2b: Download Occurrence Data ---\n")
output_exists_03 <- check_any_output_exists(config$anemone_occurrence_dir, "\\.csv$") && check_any_output_exists(config$anemonefish_occurrence_dir, "\\.csv$")
script_path_03 <- file.path(config$scripts_dir, "03_download_occurrence_data.R")
if (!file.exists(script_path_03)) stop("Script not found: ", script_path_03)
if (config$force_rerun$download_occurrences || !output_exists_03) {
  cat("Running 03_download_occurrence_data.R...\n")
  tryCatch({ source(script_path_03, echo = TRUE); cat("Finished 03.\n")}, error = function(e) {stop("Error in 03: ", e$message)})
} else {cat("Skipping 03.\n")}
cat("--- Step 2b: Finished ---\n")


# --- 2c. Download Environmental Data ---
cat("\n--- Step 2c: Download Environmental Data ---\n")
output_exists_04 <- dir.exists(config$scenario_folder_map$current) && dir.exists(config$terrain_folder) && file.exists(file.path(config$terrain_folder, "distcoast.tif"))
script_path_04 <- file.path(config$scripts_dir, "04_download_env_data.R")
if (!file.exists(script_path_04)) stop("Script not found: ", script_path_04)
if (config$force_rerun$download_env || !output_exists_04) {
  cat("Running 04_download_env_data.R...\n")
  tryCatch({ source(script_path_04, echo = TRUE); cat("Finished 04.\n")}, error = function(e) {stop("Error in 04: ", e$message)})
} else { cat("Skipping 04.\n") }
cat("--- Step 2c: Finished ---\n")


# --- 3. Preprocess Env & Occurrences (VIF/PCA) ---
cat("\n--- Step 3: Preprocess Env & Occurrences (VIF/PCA) ---\n")
output_files_05 <- c(config$selected_vars_rds_path, config$pca_raster_paths_rds_path, config$pca_models_rds_path)
script_path_05 <- file.path(config$scripts_dir, "05_preprocess_env_occurrence.R")
if (!file.exists(script_path_05)) stop("Script not found: ", script_path_05)
if (config$force_rerun$preprocess_env_occurrence || !check_output_exists(output_files_05)) {
  cat("Running 05_preprocess_env_occurrence.R...\n")
  cat("**********************************************************************\n")
  cat("*** NOTE: Script 05 requires MANUAL INTERVENTION! ***\n")
  cat("1. It will first run VIF/Correlation analysis for each group/scenario.\n")
  cat("2. It will save plots and suggested variables to:", config$log_dir, "\n")
  cat("3. REVIEW these outputs carefully.\n")
  cat("4. DEFINE your final variable selections for PCA.\n")
  cat("   (Edit script 05 directly or manage the 'selected_variables_for_pca.rds' file).\n")
  cat("5. RE-RUN this master script (or script 05 directly).\n")
  cat("   Script 05 will then detect your selections and proceed with PCA & Projection.\n")
  cat("**********************************************************************\n")
  Sys.sleep(3) # Pause for 3 seconds
  tryCatch({ source(script_path_05, echo = TRUE); cat("Finished attempt of 05.\n") }, error = function(e) { stop("Error in 05: ", e$message) })
  if (!check_output_exists(output_files_05)) { warning("Script 05 finished, but final output RDS files missing. Manual step likely needed.") }
} else { cat("Skipping 05.\n") }
# --- Critical Check for next step ---
if (!check_output_exists(output_files_05)) { stop("Cannot proceed. Output files from script 05 (VIF/PCA) are missing.") } else { cat("Verified outputs from script 05 exist.\n") }
cat("--- Step 3: Finished ---\n")


# --- 4. Run Standard SDMs ---
cat("\n--- Step 4: Run Standard SDMs ---\n")
output_exists_06 <- check_any_output_exists(config$predictions_dir, "^sdm_prediction_.*\\.tif$") && !check_any_output_exists(config$predictions_dir, "^sdm_prediction_biotic.*\\.tif$") # Check standard, exclude biotic
script_path_06 <- file.path(config$scripts_dir, "06_run_standard_sdms.R")
if (!file.exists(script_path_06)) stop("Script not found: ", script_path_06)
if (config$force_rerun$run_standard_sdms || !output_exists_06) {
  cat("Running 06_run_standard_sdms.R...\n")
  tryCatch({ source(script_path_06, echo = TRUE); cat("Finished 06.\n") }, error = function(e) { stop("Error in 06: ", e$message) })
} else { cat("Skipping 06.\n") }
# --- Critical Check for next step ---
if (!check_any_output_exists(config$predictions_dir, "^sdm_prediction_anemone.*\\.tif$")) {
  stop("Cannot proceed to Biotic SDMs. Standard SDM prediction files for ANEMONES are missing from: ", config$predictions_dir)
} else { cat("Verified standard anemone prediction outputs from script 06 exist.\n")}
cat("--- Step 4: Finished ---\n")


# --- 5. Run Anemonefish Biotic SDMs ---
cat("\n--- Step 5: Run Anemonefish Biotic SDMs ---\n")
# Check for *any* biotic prediction file
output_exists_07 <- check_any_output_exists(config$predictions_dir, "^sdm_prediction_biotic.*\\.tif$")
script_path_07 <- file.path(config$scripts_dir, "07_run_anemonefish_biotic_sdms.R")
if (!file.exists(script_path_07)) stop("Script not found: ", script_path_07)

if (config$force_rerun$run_biotic_sdms || !output_exists_07) {
  cat("Running 07_run_anemonefish_biotic_sdms.R...\n")
  # Add check for association file existence before running
  if(!file.exists(config$anemone_fish_association_file)) stop("Association file missing, cannot run biotic SDMs: ", config$anemone_fish_association_file)
  tryCatch({
    source(script_path_07, echo = TRUE)
    cat("Finished 07.\n")
  }, error = function(e) {
    stop("Error running 07_run_anemonefish_biotic_sdms.R: ", e$message)
  })
} else {
  cat("Skipping 07_run_anemonefish_biotic_sdms.R (Biotic prediction files exist).\n")
}
cat("--- Step 5: Finished ---\n")


cat("\n--- Master script finished execution for steps 1-5 ---\n")
#-------------------------------------------------------------------------------