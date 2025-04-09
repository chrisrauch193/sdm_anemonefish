# scripts/config.R
#-------------------------------------------------------------------------------
# Configuration Settings for Anemone/Anemonefish SDM Project
# v2 - Includes paths for paper-like output structure
#-------------------------------------------------------------------------------
# Using pacman for streamlined package management
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, dplyr, terra, sf, stringr, ggplot2, readr, readxl,
               tidyr, vegan, ggvegan, corrplot, car, usdm, # For VIF/PCA
               robis, rgbif, # For occurrence download
               worrms, obistools, # For species metadata
               SDMtune, dismo, raster, # For SDM (raster maybe needed by SDMtune internals/plotting)
               devtools, # For biooracler if needed
               tools, # For file path manipulation
               ggtext, # Potentially useful for plotting
               parallel, lubridate, purrr, Hmisc, # Added Hmisc
               # --- Added for Parallel, Logging, Progress ---
               future, furrr, progressr, log4r
)

# --- Project Structure ---
base_dir <- here::here()

# --- Input Data Paths ---
data_dir         <- file.path(base_dir, "data")
occurrence_dir   <- file.path(data_dir, "occurrence")
env_data_dir     <- file.path(data_dir, "env") # Base env directory
shapefile_dir    <- file.path(data_dir, "shapefiles")

anemone_occurrence_dir <- file.path(occurrence_dir, "anemone")
anemonefish_occurrence_dir <- file.path(occurrence_dir, "anemonefish")

anemone_list_xlsx <- file.path(base_dir, "Anemone List.xlsx")
anemonefish_list_xlsx <- file.path(base_dir, "Anemonefish List.xlsx")

anemone_species_list_file <- file.path(base_dir, "final_anemone_species_list.csv")
anemonefish_species_list_file <- file.path(base_dir, "final_anemonefish_species_list.csv")

anemone_fish_association_file <- file.path(base_dir, "anemone_anemonefish_associations.csv")
processed_association_file <- file.path(data_dir, "processed_anemonefish_host_associations.csv") # Path for processed version

coral_shapefile <- file.path(shapefile_dir, "WCMC008_CoralReef2018_Py_v4_1.shp")

# --- Standard Output Data Paths (Intermediate/Backup) ---
sdm_output_dir   <- file.path(data_dir, "sdm_output")
log_dir_base     <- file.path(sdm_output_dir, "logs_and_analysis") # Base log dir
predictions_dir  <- file.path(sdm_output_dir, "predictions") # YOUR intermediate preds
results_dir      <- file.path(sdm_output_dir, "results")     # YOUR intermediate results
models_dir       <- file.path(sdm_output_dir, "models")      # YOUR intermediate models
species_log_dir  <- file.path(log_dir_base, "species_logs")

# --- Paper Structure Output Paths (Target for final analysis) ---
paper_output_base_dir <- base_dir # Set to your project base

# Predictions (Mimicking 'data/output/predictions/')
paper_pred_base_dir        <- file.path(paper_output_base_dir, "data", "output", "predictions")
paper_pred_future_base_dir <- file.path(paper_pred_base_dir, "Future") # For Future/sspXXX folders

# CV Results (Mimicking 'predictions/')
paper_results_base_dir <- file.path(paper_output_base_dir, "predictions") # For env-only models (plants in paper)
paper_results_biotic_dir <- file.path(paper_results_base_dir, "antcontrasts") # Mimicking 'antcontrasts' for biotic/combined CV

# Variable Importance (Mimicking 'data/output/vi/' and 'data/output/vi_ants/')
paper_vi_base_dir <- file.path(paper_output_base_dir, "data", "output", "vi") # For env-only VI (plants in paper)
paper_vi_biotic_dir <- file.path(paper_output_base_dir, "data", "output", "vi_ants") # For biotic/combined VI (ants in paper)

# --- Script Paths ---
scripts_dir      <- file.path(base_dir, "scripts")
helpers_dir      <- file.path(scripts_dir, "helpers")

# --- Processing Parameters ---
force_rerun <- list(
  get_species = FALSE,
  download_occurrences = FALSE,
  download_env = FALSE,
  preprocess_env_occurrence = FALSE, # Controls VIF/PCA reprocessing
  run_standard_sdms = TRUE,         # Controls 06a, 06b
  run_biotic_sdms = TRUE           # Controls 06c, 06d
)

occurrence_crs  <- "EPSG:4326"
env_crs         <- "EPSG:4326"

# --- Variables for Environmental Data Download (Script 04) ---
# These define what `04_download_env_data.R` will attempt to download if run.
# Ensure these match the actual datasets used if you modify the download process.
original_datasets_mean <- c(
  "par_mean_baseline_2000_2020_depthsurf", "sws_baseline_2000_2019_depthsurf",
  "thetao_baseline_2000_2019_depthsurf", "thetao_baseline_2000_2019_depthmax",
  "so_baseline_2000_2019_depthmax", "no3_baseline_2000_2018_depthmax",
  "chl_baseline_2000_2018_depthmax", "phyc_baseline_2000_2020_depthmax",
  "o2_baseline_2000_2018_depthmax", "ph_baseline_2000_2018_depthmax"
)
original_range_datasets <- c(
  "thetao_baseline_2000_2019_depthmax", "o2_baseline_2000_2018_depthmax",
  "no3_baseline_2000_2018_depthmax"
)
original_future_scenarios <- c("ssp119", "ssp585") # Used by script 04
original_time_steps <- list(
  dec50 = c("2030-01-01", "2040-01-01"),
  dec100 = c("2080-01-01", "2090-01-01")
)
original_terrain_vars_for_download <- c(
  "bathymetry_mean", "slope", "terrain_ruggedness_index" # Note: rugosity is calculated later
)
# Terrain Variables (final names after processing in script 04)
terrain_variables_final <- c("bathymetry_mean", "slope", "rugosity", "distcoast")
# Variables copied from current to future (needed by env loading helper)
vars_to_copy_to_future <- c(
  "chl_baseline_2000_2018_depthmax_mean",
  "par_mean_baseline_2000_2020_depthsurf"
  # Add other variables here if they are missing from future predictions
)

# --- Scenarios & Time Steps for SDM Runs ---
# Define the scenarios to actually model and predict (used by scripts 05-07)
env_scenarios <- c(
  "current" # Generally needed for tuning, can be commented out for final prediction runs if not needed
  # "ssp119_2050", "ssp119_2100",
  # "ssp585_2050", "ssp585_2100"
)

# Maps internal scenario names to the *directory* containing the files.
scenario_folder_map <- list(
  current     = file.path(env_data_dir, "current"),
  ssp119_2050 = file.path(env_data_dir, "future", "ssp119"),
  ssp119_2100 = file.path(env_data_dir, "future", "ssp119"),
  ssp585_2050 = file.path(env_data_dir, "future", "ssp585"),
  ssp585_2100 = file.path(env_data_dir, "future", "ssp585")
)
terrain_folder <- file.path(env_data_dir, "terrain")

# --- Predictor Selection Switch ---
use_pca_predictors <- TRUE # <<< SET THIS TO TRUE (for paper method) or FALSE
n_pca_components <- 4      # Number of PCA components to retain

# --- Define final VIF lists IF use_pca_predictors is FALSE ---
# *** ONLY NEEDED IF use_pca_predictors = FALSE ***
# *** These should be the final non-collinear lists based on 05a/b/c runs ***
final_vars_vif_anemone <- c(
  "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_range",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  "chl_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  # Terrain Vars
  "bathymetry_mean", "distcoast", "rugosity"
)
final_vars_vif_anemonefish_env <- c( # From 05b analysis
  "thetao_baseline_depthmax_mean",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "chl_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
final_vars_vif_anemonefish_combined_env <- c( # From 05c analysis - ENV vars ONLY
  "thetao_baseline_depthmax_mean",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "chl_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
vif_threshold <- 5             # VIF threshold if used
correlation_threshold <- 0.7   # Correlation threshold if used

# --- PCA/VIF Results Paths (Intermediate) ---
pca_model_save_path <- file.path(log_dir_base, "pca_model_current.rds") # PCA model based on current env
pca_raster_paths_rds_path <- file.path(log_dir_base, "pca_raster_paths.rds") # List of paths to projected PCA rasters
pca_models_rds_path <- file.path(log_dir_base, "pca_models.rds") # If storing PCA model per scenario (optional)
selected_vars_rds_path <- file.path(log_dir_base, "selected_variables_for_sdm.rds") # If saving VIF selections

# --- SDM Settings (SDMtune) ---
sdm_method <- "Maxnet"
sdm_partitions <- "randomkfold"
sdm_n_folds <- 5
sdm_tune_grid <- list(fc = c("l", "lq", "lh", "lp", "lqp"), reg = seq(0.5, 4, 0.5))
sdm_evaluation_metric <- "auc" # Metric for hyperparameter tuning selection
sdm_predict_type <- "cloglog"  # Prediction output type (raw, logistic, cloglog)
background_points_n <- 100
thinning_method <- "none" # Set to "cell" for cell-based thinning, "none" to skip

# Study Area / Environmental Constraints
apply_coral_mask <- TRUE
apply_depth_filter <- TRUE
depth_min <- -50 # meters (negative for depth below surface)
depth_max <- 0   # meters
min_occurrences_sdm <- 15 # Minimum occurrences AFTER cleaning/thinning

# --- Variable Importance Settings ---
run_variable_importance <- TRUE
vi_permutations <- 10 # Number of permutations for varImp (paper used 50, 10 is faster for testing)

# --- Parallel Processing ---
use_parallel <- TRUE # Set to FALSE to run sequentially
num_cores <- parallel::detectCores() - 1
if (num_cores < 1) num_cores <- 1
if (!use_parallel) num_cores <- 1 # Force 1 core if parallel is disabled

# --- Logging Configuration ---
log_file_path    <- file.path(log_dir_base, paste0("sdm_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_level        <- "INFO" # Options: "DEBUG", "INFO", "WARN", "ERROR", "FATAL"
log_append       <- TRUE
log_to_console   <- TRUE
log_console_level<- "INFO" # Use string here too

# --- Create Output Directories (Both standard and paper structure) ---
dir.create(sdm_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir_base, showWarnings = FALSE, recursive = TRUE)
dir.create(predictions_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(species_log_dir, showWarnings = FALSE, recursive = TRUE)
# Paper structure dirs
dir.create(paper_pred_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paper_pred_future_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paper_results_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paper_results_biotic_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paper_vi_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paper_vi_biotic_dir, showWarnings = FALSE, recursive = TRUE)
# Create SSP subdirs for future predictions
for(ssp in original_future_scenarios) {
  dir.create(file.path(paper_pred_future_base_dir, ssp), showWarnings=FALSE, recursive=TRUE)
}


# --- Display Name Function (Keep as is or adapt) ---
core_var_display_names <- c(
  "par_mean_baseline_2000_2020_depthsurf"  = "Surface PAR",
  "sws_baseline_2000_2019_depthsurf"  = "Surface Wind Stress",
  "thetao_baseline_2000_2019_depthsurf"= "Surface Temp.",
  "thetao_baseline_2000_2019_depthmax"= "Bottom Temp. Mean", # Assuming depthmax = mean if only one depthmax file exists
  "thetao_baseline_depthmax_mean"= "Bottom Temp. Mean",
  "thetao_baseline_depthmax_range"= "Bottom Temp. Range",
  "thetao_baseline_depthmax_ltmin" = "Bottom Temp. LT Min",
  "thetao_baseline_depthmax_ltmax" = "Bottom Temp. LT Max",
  "so_baseline_2000_2019_depthmax" = "Bottom Salinity", # Assuming depthmax = mean
  "so_baseline_depthmax_mean"    = "Bottom Salinity",
  "no3_baseline_2000_2018_depthmax"= "Bottom Nitrate Mean", # Assuming depthmax = mean
  "no3_baseline_depthmax_mean"   = "Bottom Nitrate Mean",
  "no3_baseline_depthmax_range"  = "Bottom Nitrate Range",
  "no3_baseline_depthmax_ltmin"  = "Bottom Nitrate LT Min",
  "no3_baseline_depthmax_ltmax"  = "Bottom Nitrate LT Max",
  "chl_baseline_2000_2018_depthmax"= "Bottom Chlorophyll", # Assuming depthmax = mean
  "chl_baseline_depthmax_mean"   = "Bottom Chlorophyll",
  "phyc_baseline_2000_2020_depthmax"= "Bottom Phytoplankton", # Assuming depthmax = mean
  "phyc_baseline_depthmax_mean"  = "Bottom Phytoplankton",
  "o2_baseline_2000_2018_depthmax" = "Bottom Oxygen Mean", # Assuming depthmax = mean
  "o2_baseline_depthmax_mean"    = "Bottom Oxygen Mean",
  "o2_baseline_depthmax_range"   = "Bottom Oxygen Range",
  "o2_baseline_depthmax_ltmin"   = "Bottom Oxygen LT Min",
  "o2_baseline_depthmax_ltmax"   = "Bottom Oxygen LT Max",
  "ph_baseline_2000_2018_depthmax" = "Bottom pH", # Assuming depthmax = mean
  "ph_baseline_depthmax_mean"    = "Bottom pH",
  "bathymetry_mean"     = "Bathymetry",
  "distcoast"           = "Distance to Coast",
  "rugosity"            = "Rugosity",
  "slope"               = "Slope",
  # PCA and Host names
  "PC1" = "PC1", "PC2" = "PC2", "PC3" = "PC3", "PC4" = "PC4", # Max 4 for now
  "host_suitability_max" = "Max Host Suitability"
)

get_display_name <- function(technical_name, lookup = core_var_display_names) {
  # Direct match first
  if (technical_name %in% names(lookup)) {
    return(lookup[technical_name])
  }
  # Try removing scenario/time tags for environmental vars
  core_name_cleaned <- gsub("_ssp\\d{3}_depth(surf|max)_dec\\d{3,3}", "_depth\\1", technical_name)
  core_name_cleaned <- gsub("_baseline(_\\d{4}_\\d{4})?", "", core_name_cleaned)
  # Add a check to remove trailing stats if they weren't in the original core name
  core_name_cleaned_no_stat <- gsub("(_mean|_range|_ltmin|_ltmax)$", "", core_name_cleaned)
  
  if (core_name_cleaned %in% names(lookup)) {
    return(lookup[core_name_cleaned])
  } else if (core_name_cleaned_no_stat %in% names(lookup)) {
    return(lookup[core_name_cleaned_no_stat])
  } else {
    # Fallback to original name if no match found
    # warning("No display name mapping found for: ", technical_name, call. = FALSE)
    return(technical_name)
  }
}


# --- Bundle settings into a list named 'config' ---
config <- list(
  # Paths (Input)
  base_dir = base_dir, data_dir = data_dir, occurrence_dir = occurrence_dir,
  env_data_dir = env_data_dir, shapefile_dir = shapefile_dir,
  anemone_occurrence_dir = anemone_occurrence_dir, anemonefish_occurrence_dir = anemonefish_occurrence_dir,
  anemone_list_xlsx = anemone_list_xlsx, anemonefish_list_xlsx = anemonefish_list_xlsx,
  anemone_species_list_file = anemone_species_list_file, anemonefish_species_list_file = anemonefish_species_list_file,
  anemone_fish_association_file = anemone_fish_association_file, processed_association_file = processed_association_file,
  coral_shapefile = coral_shapefile,
  
  # Paths (Standard Output - Intermediate/Backup)
  sdm_output_dir = sdm_output_dir, log_dir_base = log_dir_base, predictions_dir = predictions_dir,
  results_dir = results_dir, models_dir = models_dir, species_log_dir = species_log_dir,
  
  # Paths (Paper Structure Output - Target)
  paper_output_base_dir = paper_output_base_dir,
  paper_pred_base_dir = paper_pred_base_dir, paper_pred_future_base_dir = paper_pred_future_base_dir,
  paper_results_base_dir = paper_results_base_dir, paper_results_biotic_dir = paper_results_biotic_dir,
  paper_vi_base_dir = paper_vi_base_dir, paper_vi_biotic_dir = paper_vi_biotic_dir,
  
  # Paths (Scripts)
  scripts_dir = scripts_dir, helpers_dir = helpers_dir,
  
  # Processing Parameters
  force_rerun = force_rerun, occurrence_crs = occurrence_crs, env_crs = env_crs,
  original_datasets_mean = original_datasets_mean, original_range_datasets = original_range_datasets,
  original_future_scenarios = original_future_scenarios, original_time_steps = original_time_steps,
  original_terrain_vars_for_download = original_terrain_vars_for_download,
  terrain_variables_final = terrain_variables_final, vars_to_copy_to_future = vars_to_copy_to_future,
  env_scenarios = env_scenarios, scenario_folder_map = scenario_folder_map, terrain_folder = terrain_folder,
  
  # Predictor Selection
  use_pca_predictors = use_pca_predictors, n_pca_components = n_pca_components,
  final_vars_vif_anemone = final_vars_vif_anemone,
  final_vars_vif_anemonefish_env = final_vars_vif_anemonefish_env,
  final_vars_vif_anemonefish_combined_env = final_vars_vif_anemonefish_combined_env,
  vif_threshold = vif_threshold, correlation_threshold = correlation_threshold,
  
  # Results Paths (Intermediate)
  pca_model_save_path = pca_model_save_path, pca_raster_paths_rds_path = pca_raster_paths_rds_path,
  pca_models_rds_path = pca_models_rds_path, selected_vars_rds_path = selected_vars_rds_path,
  
  # SDM Settings
  sdm_method = sdm_method, sdm_partitions = sdm_partitions, sdm_n_folds = sdm_n_folds,
  sdm_tune_grid = sdm_tune_grid, sdm_evaluation_metric = sdm_evaluation_metric,
  sdm_predict_type = sdm_predict_type,
  background_points_n = background_points_n, thinning_method = thinning_method,
  apply_coral_mask = apply_coral_mask, apply_depth_filter = apply_depth_filter,
  depth_min = depth_min, depth_max = depth_max, min_occurrences_sdm = min_occurrences_sdm,
  
  # Variable Importance
  run_variable_importance = run_variable_importance, vi_permutations = vi_permutations,
  
  # Parallel & Logging
  num_cores = num_cores, use_parallel = use_parallel,
  log_file_path = log_file_path, log_level = log_level, log_append = log_append,
  log_to_console = log_to_console, log_console_level = log_console_level,
  
  # Display Names
  get_display_name = get_display_name, core_var_display_names = core_var_display_names
)

cat("Configuration loaded and bundled into 'config' list.\n")
cat("Base directory:", config$base_dir, "\n")
cat("Logging to file:", config$log_file_path, "(Level:", config$log_level, "Append:", config$log_append, ")\n")
cat("Logging to console:", config$log_to_console, "(Level:", config$log_console_level,")\n")
cat("Parallel execution:", config$use_parallel, "with", config$num_cores, "cores.\n")
cat("Using predictors:", ifelse(config$use_pca_predictors, "PCA", "VIF"), "\n")

#-------------------------------------------------------------------------------