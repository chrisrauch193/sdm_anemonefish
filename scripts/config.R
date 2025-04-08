# scripts/config.R
#-------------------------------------------------------------------------------
# Configuration Settings for Anemone/Anemonefish SDM Project
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

coral_shapefile <- file.path(shapefile_dir, "WCMC008_CoralReef2018_Py_v4_1.shp")

# --- Output Data Paths ---
sdm_output_dir   <- file.path(data_dir, "sdm_output")
log_dir_base     <- file.path(sdm_output_dir, "logs_and_analysis") # Base log dir
predictions_dir  <- file.path(sdm_output_dir, "predictions")
results_dir      <- file.path(sdm_output_dir, "results")
models_dir       <- file.path(sdm_output_dir, "models")
species_log_dir  <- file.path(log_dir_base, "species_logs") # <--- ADD THIS

# --- Script Paths ---
scripts_dir      <- file.path(base_dir, "scripts")
helpers_dir      <- file.path(scripts_dir, "helpers")

# --- Processing Parameters ---
force_rerun <- list(
  get_species = FALSE,
  download_occurrences = FALSE,
  download_env = FALSE,
  preprocess_env_occurrence = FALSE,
  run_standard_sdms = FALSE, # Set to TRUE to force rerun of SDMs
  run_biotic_sdms = FALSE
)

occurrence_crs  <- "EPSG:4326"
env_crs         <- "EPSG:4326"

# --- !!! Original Download Lists !!! ---
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
original_future_scenarios <- c("ssp119", "ssp585")
original_time_steps <- list(
  dec50 = c("2030-01-01", "2040-01-01"), # Corresponds to "_dec50_" in filename
  dec100 = c("2080-01-01", "2090-01-01") # Corresponds to "_dec100_" in filename
)
original_terrain_vars_for_download <- c(
  "bathymetry_mean", "slope", "terrain_ruggedness_index"
)
# --- End Original Download Lists ---

terrain_variables_final <- c("bathymetry_mean", "slope", "rugosity", "distcoast")

# --- !!! Variables Missing in Future Data !!! ---
vars_to_copy_to_future <- c(
  "chl_baseline_2000_2018_depthmax_mean", # Assuming this is the final CHL name
  "par_mean_baseline_2000_2020_depthsurf"  # Assuming this is the final PAR name
)

# Scenarios and Time Steps (Internal logical names used in loops)
env_scenarios <- c(
  "current"
  # Add your future scenarios back here when needed:
  # "ssp119_2050", "ssp119_2100",
  # "ssp585_2050", "ssp585_2100"
)

# Maps internal scenario names to the *directory* containing the files.
scenario_folder_map <- list(
  current     = file.path(env_data_dir, "current"),
  ssp119_2050 = file.path(env_data_dir, "future", "ssp119"), # Folder for ALL ssp119 files
  ssp119_2100 = file.path(env_data_dir, "future", "ssp119"),
  ssp585_2050 = file.path(env_data_dir, "future", "ssp585"),
  ssp585_2100 = file.path(env_data_dir, "future", "ssp585")
)
terrain_folder <- file.path(env_data_dir, "terrain")

# --- Predictor Selection Switch ---
# Set to TRUE to use PCA predictors (Braun et al. method, script 05b)
# Set to FALSE to use VIF-selected raw predictors (defined below, requires 05a/c run)
use_pca_predictors <- TRUE # <<< SET THIS TO TRUE (default) or FALSE

# --- Define final VIF lists IF use_pca_predictors is FALSE ---
# These should be the *final* lists determined after running 05a/05b/05c
# *** ONLY NEEDED IF use_pca_predictors = FALSE ***
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
  "bathymetry_mean", "distcoast", "rugosity"
)
final_vars_vif_anemonefish_env <- c( # From 05b
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
final_vars_vif_anemonefish_combined_env <- c( # From 05c - ENV vars ONLY
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
# Note: The host predictor is added dynamically in script 07

# --- PCA/VIF Results Paths ---
# Path for the PCA model object (needed for projecting onto future data in 05b)
pca_model_save_path <- file.path(log_dir_base, "pca_model_current.rds")
# Path for the list of saved PCA raster stack paths (output of 05b)
pca_raster_paths_rds_path <- file.path(log_dir_base, "pca_raster_paths.rds")

# SDM Settings (SDMtune)
sdm_method <- "Maxnet"
sdm_partitions <- "randomkfold"
sdm_n_folds <- 5
sdm_tune_grid <- list(fc = c("l", "lq", "lh", "lp", "lqp"), reg = seq(0.5, 4, 0.5))
sdm_evaluation_metric <- "auc" # Metric for hyperparameter tuning selection
background_points_n <- 10000
thinning_method <- "none" # Set to "cell" for cell-based thinning, "none" to skip

# Study Area / Environmental Constraints
apply_coral_mask <- TRUE
apply_depth_filter <- TRUE
depth_min <- -50
depth_max <- 0
min_occurrences_sdm <- 15 # Minimum occurrences AFTER cleaning/thinning to run SDM

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

# --- Create Output Directories ---
dir.create(sdm_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir_base, showWarnings = FALSE, recursive = TRUE)
dir.create(predictions_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(env_data_dir, showWarnings = FALSE, recursive = TRUE)

# Function to get display name (kept from previous version)
core_var_display_names <- c(
  "par_baseline_depthsurf_mean"  = "Surface PAR",
  "sws_baseline_depthsurf_mean"  = "Surface Wind Stress",
  "thetao_baseline_depthsurf_mean"= "Surface Temp.",
  "thetao_baseline_depthmax_mean"= "Bottom Temp. Mean",
  "thetao_baseline_depthmax_range"= "Bottom Temp. Range",
  "thetao_baseline_depthmax_ltmin" = "Bottom Temp. LT Min",
  "thetao_baseline_depthmax_ltmax" = "Bottom Temp. LT Max",
  "so_baseline_depthmax_mean"    = "Bottom Salinity",
  "no3_baseline_depthmax_mean"   = "Bottom Nitrate Mean",
  "no3_baseline_depthmax_range"  = "Bottom Nitrate Range",
  "no3_baseline_depthmax_ltmin"  = "Bottom Nitrate LT Min",
  "no3_baseline_depthmax_ltmax"  = "Bottom Nitrate LT Max",
  "chl_baseline_depthmax_mean"   = "Bottom Chlorophyll",
  "phyc_baseline_depthmax_mean"  = "Bottom Phytoplankton",
  "o2_baseline_depthmax_mean"    = "Bottom Oxygen Mean",
  "o2_baseline_depthmax_range"   = "Bottom Oxygen Range",
  "o2_baseline_depthmax_ltmin"   = "Bottom Oxygen LT Min",
  "o2_baseline_depthmax_ltmax"   = "Bottom Oxygen LT Max",
  "ph_baseline_depthmax_mean"    = "Bottom pH",
  "bathymetry_mean"     = "Bathymetry",
  "distcoast"           = "Distance to Coast",
  "rugosity"            = "Rugosity",
  "slope"               = "Slope"
)

get_display_name <- function(technical_name, lookup = core_var_display_names) {
  core_name <- technical_name # Default
  if (technical_name %in% names(lookup)) {
    core_name <- technical_name
  } else {
    core_name_cleaned <- gsub("_ssp\\d{3}_depth(surf|max)_dec\\d{3,3}", "_depth\\1", technical_name)
    core_name_cleaned <- gsub("_baseline(_\\d{4}_\\d{4})?", "", core_name_cleaned)
    if (core_name_cleaned %in% names(lookup)) {
      core_name <- core_name_cleaned
    } else {
      core_name_alt <- gsub("(_mean|_range|_ltmin|_ltmax)$", "", core_name_cleaned)
      if (core_name_alt %in% names(lookup)) {
        core_name <- core_name_alt
      } else {
        # Warning is now generated within the plotting functions if needed
        # warning("No display name mapping found for: ", technical_name, call. = FALSE)
        return(technical_name) # Return original if no match
      }
    }
  }
  display_name <- lookup[core_name]
  return(ifelse(is.na(display_name), technical_name, display_name))
}

# --- Bundle settings into a list named 'config' ---
config <- list(
  # Paths
  base_dir = base_dir, data_dir = data_dir, occurrence_dir = occurrence_dir,
  env_data_dir = env_data_dir, shapefile_dir = shapefile_dir,
  anemone_occurrence_dir = anemone_occurrence_dir, anemonefish_occurrence_dir = anemonefish_occurrence_dir,
  anemone_list_xlsx = anemone_list_xlsx, anemonefish_list_xlsx = anemonefish_list_xlsx,
  anemone_species_list_file = anemone_species_list_file, anemonefish_species_list_file = anemonefish_species_list_file,
  anemone_fish_association_file = anemone_fish_association_file, coral_shapefile = coral_shapefile,
  sdm_output_dir = sdm_output_dir, log_dir_base = log_dir_base, predictions_dir = predictions_dir,
  species_log_dir = species_log_dir,
  results_dir = results_dir, models_dir = models_dir, scripts_dir = scripts_dir,
  helpers_dir = helpers_dir, terrain_folder = terrain_folder,
  pca_model_save_path = pca_model_save_path, pca_raster_paths_rds_path = pca_raster_paths_rds_path,
  
  # Processing Parameters
  force_rerun = force_rerun, occurrence_crs = occurrence_crs, env_crs = env_crs,
  original_datasets_mean = original_datasets_mean, original_range_datasets = original_range_datasets,
  original_future_scenarios = original_future_scenarios, original_time_steps = original_time_steps,
  original_terrain_vars_for_download = original_terrain_vars_for_download,
  terrain_variables_final = terrain_variables_final, vars_to_copy_to_future = vars_to_copy_to_future,
  env_scenarios = env_scenarios, scenario_folder_map = scenario_folder_map,
  use_pca_predictors = use_pca_predictors,
  # VIF Lists (only used if use_pca_predictors = FALSE)
  final_vars_vif_anemone = final_vars_vif_anemone,
  final_vars_vif_anemonefish_env = final_vars_vif_anemonefish_env,
  final_vars_vif_anemonefish_combined_env = final_vars_vif_anemonefish_combined_env,
  
  # SDM Settings
  sdm_method = sdm_method, sdm_partitions = sdm_partitions, sdm_n_folds = sdm_n_folds,
  sdm_tune_grid = sdm_tune_grid, sdm_evaluation_metric = sdm_evaluation_metric,
  background_points_n = background_points_n, thinning_method = thinning_method,
  apply_coral_mask = apply_coral_mask, apply_depth_filter = apply_depth_filter,
  depth_min = depth_min, depth_max = depth_max, min_occurrences_sdm = min_occurrences_sdm,
  
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