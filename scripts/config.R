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
               ENMeval, dismo, raster, # For SDM (raster maybe needed by ENMeval internals/plotting)
               devtools, # For biooracler if needed
               tools, # For file path manipulation
               ggtext, # Potentially useful for plotting
               parallel, lubridate, purrr, Hmisc) # Added Hmisc

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
log_dir          <- file.path(sdm_output_dir, "logs_and_analysis")
predictions_dir  <- file.path(sdm_output_dir, "predictions")
results_dir      <- file.path(sdm_output_dir, "results")
models_dir       <- file.path(sdm_output_dir, "models")

# --- Script Paths ---
scripts_dir      <- file.path(base_dir, "scripts")
helpers_dir      <- file.path(scripts_dir, "helpers")

# --- Processing Parameters ---
force_rerun <- list(
  get_species = FALSE,
  download_occurrences = FALSE,
  download_env = FALSE,
  preprocess_env_occurrence = FALSE,
  run_standard_sdms = FALSE,
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
# List the EXACT FILENAMES (without extension) of variables in the 'current' folder
# that need to be copied to future scenario folders because they lack future versions.
vars_to_copy_to_future <- c(
  "chl_baseline_depthmax_mean", # Assuming this is the final CHL name
  "par_baseline_depthsurf_mean"  # Assuming this is the final PAR name
  # Add any others if necessary
)

# Scenarios and Time Steps (Internal logical names used in loops)
env_scenarios <- c(
  "current",
  "ssp119_2050", # Represents the dec50 time step
  "ssp119_2100", # Represents the dec100 time step
  "ssp585_2050", # Represents the dec50 time step
  "ssp585_2100"  # Represents the dec100 time step
)

# --- !!! UPDATED scenario_folder_map !!! ---
# Maps internal scenario names to the *directory* containing the files.
# For future, points to the SSP directory, NOT the time step subdir.
scenario_folder_map <- list(
  current     = file.path(env_data_dir, "current"),
  ssp119_2050 = file.path(env_data_dir, "future", "ssp119"), # Folder for ALL ssp119 files
  ssp119_2100 = file.path(env_data_dir, "future", "ssp119"), # Folder for ALL ssp119 files
  ssp585_2050 = file.path(env_data_dir, "future", "ssp585"), # Folder for ALL ssp585 files
  ssp585_2100 = file.path(env_data_dir, "future", "ssp585")  # Folder for ALL ssp585 files
)
# Define the terrain folder path explicitly
terrain_folder <- file.path(env_data_dir, "terrain")

# Define Chlorophyll details
chl_variable_stem <- "chl_baseline_2000_2018_depthmax"
chl_file_path <- file.path(scenario_folder_map[["current"]], paste0(chl_variable_stem, "_mean.tif"))

# VIF/Correlation Settings
vif_threshold <- 10
correlation_threshold <- 0.8
n_pca_components <- 4

# SDM Settings (Updated for sdmtune)
sdm_method <- "Maxnet" # sdmtune uses method names like "Maxnet", "RF", etc.
sdm_partitions <- "randomkfold" # sdmtune uses this term
sdm_n_folds <- 5
# sdmtune uses 'reg' for regularization multiplier and 'fc' for feature class
sdm_tune_grid <- list(fc = c("L", "LQ", "H", "LQH", "LQP"), reg = seq(0.5, 4, 0.5))
sdm_evaluation_metric <- "AUC" # Common metric used by sdmtune tuning/results
background_points_n <- 10000
thinning_method <- "cell" # Or your preferred method

# Study Area / Environmental Constraints
apply_coral_mask <- TRUE
apply_depth_filter <- TRUE
depth_min <- -50
depth_max <- 0
min_occurrences_sdm <- 10

# --- Placeholders for VIF/PCA results ---
selected_variables_for_pca <- list()
selected_vars_rds_path <- file.path(log_dir, "selected_variables_for_pca.rds")
pca_raster_paths_rds_path <- file.path(log_dir, "pca_raster_paths.rds")
pca_models_rds_path <- file.path(log_dir, "pca_models.rds")

# --- Parallel Processing ---
num_cores <- parallel::detectCores() - 1
if (num_cores < 1) num_cores <- 1
use_parallel <- TRUE

# --- Create Output Directories ---
dir.create(sdm_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(predictions_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(env_data_dir, showWarnings = FALSE, recursive = TRUE)

# Define CORE variable stems and their display names
core_var_display_names <- c(
  "par_depthsurf_mean"  = "Surface PAR",
  "sws_depthsurf_mean"  = "Surface Wind Stress",
  "thetao_depthsurf_mean"= "Surface Temp.", # Added Surf Temp
  "thetao_depthmax_mean"= "Bottom Temp. Mean",
  "thetao_depthmax_range"= "Bottom Temp. Range",
  "thetao_depthmax_ltmin" = "Bottom Temp. LT Min", # Added LTMin/Max
  "thetao_depthmax_ltmax" = "Bottom Temp. LT Max",
  "so_depthmax_mean"    = "Bottom Salinity",
  "no3_depthmax_mean"   = "Bottom Nitrate Mean",
  "no3_depthmax_range"  = "Bottom Nitrate Range",
  "no3_depthmax_ltmin"  = "Bottom Nitrate LT Min", # Added LTMin/Max
  "no3_depthmax_ltmax"  = "Bottom Nitrate LT Max",
  "chl_depthmax_mean"   = "Bottom Chlorophyll",
  "phyc_depthmax_mean"  = "Bottom Phytoplankton", # Added Phyc
  "o2_depthmax_mean"    = "Bottom Oxygen Mean", # Added O2 Mean
  "o2_depthmax_range"   = "Bottom Oxygen Range",
  "o2_depthmax_ltmin"   = "Bottom Oxygen LT Min", # Added LTMin/Max
  "o2_depthmax_ltmax"   = "Bottom Oxygen LT Max",
  "ph_depthmax_mean"    = "Bottom pH", # Added pH
  "bathymetry_mean"     = "Bathymetry",
  "distcoast"           = "Distance to Coast",
  "rugosity"            = "Rugosity",
  "slope"               = "Slope" # Added Slope
  # Add other CORE variable stems as needed
)

final_vars_anemone <- c(
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

# Function to get display name from a potentially scenario-specific technical name
get_display_name <- function(technical_name, lookup = core_var_display_names) {
  # Remove scenario/time tags to find the core name
  # Regex tries to remove _baseline_YYYY_YYYY, _sspXXX_decYY, etc.
  core_name <- gsub("_ssp\\d{3}_depth(surf|max)_dec\\d{3,3}_", "_depth\\1_", technical_name) # Future pattern
  core_name <- gsub("_baseline_\\d{4}_\\d{4}_", "_", core_name) # Baseline pattern
  core_name <- gsub("_baseline_", "_", core_name) # Simpler baseline if years missing
  
  # Special case for terrain vars that don't have scenario tags
  if (technical_name %in% names(lookup)) {
    core_name <- technical_name
  } else {
    # More robust: find which key *starts* the technical name
    possible_keys <- names(lookup)[sapply(names(lookup), function(key) startsWith(technical_name, key))]
    if(length(possible_keys) == 1) {
      core_name <- possible_keys[1]
    } else {
      # Fallback if pattern matching failed or multiple matches (less likely with good keys)
      core_name_alt <- gsub("_ssp\\d{3}_dec\\d{3,3}", "", technical_name) # Try simpler remove
      core_name_alt <- gsub("_baseline_\\d{4}_\\d{4}", "", core_name_alt)
      core_name_alt <- gsub("(_mean|_range|_ltmin|_ltmax)$", "", core_name_alt) # Remove stat suffix? Maybe too risky
      # Find the best match based on what's left if primary didn't work
      if (core_name_alt %in% names(lookup)) {
        core_name <- core_name_alt
      } else {
        # Last resort: return original if no match found after cleaning
        warning("No display name mapping found for: ", technical_name, call. = FALSE)
        return(technical_name)
      }
    }
  }
  
  
  display_name <- lookup[core_name]
  # Return original technical name if lookup failed
  return(ifelse(is.na(display_name), technical_name, display_name))
}
# labels = sapply(technical_names_vector, get_display_name)

# --- Bundle settings into a list named 'config' ---
config <- list(
  # Paths
  base_dir = base_dir, data_dir = data_dir, occurrence_dir = occurrence_dir,
  env_data_dir = env_data_dir, shapefile_dir = shapefile_dir,
  anemone_occurrence_dir = anemone_occurrence_dir, anemonefish_occurrence_dir = anemonefish_occurrence_dir,
  anemone_list_xlsx = anemone_list_xlsx, anemonefish_list_xlsx = anemonefish_list_xlsx,
  anemone_species_list_file = anemone_species_list_file, anemonefish_species_list_file = anemonefish_species_list_file,
  anemone_fish_association_file = anemone_fish_association_file, coral_shapefile = coral_shapefile,
  sdm_output_dir = sdm_output_dir, log_dir = log_dir, predictions_dir = predictions_dir,
  results_dir = results_dir, models_dir = models_dir, scripts_dir = scripts_dir,
  helpers_dir = helpers_dir, terrain_folder = terrain_folder, chl_file_path = chl_file_path,
  selected_vars_rds_path = selected_vars_rds_path, pca_raster_paths_rds_path = pca_raster_paths_rds_path,
  pca_models_rds_path = pca_models_rds_path,
  
  # Processing Parameters
  force_rerun = force_rerun, occurrence_crs = occurrence_crs, env_crs = env_crs,
  # --- Include Original Download Lists in Config ---
  original_datasets_mean = original_datasets_mean,
  original_range_datasets = original_range_datasets,
  original_future_scenarios = original_future_scenarios,
  original_time_steps = original_time_steps,
  original_terrain_vars_for_download = original_terrain_vars_for_download,
  # --- End Original Download Lists ---
  terrain_variables_final = terrain_variables_final,
  vars_to_copy_to_future = vars_to_copy_to_future,
  env_scenarios = env_scenarios, scenario_folder_map = scenario_folder_map, # Use updated map
  chl_variable_stem = chl_variable_stem,
  vif_threshold = vif_threshold, correlation_threshold = correlation_threshold,
  n_pca_components = n_pca_components, sdm_method = sdm_method,
  sdm_partitions = sdm_partitions, sdm_n_folds = sdm_n_folds,
  sdm_tune_grid = sdm_tune_grid, sdm_evaluation_metric = sdm_evaluation_metric,
  background_points_n = background_points_n, thinning_method = thinning_method,
  apply_coral_mask = apply_coral_mask, apply_depth_filter = apply_depth_filter,
  depth_min = depth_min, depth_max = depth_max, min_occurrences_sdm = min_occurrences_sdm,
  selected_variables_for_pca = selected_variables_for_pca, num_cores = num_cores,
  use_parallel = use_parallel,
  get_display_name = get_display_name,
  final_vars_anemone = final_vars_anemone
)

cat("Configuration loaded and bundled into 'config' list.\n")
cat("Base directory:", config$base_dir, "\n")
cat("Using", config$num_cores, "cores for parallel tasks (if enabled).\n")

#-------------------------------------------------------------------------------