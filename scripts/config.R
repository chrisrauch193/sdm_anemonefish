# scripts/config.R
#-------------------------------------------------------------------------------
# Configuration Settings for Anemone/Anemonefish SDM Project (v5 - Correct Path Definitions)
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

anemone_fish_association_file <- file.path(data_dir, "processed_anemonefish_host_associations.csv")

coral_shapefile <- file.path(shapefile_dir, "WCMC008_CoralReef2018_Py_v4_1.shp")

# --- Cropping and Masking Settings ---
apply_indo_pacific_crop <- TRUE # Set to TRUE to crop all rasters to the defined bbox
indo_pacific_bbox <- c(xmin=30, xmax=180, ymin=-50, ymax=50) # Define Lon/Lat bounding box

mask_background_points_to_coral <- TRUE # Set to TRUE to sample BG points ONLY from coral areas

# --- Intermediate Output Paths (for temp files, RDS models/tuning) ---
# Define these *before* the final config list is created
sdm_output_dir_intermediate   <- file.path(data_dir, "sdm_output_intermediate") # Base for intermediates
log_dir_base     <- file.path(sdm_output_dir_intermediate, "logs_and_analysis")
species_log_dir  <- file.path(log_dir_base, "species_logs")
models_dir_intermediate       <- file.path(sdm_output_dir_intermediate, "models")
results_dir_intermediate      <- file.path(sdm_output_dir_intermediate, "results")

# --- Target Output Paths (Mirroring desertantplantdistributions for final analysis) ---
target_analysis_output_base <- file.path(base_dir, "data", "output")
target_predictions_current_dir <- file.path(target_analysis_output_base, "predictions")
target_predictions_future_base <- file.path(target_predictions_current_dir, "Future")
target_results_base <- file.path(base_dir, "predictions") # For CV CSVs
target_vi_base <- file.path(target_analysis_output_base) # For VI CSVs

# Create the base target directories now (safe, as paths are defined above)
dir.create(target_analysis_output_base, recursive = TRUE, showWarnings = FALSE)
dir.create(target_predictions_current_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(target_predictions_future_base, recursive = TRUE, showWarnings = FALSE)
dir.create(target_results_base, recursive = TRUE, showWarnings = FALSE)
# Base intermediate dirs (logs etc.) will be created by scripts or logger setup now
dir.create(log_dir_base, recursive = TRUE, showWarnings = FALSE)
dir.create(species_log_dir, recursive = TRUE, showWarnings = FALSE) # Ensure species log dir exists

# --- Script Paths ---
scripts_dir      <- file.path(base_dir, "scripts")
helpers_dir      <- file.path(scripts_dir, "helpers")

# --- Processing Parameters ---
force_rerun <- list(
  get_species = FALSE,
  download_occurrences = FALSE,
  download_env = FALSE,
  preprocess_env_occurrence = FALSE,
  run_standard_sdms = TRUE, # Set to TRUE to force rerun of SDMs
  run_biotic_sdms = TRUE,
  run_enmeval_sdms = TRUE
)

occurrence_crs  <- "EPSG:4326"
env_crs         <- "EPSG:4326"

# --- !!! Original Download Lists !!! ---
# (Keep these as they were)
original_datasets_mean <- c("par_mean_baseline_2000_2020_depthsurf", "sws_baseline_2000_2019_depthsurf", "thetao_baseline_2000_2019_depthsurf", "thetao_baseline_2000_2019_depthmax", "so_baseline_2000_2019_depthmax", "no3_baseline_2000_2018_depthmax", "chl_baseline_2000_2018_depthmax", "phyc_baseline_2000_2020_depthmax", "o2_baseline_2000_2018_depthmax", "ph_baseline_2000_2018_depthmax")
original_range_datasets <- c("thetao_baseline_2000_2019_depthmax", "o2_baseline_2000_2018_depthmax", "no3_baseline_2000_2018_depthmax")
original_future_scenarios <- c("ssp119", "ssp585")
original_time_steps <- list(dec50 = c("2030-01-01", "2040-01-01"), dec100 = c("2080-01-01", "2090-01-01"))
original_terrain_vars_for_download <- c("bathymetry_mean", "slope", "terrain_ruggedness_index")

terrain_variables_final <- c("bathymetry_mean", "slope", "rugosity", "distcoast")
vars_to_copy_to_future <- c("chl_baseline_2000_2018_depthmax_mean", "par_mean_baseline_2000_2020_depthsurf")

# Scenarios and Time Steps
env_scenarios <- c("current", "ssp119_2050", "ssp119_2100", "ssp585_2050", "ssp585_2100")
scenario_folder_map <- list(current = file.path(env_data_dir, "current"), ssp119_2050 = file.path(env_data_dir, "future", "ssp119"), ssp119_2100 = file.path(env_data_dir, "future", "ssp119"), ssp585_2050 = file.path(env_data_dir, "future", "ssp585"), ssp585_2100 = file.path(env_data_dir, "future", "ssp585"))
terrain_folder <- file.path(env_data_dir, "terrain")
ssp_scenario_map <- list(ssp119_2050 = "ssp119", ssp119_2100 = "ssp119", ssp585_2050 = "ssp585", ssp585_2100 = "ssp585")
model_output_subdir_map <- list(`_pca` = "", `_biotic_pc1` = "biotic_pc1", `_combined_pca` = "combined_pca")

# Predictor Selection Switch
use_pca_predictors <- TRUE
n_pca_components <- 4

# VIF Lists (ONLY used if use_pca_predictors = FALSE)
final_vars_vif_anemone <- c("chl_baseline_2000_2018_depthmax_mean", "no3_baseline_2000_2018_depthmax_ltmin", "o2_baseline_2000_2018_depthmax_range", "par_mean_baseline_2000_2020_depthsurf", "so_baseline_2000_2019_depthmax_mean", "bathymetry_mean", "distcoast", "rugosity")
final_vars_vif_anemonefish_env <- c("chl_baseline_2000_2018_depthmax_mean", "no3_baseline_2000_2018_depthmax_ltmin", "o2_baseline_2000_2018_depthmax_range", "par_mean_baseline_2000_2020_depthsurf", "so_baseline_2000_2019_depthmax_mean", "bathymetry_mean", "distcoast", "rugosity")
final_vars_vif_anemonefish_combined_env <- c("chl_baseline_2000_2018_depthmax_mean", "no3_baseline_2000_2018_depthmax_ltmin", "o2_baseline_2000_2018_depthmax_range", "par_mean_baseline_2000_2020_depthsurf", "so_baseline_2000_2019_depthmax_mean", "bathymetry_mean", "distcoast", "rugosity")

# Paths for Intermediate PCA/VIF Results (saved in log dir)
pca_model_save_path <- file.path(log_dir_base, "pca_model_current.rds")
pca_raster_paths_rds_path <- file.path(log_dir_base, "pca_raster_paths.rds")
pca_models_rds_path <- file.path(log_dir_base, "pca_models.rds")

# SDM Settings
sdm_method <- "Maxnet"; sdm_partitions <- "randomkfold"; sdm_n_folds <- 5
sdm_tune_grid <- list(reg = seq(0.5, 4, 0.5), fc = c("l", "lq", "lh", "lp", "lqp"))
sdm_evaluation_metric <- "auc"; pca_background_points_n <- 100000; background_points_n <- 10000; thinning_method <- "cell"
apply_coral_mask <- TRUE; apply_depth_filter <- TRUE; depth_min <- -50; depth_max <- 0; min_occurrences_sdm <- 15

# --- Spatial Cross-Validation Settings (Simplified blockCV) ---
# ("spatial_grid" or "spatial_lat" or "random")
sdm_spatial_cv_type_to_use <- "spatial_grid"
blockcv_auto_range <- TRUE
blockcv_range_m <- 300000
blockcv_hexagon <- FALSE
# ("systematic", "random")
blockcv_selection <- "random"
blockcv_n_iterate <- 300
blockcv_lat_blocks <- 20000

# --- ENMeval Specific Settings ---
enmeval_algorithms <- c("maxnet") # Algorithm(s) to run (e.g., c("maxnet", "bioclim"))
enmeval_partitions <- "checkerboard1" # Partitioning method ("block", "checkerboard1", "checkerboard2", "jackknife", "randomkfold", "none", "testing")
enmeval_partition_settings <- list( # Settings specific to the chosen partition method
  # For checkerboard1/2: Define aggregation factor(s)
  aggregation.factor = 25 # For checkerboard1 (one value) or c(val1, val2) for checkerboard2
  # For block: Define orientation
  # orientation = "lat_lon" # Options: "lat_lon", "lon_lat", "lat_lat", "lon_lon"
  # For randomkfold: Define k
  # kfolds = 5
)
enmeval_tuning_settings <- list( # Define the tuning parameters grid
  fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), # Maxnet feature classes
  rm = seq(0.5, 4, 0.5)                             # Maxnet regularization multipliers
  # Add grids for other algorithms if enmeval_algorithms includes them
  # e.g., bioclim.threshold = c(0, 0.1, 0.3) # for bioclim
)
enmeval_num_cores <- 1 # Set > 1 for internal ENMeval parallel, but often better handled by future/furrr
enmeval_pred_type <- "cloglog" # Prediction output type ("cloglog", "raw", "logistic")
enmeval_clamp <- TRUE # Apply clamping during prediction?
enmeval_selection_metric <- "AICc" # Metric for selecting best model ('AICc', 'or.mtp.avg', 'auc.val.avg', 'cbi.val.avg', or define a custom sequence)
enmeval_delta_aicc_threshold <- 2 # Threshold for considering models equivalent based on delta AICc
enmeval_omit_na_models_from_selection <- TRUE # Remove models with NA evaluation stats before selection?

# Parallel & Logging
use_parallel <- TRUE; num_cores <- parallel::detectCores() - 1; if (num_cores < 1) num_cores <- 1; if (!use_parallel) num_cores <- 1
log_file_path <- file.path(log_dir_base, paste0("sdm_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_level <- "INFO"; log_append <- TRUE; log_to_console <- TRUE; log_console_level <- "INFO"

# Display Names
core_var_display_names <- c(par_baseline_depthsurf_mean="Surface PAR", sws_baseline_depthsurf_mean="Surface Wind Stress", thetao_baseline_depthsurf_mean="Surface Temp.", thetao_baseline_depthmax_mean="Bottom Temp. Mean", thetao_baseline_depthmax_range="Bottom Temp. Range", thetao_baseline_depthmax_ltmin="Bottom Temp. LT Min", thetao_baseline_depthmax_ltmax="Bottom Temp. LT Max", so_baseline_depthmax_mean="Bottom Salinity", no3_baseline_depthmax_mean="Bottom Nitrate Mean", no3_baseline_depthmax_range="Bottom Nitrate Range", no3_baseline_depthmax_ltmin="Bottom Nitrate LT Min", no3_baseline_depthmax_ltmax="Bottom Nitrate LT Max", chl_baseline_depthmax_mean="Bottom Chlorophyll", phyc_baseline_depthmax_mean="Bottom Phytoplankton", o2_baseline_depthmax_mean="Bottom Oxygen Mean", o2_baseline_depthmax_range="Bottom Oxygen Range", o2_baseline_depthmax_ltmin="Bottom Oxygen LT Min", o2_baseline_depthmax_ltmax="Bottom Oxygen LT Max", ph_baseline_depthmax_mean="Bottom pH", bathymetry_mean="Bathymetry", distcoast="Distance to Coast", rugosity="Rugosity", slope="Slope", PC1="PC1", PC2="PC2", PC3="PC3", PC4="PC4", host_suitability_max="Max Host Suitability")
get_display_name <- function(technical_name, lookup = NULL) { if (is.null(lookup)) lookup <- config$core_var_display_names; if (technical_name %in% names(lookup)) return(lookup[technical_name]); core_name_cleaned <- gsub("_ssp\\d{3}_depth(surf|max)_dec\\d{3,3}", "_depth\\1", technical_name); core_name_cleaned <- gsub("_baseline(_\\d{4}_\\d{4})?", "", core_name_cleaned); if (core_name_cleaned %in% names(lookup)) return(lookup[core_name_cleaned]); core_name_alt <- gsub("(_mean|_range|_ltmin|_ltmax)$", "", core_name_cleaned); if (core_name_alt %in% names(lookup)) return(lookup[core_name_alt]); return(technical_name) }


# --- Bundle settings into a list named 'config' ---
# *** Make sure intermediate paths are included here ***
config <- list(
  # Base Paths
  base_dir = base_dir, data_dir = data_dir, scripts_dir = scripts_dir, helpers_dir = helpers_dir,
  # Input Paths
  occurrence_dir = occurrence_dir, env_data_dir = env_data_dir, shapefile_dir = shapefile_dir,
  anemone_occurrence_dir = anemone_occurrence_dir, anemonefish_occurrence_dir = anemonefish_occurrence_dir,
  anemone_list_xlsx = anemone_list_xlsx, anemonefish_list_xlsx = anemonefish_list_xlsx,
  anemone_species_list_file = anemone_species_list_file, anemonefish_species_list_file = anemonefish_species_list_file,
  anemone_fish_association_file = anemone_fish_association_file, coral_shapefile = coral_shapefile,
  terrain_folder = terrain_folder,
  # Intermediate Output Paths <- Added Explicitly
  sdm_output_dir_intermediate = sdm_output_dir_intermediate,
  log_dir_base = log_dir_base,
  species_log_dir = species_log_dir,
  models_dir_intermediate = models_dir_intermediate, # <<< ADDED
  results_dir_intermediate = results_dir_intermediate, # <<< ADDED
  # Target Output Paths (for final analysis)
  target_analysis_output_base = target_analysis_output_base,
  target_predictions_current_dir = target_predictions_current_dir,
  target_predictions_future_base = target_predictions_future_base,
  target_results_base = target_results_base,
  target_vi_base = target_vi_base,
  # Scenario/Model Mapping
  env_scenarios = env_scenarios, scenario_folder_map = scenario_folder_map, ssp_scenario_map = ssp_scenario_map,
  model_output_subdir_map = model_output_subdir_map,
  # Processing Parameters
  force_rerun = force_rerun, occurrence_crs = occurrence_crs, env_crs = env_crs,
  original_datasets_mean = original_datasets_mean, original_range_datasets = original_range_datasets,
  original_future_scenarios = original_future_scenarios, original_time_steps = original_time_steps,
  original_terrain_vars_for_download = original_terrain_vars_for_download,
  terrain_variables_final = terrain_variables_final, vars_to_copy_to_future = vars_to_copy_to_future,
  use_pca_predictors = use_pca_predictors, n_pca_components = n_pca_components,
  final_vars_vif_anemone = final_vars_vif_anemone,
  final_vars_vif_anemonefish_env = final_vars_vif_anemonefish_env,
  final_vars_vif_anemonefish_combined_env = final_vars_vif_anemonefish_combined_env,
  # PCA/VIF Intermediate Results Paths
  pca_model_save_path = pca_model_save_path, pca_raster_paths_rds_path = pca_raster_paths_rds_path, pca_models_rds_path = pca_models_rds_path,
  # SDM Settings
  sdm_method = sdm_method, sdm_partitions = sdm_partitions, sdm_n_folds = sdm_n_folds,
  sdm_tune_grid = sdm_tune_grid, sdm_evaluation_metric = sdm_evaluation_metric,
  pca_background_points_n = pca_background_points_n, background_points_n = background_points_n, thinning_method = thinning_method,
  apply_coral_mask = apply_coral_mask, apply_indo_pacific_crop = apply_indo_pacific_crop,
  indo_pacific_bbox = indo_pacific_bbox, mask_background_points_to_coral = mask_background_points_to_coral,
  apply_depth_filter = apply_depth_filter, depth_min = depth_min, depth_max = depth_max,
  min_occurrences_sdm = min_occurrences_sdm,
  
  # --- Spatial Cross-Validation Settings (Simplified blockCV) ---
  sdm_spatial_cv_type_to_use = sdm_spatial_cv_type_to_use, # Which generated block type to use? ("spatial_grid" or "spatial_lat" or "random")
  blockcv_auto_range = blockcv_auto_range,       # TRUE: Calculate range based on autocorrelation
  blockcv_range_m = blockcv_range_m,        # Fixed range in METERS (used only if blockcv_auto_range = FALSE)
  blockcv_hexagon = blockcv_hexagon,          # Use hexagonal blocks for spatial_grid?
  blockcv_selection = blockcv_selection, # Fold assignment method ("systematic" or "random")
  blockcv_n_iterate = blockcv_n_iterate,          # Iterations for blockCV fold assignment (more relevant for 'random' selection)
  blockcv_lat_blocks = blockcv_lat_blocks,       # Only needed if using "spatial_lat"
  
  # --- ENMeval Settings ---
  enmeval_algorithms = enmeval_algorithms, enmeval_partitions = enmeval_partitions,
  enmeval_partition_settings = enmeval_partition_settings,
  enmeval_tuning_settings = enmeval_tuning_settings, enmeval_num_cores = enmeval_num_cores,
  enmeval_pred_type = enmeval_pred_type, enmeval_clamp = enmeval_clamp,
  enmeval_selection_metric = enmeval_selection_metric,
  enmeval_delta_aicc_threshold = enmeval_delta_aicc_threshold,
  enmeval_omit_na_models_from_selection = enmeval_omit_na_models_from_selection,
  
  # Parallel & Logging
  num_cores = num_cores, use_parallel = use_parallel,
  log_file_path = log_file_path, log_level = log_level, log_append = log_append,
  log_to_console = log_to_console, log_console_level = log_console_level,
  # Display Names
  get_display_name = get_display_name, core_var_display_names = core_var_display_names
  

  
)


# --- Final Check and Print Key Paths ---
cat("Configuration loaded and bundled into 'config' list.\n")
cat("Base directory:", config$base_dir, "\n")
cat("ENMeval scripts directory:", config$enmeval_scripts_dir, "\n") #<-- Added
cat("Intermediate SDM Output Dir:", config$sdm_output_dir_intermediate, "\n")
cat("Intermediate Models Dir:", config$models_dir_intermediate, "\n")
cat("Intermediate Results Dir:", config$results_dir_intermediate, "\n")
cat("Target Prediction Base:", config$target_predictions_current_dir, "\n")
cat("Target Results Base:", config$target_results_base, "\n")
cat("Logging to file:", config$log_file_path, "(Level:", config$log_level, "Append:", config$log_append, ")\n")
cat("Logging to console:", config$log_to_console, "(Level:", config$log_console_level,")\n")
cat("Parallel execution:", config$use_parallel, "with", config$num_cores, "cores.\n")
cat("Using predictors:", ifelse(config$use_pca_predictors, "PCA", "VIF"), "\n")
cat("ENMeval algorithm(s):", paste(config$enmeval_algorithms, collapse=", "), "\n")
cat("ENMeval partition method:", config$enmeval_partitions, "\n")
#-------------------------------------------------------------------------------