# scripts/config.R
#-------------------------------------------------------------------------------
# Configuration Settings for Anemone/Anemonefish SDM Project (v5 - Correct Path Definitions)
#-------------------------------------------------------------------------------

# setwd("~/a0236995/sdm_anemonefish")

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

# Set working directory
wd = "/home/bi-server-kyoto/a0236995/sdm_anemonefish"
setwd(wd)

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
indo_pacific_bbox <- c(xmin=-180, xmax=180, ymin=-40, ymax=45) # Define Lon/Lat bounding box

mask_background_points_to_coral <- FALSE # Set to TRUE to sample BG points ONLY from coral areas

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
  run_biotic_sdms = TRUE
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
model_output_subdir_map <- list(`_pca` = "_pca", `_biotic_only` = "_biotic_only", `_combined_pca` = "_combined_pca")

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
apply_coral_mask <- FALSE; depth_min <- -200; pca_temp_min_threshold <- 20; depth_max <- 0; min_occurrences_sdm <- 15

# Parallel & Logging
use_parallel <- TRUE; num_cores <- parallel::detectCores() - 1; if (num_cores < 1) num_cores <- 1; if (!use_parallel) num_cores <- 1
log_file_path <- file.path(log_dir_base, paste0("sdm_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_level <- "INFO"; log_append <- TRUE; log_to_console <- TRUE; log_console_level <- "INFO"

# Display Names
core_var_display_names <- c(par_baseline_depthsurf_mean="Surface PAR", sws_baseline_depthsurf_mean="Surface Wind Stress", thetao_baseline_depthsurf_mean="Surface Temp.", thetao_baseline_depthmax_mean="Bottom Temp. Mean", thetao_baseline_depthmax_range="Bottom Temp. Range", thetao_baseline_depthmax_ltmin="Bottom Temp. LT Min", thetao_baseline_depthmax_ltmax="Bottom Temp. LT Max", so_baseline_depthmax_mean="Bottom Salinity", no3_baseline_depthmax_mean="Bottom Nitrate Mean", no3_baseline_depthmax_range="Bottom Nitrate Range", no3_baseline_depthmax_ltmin="Bottom Nitrate LT Min", no3_baseline_depthmax_ltmax="Bottom Nitrate LT Max", chl_baseline_depthmax_mean="Bottom Chlorophyll", phyc_baseline_depthmax_mean="Bottom Phytoplankton", o2_baseline_depthmax_mean="Bottom Oxygen Mean", o2_baseline_depthmax_range="Bottom Oxygen Range", o2_baseline_depthmax_ltmin="Bottom Oxygen LT Min", o2_baseline_depthmax_ltmax="Bottom Oxygen LT Max", ph_baseline_depthmax_mean="Bottom pH", bathymetry_mean="Bathymetry", distcoast="Distance to Coast", rugosity="Rugosity", slope="Slope", PC1="PC1", PC2="PC2", PC3="PC3", PC4="PC4", host_suitability_max="Max Host Suitability")
get_display_name <- function(technical_name, lookup = NULL) { if (is.null(lookup)) lookup <- config$core_var_display_names; if (technical_name %in% names(lookup)) return(lookup[technical_name]); core_name_cleaned <- gsub("_ssp\\d{3}_depth(surf|max)_dec\\d{3,3}", "_depth\\1", technical_name); core_name_cleaned <- gsub("_baseline(_\\d{4}_\\d{4})?", "", core_name_cleaned); if (core_name_cleaned %in% names(lookup)) return(lookup[core_name_cleaned]); core_name_alt <- gsub("(_mean|_range|_ltmin|_ltmax)$", "", core_name_cleaned); if (core_name_alt %in% names(lookup)) return(lookup[core_name_alt]); return(technical_name) }

# --- Spatial Cross-Validation Settings (Simplified blockCV) ---
# ("spatial_grid" or "spatial_lat" or "random")
sdm_spatial_cv_type_to_use <- "spatial_grid"
blockcv_auto_range <- TRUE
blockcv_range_default <- 300000 # 20000 300000
blockcv_range_max <- 1000000 # 1000000
blockcv_hexagon <- TRUE
# ("systematic", "random")
blockcv_selection <- "systematic"
blockcv_n_iterate <- 300
blockcv_lat_blocks <- 10

# OBIS stuff
ecoregion_shapefile <- file.path(shapefile_dir, "MarineRealms_BO.shp")
bathymetry_file     <- file.path(terrain_folder, "bathymetry_mean.tif")
limit_by_depth_obis <- TRUE  # Or FALSE
poly_buffer_obis    <- 0.2 # Small degree buffer for adjacency
poly_buffer_final   <- 0.5 # Larger degree buffer for final extent (optional, replicating OBIS)


# Spatial Autocorrelation (SAC) Thinning Settings
apply_sac_thinning <- TRUE     # Apply Mantel test + spThin thinning
autocor_classdist <- 50000      # 50 km steps for Mantel correlogram
autocor_maxdist <- 1000000   # 1000 km max distance for correlogram
autocor_signif <- 0.05      # Significance level for non-correlation
sac_prune_threshold <- 20000   # Thin if non-sig distance >= 20 km (based on BlockCV results)

do_final_prediction <- TRUE

global_seed = 1

meow_provinces_shapefile <- file.path(shapefile_dir, "meow_ecos.shp")

depth_raster_path <- file.path(env_data_dir, "terrain/bathymetry_mean.tif")

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
  depth_min = depth_min, pca_temp_min_threshold = pca_temp_min_threshold, depth_max = depth_max,
  min_occurrences_sdm = min_occurrences_sdm,
  # Parallel & Logging
  num_cores = num_cores, use_parallel = use_parallel,
  log_file_path = log_file_path, log_level = log_level, log_append = log_append,
  log_to_console = log_to_console, log_console_level = log_console_level,
  # Display Names
  get_display_name = get_display_name, core_var_display_names = core_var_display_names,
  
  # --- Spatial Cross-Validation Settings (Simplified blockCV) ---
  sdm_spatial_cv_type_to_use = sdm_spatial_cv_type_to_use, # Which generated block type to use? ("spatial_grid" or "spatial_lat" or "random")
  blockcv_auto_range = blockcv_auto_range,       # TRUE: Calculate range based on autocorrelation
  blockcv_range_default = blockcv_range_default,        # Fixed range in METERS (used only if blockcv_auto_range = FALSE)
  blockcv_range_max = blockcv_range_max,
  blockcv_hexagon = blockcv_hexagon,          # Use hexagonal blocks for spatial_grid?
  blockcv_selection = blockcv_selection, # Fold assignment method ("systematic" or "random")
  blockcv_n_iterate = blockcv_n_iterate,          # Iterations for blockCV fold assignment (more relevant for 'random' selection)
  blockcv_lat_blocks = blockcv_lat_blocks,       # Only needed if using "spatial_lat"
  
  # OBIS Stuff
  ecoregion_shapefile = ecoregion_shapefile,
  bathymetry_file = bathymetry_file,
  limit_by_depth_obis = limit_by_depth_obis,
  poly_buffer_obis = poly_buffer_obis,
  poly_buffer_final = poly_buffer_final,
  
  # Spatial Autocorrelation (SAC) Thinning Settings
  apply_sac_thinning = apply_sac_thinning,
  autocor_classdist = autocor_classdist,
  autocor_maxdist = autocor_maxdist,
  autocor_signif = autocor_signif,
  sac_prune_threshold = sac_prune_threshold,
  
  do_final_prediction = do_final_prediction,
  
  global_seed = global_seed,
  
  meow_provinces_shapefile = meow_provinces_shapefile,
  
  depth_raster_path = depth_raster_path
)


# --- Final Check and Print Key Paths ---
cat("Configuration loaded and bundled into 'config' list.\n")
cat("Base directory:", config$base_dir, "\n")
cat("Intermediate SDM Output Dir:", config$sdm_output_dir_intermediate, "\n") # Added intermediate
cat("Intermediate Models Dir:", config$models_dir_intermediate, "\n")      # Added intermediate
cat("Intermediate Results Dir:", config$results_dir_intermediate, "\n")     # Added intermediate
cat("Target Prediction Base:", config$target_predictions_current_dir, "\n")
cat("Target Results Base:", config$target_results_base, "\n")
cat("Logging to file:", config$log_file_path, "(Level:", config$log_level, "Append:", config$log_append, ")\n")
cat("Logging to console:", config$log_to_console, "(Level:", config$log_console_level,")\n")
cat("Parallel execution:", config$use_parallel, "with", config$num_cores, "cores.\n")
cat("Using predictors:", ifelse(config$use_pca_predictors, "PCA", "VIF"), "\n")

#-------------------------------------------------------------------------------