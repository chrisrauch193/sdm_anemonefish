# scripts/vif/anemone_vif.R
# Driver script for running VIF analysis on Anemone data.

library(tools) # Used for file_path_as_absolute

# --- Scenario Specific Parameters ---

# Define paths relative to the project root or use absolute paths
project_root <- getwd() # Or set explicitly: "/path/to/your/project"

env_folder      <- file.path(project_root, "data/env/current")
terrain_folder  <- file.path(project_root, "data/env/terrain")
occurrence_folder_anemone <- file.path(project_root, "data/occurrence/anemone")
save_location   <- file.path(project_root, "data/log/vif_analysis") # Specific subfolder for VIF results
output_prefix_anemone   <- "anemone_current"

# Define the initial list of variables to consider for *this* scenario
# Note: Variables not found in loaded data or with zero variance will be dropped.
initial_vars_anemone <- c(
  # "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_mean",
  # "thetao_baseline_depthmax_range",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  "chl_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  # "ph_baseline_depthmax_mean",
  # Terrain Vars
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)


# --- Global Parameters ---

coral_shapefile <- file.path(project_root, "data/shapefiles/WCMC008_CoralReef2018_Py_v4_1.shp")
occurrence_crs  <- "EPSG:4326"
env_crs         <- "EPSG:4326"
vif_threshold   <- 10 # For plotting interpretation line


# --- Source the Core Function ---
# Use absolute path for sourcing if running script from different directories
core_script_path <- file.path(project_root, "R/vif_analysis_core.R")
source(core_script_path)


# --- Run the Analysis ---
run_vif_analysis(
  env_folder             = env_folder,
  terrain_folder         = terrain_folder,
  coral_shapefile        = coral_shapefile,
  occurrence_folder      = occurrence_folder_anemone,
  initial_selected_vars  = initial_vars_anemone,
  save_location          = save_location,
  output_prefix          = output_prefix_anemone,
  occurrence_crs         = occurrence_crs,
  env_crs                = env_crs,
  vif_threshold          = vif_threshold
)


initial_vars_anemone <- c(
  # "par_ssp119_depthsurf_dec50_mean",
  "sws_ssp119_depthsurf_dec50_mean",
  "thetao_ssp119_depthmax_dec50_mean",
  # "thetao_ssp119_depthmax_dec50_range",
  "so_ssp119_depthmax_dec50_mean",
  "no3_ssp119_depthmax_dec50_mean",
  "no3_ssp119_depthmax_dec50_range",
  "chl_ssp119_depthmax_dec50_mean",
  "o2_ssp119_depthmax_dec50_range",
  # "ph_ssp119_depthmax_dec50_mean",
  # Terrain Vars
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
env_folder      <- file.path(project_root, "data/env/future/ssp119/dec50")
output_prefix_anemone   <- "anemone_ssp119_dec50"
run_vif_analysis(
  env_folder             = env_folder,
  terrain_folder         = terrain_folder,
  coral_shapefile        = coral_shapefile,
  occurrence_folder      = occurrence_folder_anemone,
  initial_selected_vars  = initial_vars_anemone,
  save_location          = save_location,
  output_prefix          = output_prefix_anemone,
  occurrence_crs         = occurrence_crs,
  env_crs                = env_crs,
  vif_threshold          = vif_threshold
)


initial_vars_anemone <- c(
  # "par_ssp119_depthsurf_dec100_mean",
  "sws_ssp119_depthsurf_dec100_mean",
  "thetao_ssp119_depthmax_dec100_mean",
  # "thetao_ssp119_depthmax_dec100_range",
  "so_ssp119_depthmax_dec100_mean",
  "no3_ssp119_depthmax_dec100_mean",
  "no3_ssp119_depthmax_dec100_range",
  "chl_ssp119_depthmax_dec100_mean",
  "o2_ssp119_depthmax_dec100_range",
  # "ph_ssp119_depthmax_dec100_mean",
  # Terrain Vars
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
env_folder      <- file.path(project_root, "data/env/future/ssp119/dec100")
output_prefix_anemone   <- "anemone_ssp119_dec100"
run_vif_analysis(
  env_folder             = env_folder,
  terrain_folder         = terrain_folder,
  coral_shapefile        = coral_shapefile,
  occurrence_folder      = occurrence_folder_anemone,
  initial_selected_vars  = initial_vars_anemone,
  save_location          = save_location,
  output_prefix          = output_prefix_anemone,
  occurrence_crs         = occurrence_crs,
  env_crs                = env_crs,
  vif_threshold          = vif_threshold
)


initial_vars_anemone <- c(
  # "par_ssp585_depthsurf_dec50_mean",
  "sws_ssp585_depthsurf_dec50_mean",
  "thetao_ssp585_depthmax_dec50_mean",
  # "thetao_ssp585_depthmax_dec50_range",
  "so_ssp585_depthmax_dec50_mean",
  "no3_ssp585_depthmax_dec50_mean",
  "no3_ssp585_depthmax_dec50_range",
  "chl_ssp585_depthmax_dec50_mean",
  "o2_ssp585_depthmax_dec50_range",
  # "ph_ssp585_depthmax_dec50_mean",
  # Terrain Vars
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
env_folder      <- file.path(project_root, "data/env/future/ssp585/dec50")
output_prefix_anemone   <- "anemone_ssp585_dec50"
run_vif_analysis(
  env_folder             = env_folder,
  terrain_folder         = terrain_folder,
  coral_shapefile        = coral_shapefile,
  occurrence_folder      = occurrence_folder_anemone,
  initial_selected_vars  = initial_vars_anemone,
  save_location          = save_location,
  output_prefix          = output_prefix_anemone,
  occurrence_crs         = occurrence_crs,
  env_crs                = env_crs,
  vif_threshold          = vif_threshold
)


initial_vars_anemone <- c(
  # "par_ssp585_depthsurf_dec100_mean",
  "sws_ssp585_depthsurf_dec100_mean",
  "thetao_ssp585_depthmax_dec100_mean",
  # "thetao_ssp585_depthmax_dec100_range",
  "so_ssp585_depthmax_dec100_mean",
  "no3_ssp585_depthmax_dec100_mean",
  "no3_ssp585_depthmax_dec100_range",
  "chl_ssp585_depthmax_dec100_mean",
  "o2_ssp585_depthmax_dec100_range",
  # "ph_ssp585_depthmax_dec100_mean",
  # Terrain Vars
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
env_folder      <- file.path(project_root, "data/env/future/ssp585/dec100")
output_prefix_anemone   <- "anemone_ssp585_dec100"
run_vif_analysis(
  env_folder             = env_folder,
  terrain_folder         = terrain_folder,
  coral_shapefile        = coral_shapefile,
  occurrence_folder      = occurrence_folder_anemone,
  initial_selected_vars  = initial_vars_anemone,
  save_location          = save_location,
  output_prefix          = output_prefix_anemone,
  occurrence_crs         = occurrence_crs,
  env_crs                = env_crs,
  vif_threshold          = vif_threshold
)


cat("\nAnemone VIF analysis script finished.\n")
