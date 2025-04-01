# scripts/04_download_env_data.R
#-------------------------------------------------------------------------------
# Download Environmental Data using logic from the original script.
# Uses obissdm::get_env_data (ensure package is installed).
# Uses variable lists and paths defined in config.R.
# Also prepares terrain layers (rugosity, distance to coast).
# Assumes obissdm downloads to data/env/current, data/env/future/sspXXX, data/env/terrain
#-------------------------------------------------------------------------------
cat("--- Running Script 04: Download Environmental Data (User's Original Logic) ---\n")

# Ensure config is loaded if running standalone
if (!exists("config")) {
  source("scripts/config.R")
  if (!exists("config")) {
    stop("Failed to load config object from scripts/config.R")
  }
}

# Ensure required packages are loaded
pacman::p_load(terra, fs)
if (!requireNamespace("obissdm", quietly = TRUE)) {
  warning("Package 'obissdm' not found. Attempting to install from GitHub...\n",
          "You might need to run: remotes::install_github('iobis/obissdm') manually.", call. = FALSE)
  tryCatch({
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("iobis/obissdm")
    library(obissdm)
  }, error = function(e){
    stop("Failed to install or load 'obissdm'. Please install it manually. Error: ", e$message)
  })
} else {
  library(obissdm)
}


# --- Parameters Directly from Config ---
datasets_mean <- config$original_datasets_mean
range_datasets <- config$original_range_datasets
future_scenarios <- config$original_future_scenarios
time_steps <- config$original_time_steps
terrain_vars_for_download <- config$original_terrain_vars_for_download

cat("Using datasets for mean download from config:", paste(datasets_mean, collapse=", "), "\n")
cat("Using datasets for range/min/max from config:", paste(range_datasets, collapse=", "), "\n")
cat("Using future scenarios from config:", paste(future_scenarios, collapse=", "), "\n")
cat("Using time steps for future data from config:", paste(names(time_steps), collapse=", "), "\n")
cat("Requesting terrain variables from config:", paste(terrain_vars_for_download, collapse=", "), "\n")


# Variables to download (hardcoded as in original script)
variables_mean_only <- c("mean")
variables_range_min_max <- c("range", "ltmin", "ltmax")


# --- Download Logic (Matching Original Script Calls - NO directory/skip_existing) ---
# Assumes download happens into data/env/current, data/env/future/sspXXX, data/env/terrain

# 1. Download Mean Vars + Terrain (First Call)
cat("\nAttempting initial download (Mean vars + Terrain)...\n")
tryCatch({
  obissdm::get_env_data(
    datasets = datasets_mean,
    future_scenarios = future_scenarios,
    time_steps = time_steps,
    variables = variables_mean_only,
    terrain_vars = terrain_vars_for_download,
    average_time = TRUE # As per original call
  )
  cat("Initial download call completed.\n")
}, error = function(e) {
  warning("Error during initial download call: ", e$message, call. = FALSE)
})


# 2. Download Range, Min, Max for specific datasets (Second Call)
cat("\nAttempting range/min/max download...\n")
tryCatch({
  obissdm::get_env_data(
    datasets = range_datasets,
    future_scenarios = future_scenarios,
    time_steps = time_steps,
    variables = variables_range_min_max,
    terrain_vars = NULL, # Don't download terrain again
    average_time = TRUE # As per original call
  )
  cat("Range/min/max download call completed.\n")
}, error = function(e) {
  warning("Error during range/min/max download call: ", e$message, call. = FALSE)
})

cat("\nDownloads attempted. Assuming files were saved by obissdm to the standard structure under '", config$env_data_dir, "'.\n", sep="")

# --- Post-processing Steps (Using paths from config, expecting files there) ---

# 3. Rename terrain_ruggedness
cat("\nProcessing Terrain: Renaming ruggedness to rugosity...\n")
terrain_target_folder <- config$terrain_folder # Path where files should be
ruggedness_download_name <- "terrain_ruggedness_index" # Name downloaded by obissdm
final_rugosity_name <- "rugosity" # Final desired name

if (dir.exists(terrain_target_folder)) {
  # Pattern uses the name expected from download
  pattern_to_find <- paste0(ruggedness_download_name, ".*\\.tif$")
  to_rename <- list.files(terrain_target_folder, pattern = pattern_to_find, recursive = TRUE, full.names = TRUE)
  to_rename <- to_rename[!grepl("\\.aux\\.xml$", to_rename, ignore.case = TRUE)]
  
  if (length(to_rename) > 0) {
    renamed_count <- 0
    for(f in to_rename){
      # Replace download name with desired final name (rugosity)
      new_name <- gsub(ruggedness_download_name, final_rugosity_name, f, ignore.case = TRUE)
      tryCatch({
        edit_r <- terra::rast(f)
        names(edit_r) <- final_rugosity_name # Final layer name
        terra::writeRaster(edit_r, new_name, overwrite = TRUE, gdal=c("COMPRESS=LZW"))
        if(file.exists(new_name) && file.size(new_name) > 0){
          fs::file_delete(f)
          aux_file_old <- paste0(tools::file_path_sans_ext(f), ".aux.xml")
          if(fs::file_exists(aux_file_old)) fs::file_delete(aux_file_old)
          renamed_count <- renamed_count + 1
        } else {
          warning("Failed to write/verify new rugosity file:", new_name, " keeping original:", f)
        }
      }, error = function(e){
        warning("Error processing ruggedness file ", f, ": ", e$message, call. = FALSE)
      })
    }
    cat("Processed", renamed_count, "ruggedness file(s) in target location.\n")
  } else {
    cat("No '", ruggedness_download_name, ".tif' files found in target location: ", terrain_target_folder, "\n", sep="")
  }
} else {
  warning("Target terrain folder not found:", terrain_target_folder, ". Cannot rename rugosity.", call. = FALSE)
}


# 4. Distance to coast layer
cat("\nProcessing Terrain: Calculating distance to coast...\n")
terrain_target_folder <- config$terrain_folder # Use config path
distcoast_file <- file.path(terrain_target_folder, "distcoast.tif") # Output file path

# Check if file exists in the target terrain folder
if (!file.exists(distcoast_file)) {
  # Base layer path expected in the config-defined current scenario folder
  base_layer_filename <- "thetao_baseline_2000_2019_depthsurf_mean.tif" # Check exact name if needed
  base_layer_path <- file.path(config$scenario_folder_map$current, base_layer_filename)
  
  if (file.exists(base_layer_path)) {
    tryCatch({
      base <- terra::rast(base_layer_path)
      cat("Using base raster:", base_layer_path, "\n")
      
      land_mask <- terra::classify(base[[1]], cbind(NA, 1), others = NA) # Land=1, Ocean=NA
      
      if(is.na(terra::global(land_mask, "max", na.rm=TRUE)$max) || terra::global(land_mask, "max", na.rm=TRUE)$max != 1){
        warning("Base layer seems to be all NA or all non-NA after masking. Cannot create land mask for distance.", call.=FALSE)
      } else {
        coast_agg <- terra::aggregate(land_mask, fact = 4, fun = "max", na.rm = TRUE)
        coast_dist_agg <- terra::distance(coast_agg)
        coast_dist <- terra::disagg(coast_dist_agg, fact = 4, method = "bilinear")
        coast_dist <- terra::resample(coast_dist, base, method="bilinear")
        coast_dist <- terra::mask(coast_dist, base[[1]])
        coast_dist_km <- coast_dist / 1000
        names(coast_dist_km) <- "distcoast"
        dir.create(dirname(distcoast_file), showWarnings = FALSE, recursive = TRUE) # Ensure dir exists
        terra::writeRaster(coast_dist_km, distcoast_file, overwrite = TRUE, gdal=c("COMPRESS=LZW"))
        cat("Distance to coast raster saved to:", distcoast_file, "\n")
      }
    }, error = function(e){
      warning("Error calculating distance to coast: ", e$message, call.=FALSE)
    })
  } else {
    warning("Base layer for distance calculation not found in target location: ", base_layer_path, call. = FALSE)
  }
} else {
  cat("Distance to coast file already exists in target location. Skipping calculation.\n")
}

# Verification Step (Optional)
cat("\nVerifying expected directory structure...\n")
# Check current dir
if(dir.exists(config$scenario_folder_map$current)){
  cat("  Found:", config$scenario_folder_map$current, "\n")
} else {
  cat("  MISSING:", config$scenario_folder_map$current, "\n")
}
# Check future base dir
future_base_dir = dirname(config$scenario_folder_map$ssp119_2050) # e.g., data/env/future
if(dir.exists(future_base_dir)){
  cat("  Found:", future_base_dir, "\n")
  # Optionally check subdirs if future scenarios were requested
  if(length(future_scenarios) > 0){
    # Check ssp119 dir
    ssp119_dir = file.path(future_base_dir, "ssp119")
    if(dir.exists(ssp119_dir)) cat("  Found:", ssp119_dir, "\n") else cat("  MISSING:", ssp119_dir, "\n")
    # Check ssp585 dir
    ssp585_dir = file.path(future_base_dir, "ssp585")
    if(dir.exists(ssp585_dir)) cat("  Found:", ssp585_dir, "\n") else cat("  MISSING:", ssp585_dir, "\n")
  }
} else if (length(future_scenarios) > 0) { # Only warn if future requested but base dir missing
  cat("  MISSING base future directory:", future_base_dir, "\n")
}
# Check terrain dir
if(dir.exists(config$terrain_folder)){
  cat("  Found:", config$terrain_folder, "\n")
} else {
  cat("  MISSING:", config$terrain_folder, "\n")
}


cat("\n--- Script 04 finished. ---\n")
#-------------------------------------------------------------------------------

# --- Comments from original script (Retained for reference) ---
# ... (comments remain the same) ...
# --- End Comments ---