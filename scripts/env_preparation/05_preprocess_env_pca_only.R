# scripts/env_preparation/05_preprocess_env_pca_only.R
#-------------------------------------------------------------------------------
# Preprocess Environmental Data using PCA Only - v5 (Jiménez et al. Replication)
#
# METHODOLOGY:
# 1. Loads Current Environmental Data.
# 2. Applies the "Jiménez Physiological Mask":
#    - Latitude/Longitude: Crop to Indo-Pacific BBox (-180 to 180, -40 to 45).
#    - Depth: > -200m (Epipelagic Zone).
#    - Temperature: > 20°C (Warm Water Zone).
# 3. Performs PCA on 100,000 points sampled ONLY from this masked physiological zone.
# 4. Projects this PCA onto Current and Future scenarios (keeping the same BBox).
#-------------------------------------------------------------------------------

cat("--- Running Script 05: PCA Preprocessing (Jiménez et al. Methodology) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, stats, tools, stringr)

# Load Helper Functions
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
if (!file.exists(env_helper_path)) stop("Env Helper missing at: ", env_helper_path)
source(env_helper_path)

# --- 2. Define Target Variables ---
# These must match the "Current" layer names in your raw data
target_env_vars_current_names <- c(
  "sws_baseline_depthsurf_mean", 
  "so_baseline_depthmax_mean", 
  "thetao_baseline_depthmax_mean", 
  "no3_baseline_depthmax_mean", 
  "no3_baseline_depthmax_range", 
  "chl_baseline_depthmax_mean", 
  "o2_baseline_depthmax_range", 
  "phyc_baseline_depthmax_mean", 
  "rugosity"
)

cat("--- Target variables for PCA:", paste(target_env_vars_current_names, collapse=", "), "---\n")

# --- 3. Load & Mask Current Environmental Data ---
cat("--- Loading and Masking Current Environmental Data ---\n")

# A. Load the Full Stack
current_stack_raw <- load_stack_env_data("current", config)
if(is.null(current_stack_raw)) stop("Failed to load current environmental stack.")

# B. Subset to Target Variables
if(!all(target_env_vars_current_names %in% names(current_stack_raw))) {
  missing <- setdiff(target_env_vars_current_names, names(current_stack_raw))
  stop("Missing target variables in current stack: ", paste(missing, collapse=", "))
}
pca_input_stack <- current_stack_raw[[target_env_vars_current_names]]
rm(current_stack_raw); gc()

# C. Apply Bounding Box Crop (Indo-Pacific Wide)
# defined in config: xmin=-180, xmax=180, ymin=-40, ymax=45
if(config$apply_indo_pacific_crop) {
  cat("  Cropping to Config BBox (Lat: ", config$indo_pacific_bbox['ymin'], " to ", config$indo_pacific_bbox['ymax'], ")\n")
  bbox_ext <- terra::ext(config$indo_pacific_bbox)
  pca_input_stack <- terra::crop(pca_input_stack, bbox_ext)
}

# D. Create the Physiological Mask (Jiménez Logic)
cat("  Creating Physiological Mask (Depth >", config$depth_min, "m & Temp >", config$pca_temp_min_threshold, "C)...\n")

# Load Bathymetry (Depth)
if(!file.exists(config$bathymetry_file)) stop("Bathymetry file missing: ", config$bathymetry_file)
bathy_rast <- terra::rast(config$bathymetry_file)
bathy_rast <- terra::crop(bathy_rast, pca_input_stack) # Align extent

# --- FIX FOR compareGeom ERROR ---
# We check geometry robustly. If check fails, we assume resampling is needed.
geoms_match <- tryCatch({
  # Try basic comparison. Returns TRUE if match, FALSE/Error if not.
  # We suppress warnings to avoid console clutter if attributes differ.
  terra::compareGeom(bathy_rast, pca_input_stack, stopiffalse = FALSE, warn = FALSE)
}, error = function(e) {
  return(FALSE) # Assume mismatch on error
})

if(!geoms_match) {
  cat("    Resampling bathymetry to match environmental grid...\n")
  bathy_rast <- terra::resample(bathy_rast, pca_input_stack, method="bilinear")
}
# ---------------------------------

# Get Temperature Layer (for the > 20C filter)
# usually "thetao_baseline_depthmax_mean" or "thetao_baseline_depthsurf_mean"
temp_var_name <- grep("thetao.*mean", target_env_vars_current_names, value=TRUE)[1]
if(is.na(temp_var_name)) stop("Could not auto-detect a temperature variable for masking.")
temp_rast <- pca_input_stack[[temp_var_name]]

# --- THE MASK CALCULATION ---
# 1. Depth Mask: Shallower than -200m (values > -200)
mask_depth <- bathy_rast > config$depth_min 

# 2. Temp Mask: Warmer than 20C
mask_temp <- temp_rast > config$pca_temp_min_threshold

# 3. Combine
mask_physiological <- mask_depth & mask_temp
mask_physiological[mask_physiological == 0] <- NA # Set FALSE to NA for masking

# E. Apply Mask to Stack
cat("  Applying Physiological Mask to stack...\n")
pca_stack_masked <- terra::mask(pca_input_stack, mask_physiological)

# F. Trim (Remove empty outer margins to save space)
pca_stack_final <- terra::trim(pca_stack_masked)

cat("  Final Study Area defined. Dimensions:", paste(dim(pca_stack_final), collapse="x"), "\n")

# Cleanup
rm(pca_input_stack, bathy_rast, temp_rast, mask_depth, mask_temp, mask_physiological); gc()


# --- 4. PCA Training ---
cat("--- Training PCA Model ---\n")

# Sample background points from the MASKED stack
# This ensures PCA learns from the "Available" niche, not the deep ocean.
n_points <- config$pca_background_points_n
cat("  Sampling", n_points, "points from the physiological zone...\n")

set.seed(config$global_seed)
bg_points <- terra::spatSample(pca_stack_final, size = n_points, method = "random", na.rm = TRUE, xy = FALSE, warn = FALSE)

if(nrow(bg_points) < 1000) stop("PCA Sampling failed! Too few valid pixels. Check your depth/temp thresholds.")

cat("  Performing PCA calculation on", nrow(bg_points), "points...\n")
# Scale and Run PCA
pca_model <- stats::prcomp(bg_points, center = TRUE, scale. = TRUE)

# Save Model
dir.create(dirname(config$pca_model_save_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(pca_model, config$pca_model_save_path)
cat("  PCA Model saved to:", config$pca_model_save_path, "\n")

# Save Variance Plot (Optional)
png(file.path(config$log_dir_base, "pca_variance_plot.png"))
plot(pca_model, type = "l", main = "PCA Variance Explained")
dev.off()


# --- 5. Project PCA onto Scenarios ---
cat("\n--- Projecting PCA onto Current and Future Scenarios ---\n")
pca_raster_paths_list <- list()

# Helper function for projection
project_pca_to_disk <- function(stack_in, model, scen_name, output_dir) {
  
  # Ensure names match what the model expects
  names(stack_in) <- names(model$center)
  
  # Predict
  cat("    Calculating PCA prediction...\n")
  pca_pred <- terra::predict(stack_in, model, index = 1:config$n_pca_components)
  names(pca_pred) <- paste0("PC", 1:config$n_pca_components)
  
  # Save
  fname <- paste0("pca_rasters_", scen_name, "_selected_vars.tif")
  fpath <- file.path(output_dir, fname)
  
  terra::writeRaster(pca_pred, fpath, overwrite = TRUE, gdal=c("COMPRESS=LZW"))
  return(fpath)
}

# A. Project Current (Use the masked stack we already have)
cat("  Projecting 'current'...\n")
pca_raster_paths_list[["current"]] <- project_pca_to_disk(pca_stack_final, pca_model, "current", config$log_dir_base)

# Cleanup Current
rm(pca_stack_final); gc()

# B. Project Future Scenarios
for(scen in config$env_scenarios) {
  if(scen == "current") next
  cat("  Processing scenario:", scen, "...\n")
  
  # 1. Load Future Raw Data
  fut_stack_raw <- load_stack_env_data(scen, config)
  if(is.null(fut_stack_raw)) { warning("Skipping ", scen, " (Load failed)"); next }
  
  # 2. Crop to BBox (Consistent with Current)
  if(config$apply_indo_pacific_crop) {
    bbox_ext <- terra::ext(config$indo_pacific_bbox)
    fut_stack_raw <- terra::crop(fut_stack_raw, bbox_ext)
  }
  
  # 3. Rename Variables to match PCA Model
  # We use the helper to find which future var maps to which current var
  # But we need to subset and rename strictly.
  
  # Generate expected future names based on target current names
  future_names_expected <- generate_scenario_variable_list(target_env_vars_current_names, scen, config)
  
  # Subset
  if(!all(future_names_expected %in% names(fut_stack_raw))) {
    warning("Missing variables in ", scen, ". Skipping.")
    next
  }
  fut_stack_sub <- fut_stack_raw[[future_names_expected]]
  
  # RENAME to match Current (required for predict)
  sorted_stack <- terra::rast()
  for(i in 1:length(target_env_vars_current_names)) {
    curr_name <- target_env_vars_current_names[i]
    fut_name <- generate_scenario_variable_list(curr_name, scen, config)
    layer <- fut_stack_sub[[fut_name]]
    names(layer) <- curr_name # Rename immediately
    add(sorted_stack) <- layer
  }
  
  # 4. Project
  pca_raster_paths_list[[scen]] <- project_pca_to_disk(sorted_stack, pca_model, scen, config$log_dir_base)
  
  rm(fut_stack_raw, fut_stack_sub, sorted_stack); gc()
}

# --- 6. Save Output List ---
saveRDS(pca_raster_paths_list, file = config$pca_raster_paths_rds_path)
cat("--- PCA Processing Complete. Paths saved to:", config$pca_raster_paths_rds_path, "---\n")