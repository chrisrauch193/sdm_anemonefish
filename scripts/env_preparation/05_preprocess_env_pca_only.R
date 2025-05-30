# scripts/env_preparation/05_preprocess_env_pca_only.R
#-------------------------------------------------------------------------------
# Preprocess Environmental Data using PCA Only (Braun et al. 2023 Method) - v3
#
# 1. Loads ALL relevant current environmental + terrain layers.
# 2. **Crops** current layers to Indo-Pacific extent (if configured).
# 3. **Masks** the sampling layer to coral reefs (if configured).
# 4. Samples background points from the cropped/masked current layer.
# 5. Extracts values at points & builds ONE PCA model (using stats::prcomp).
# 6. Projects this PCA model onto the CURRENT and ALL FUTURE environmental stacks
#    (cropping future stacks first, if configured).
# 7. Saves the PCA raster stacks and paths for use in SDM scripts.
#-------------------------------------------------------------------------------

cat("--- Running Script 05b: PCA Preprocessing (v3 - IndoPacific Crop & Coral BG Masking) ---\n")

# --- 1. Setup ---
# rm(list = ls()) # Optional: Clear workspace
# gc()
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, stats, tools, stringr)

# Source helpers
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
if (!file.exists(env_helper_path)) stop("Env Helper missing.")
source(env_helper_path)

# --- 2. Load CURRENT Env Data for PCA Model Building ---
cat("--- Loading ALL relevant current environmental layers for PCA ---\n")
current_scenario_name <- "current"
pca_variables_stack_current_raw <- load_stack_env_data(current_scenario_name, config) # Load raw

if(is.null(pca_variables_stack_current_raw) || terra::nlyr(pca_variables_stack_current_raw) < 2) {
  stop("Failed to load sufficient current environmental layers for PCA.")
}
cat("  Loaded raw current environmental stack with layers:", paste(names(pca_variables_stack_current_raw), collapse=", "), "\n")

# --- 2b. Apply Indo-Pacific Crop (if configured) ---
pca_variables_stack_current_cropped <- pca_variables_stack_current_raw # Start with raw
if (config$apply_indo_pacific_crop) {
  cat("  Applying Indo-Pacific crop (", paste(config$indo_pacific_bbox, collapse=", "), ")...\n")
  ip_extent <- terra::ext(config$indo_pacific_bbox)
  pca_variables_stack_current_cropped <- tryCatch({
    terra::crop(pca_variables_stack_current_raw, ip_extent)
  }, error = function(e) {
    warning("Failed to crop current stack to Indo-Pacific extent: ", e$message, call. = FALSE)
    pca_variables_stack_current_raw # Return uncropped if error
  })
  # Optional: Check if cropping actually changed the extent
  if(terra::ext(pca_variables_stack_current_cropped) != terra::ext(pca_variables_stack_current_raw)) {
    cat("  Current stack successfully cropped.\n")
  } else {
    cat("  Current stack extent unchanged after cropping (already within or error occurred).\n")
  }
} else {
  cat("  Skipping Indo-Pacific crop for current stack.\n")
}
rm(pca_variables_stack_current_raw); gc() # Remove raw stack

# Store the ORIGINAL names used to build the PCA model
pca_variable_names_original <- names(pca_variables_stack_current_cropped)

# --- 3. Create Sampling Mask & Sample Background Points ---
n_pca_points <- config$pca_background_points_n
cat("--- Preparing sampling area and sampling", n_pca_points, "points for PCA ---\n")

# Start with the first layer of the (potentially cropped) stack
sampling_layer <- pca_variables_stack_current_cropped[[1]]
sampling_mask <- NULL # Initialize mask

# Apply Coral Mask (if configured)
if (config$mask_background_points_to_coral && config$apply_coral_mask) {
  cat("  Applying coral reef mask for background point sampling...\n")
  if (!is.null(config$coral_shapefile) && file.exists(config$coral_shapefile)) {
    coral_areas_sf <- tryCatch({ sf::st_read(config$coral_shapefile, quiet = TRUE) }, error = function(e) {warning(" Failed load coral shapefile:", e$message, call.=FALSE); NULL})
    if (!is.null(coral_areas_sf)) {
      coral_areas_vect <- tryCatch({ terra::vect(coral_areas_sf) }, error = function(e) {warning(" Failed convert coral sf to vect:", e$message, call.=FALSE); NULL})
      if (!is.null(coral_areas_vect)) {
        if(terra::crs(coral_areas_vect) != terra::crs(sampling_layer)){
          cat("    Projecting coral shapefile CRS to match raster...\n")
          coral_areas_vect <- tryCatch(terra::project(coral_areas_vect, terra::crs(sampling_layer)), error = function(e){warning(" Failed project coral shapefile:", e$message, call.=FALSE); NULL})
        }
        if(!is.null(coral_areas_vect)){
          # Mask the sampling layer itself
          sampling_mask <- tryCatch(terra::mask(sampling_layer, coral_areas_vect), error=function(e){warning(" Failed mask sampling layer:", e$message, call.=FALSE); NULL})
          if(!is.null(sampling_mask)) {
            cat("  Coral reef mask applied for sampling.\n")
          } else {
            warning("  Failed to create sampling mask. Sampling from unmasked layer.", call.=FALSE)
            sampling_mask <- sampling_layer # Fallback
          }
        } else { warning("  CRS projection failed. Sampling from unmasked layer.", call.=FALSE); sampling_mask <- sampling_layer }
      } else { warning("  sf to vect conversion failed. Sampling from unmasked layer.", call.=FALSE); sampling_mask <- sampling_layer }
    } else { warning("  Coral shapefile loading failed. Sampling from unmasked layer.", call.=FALSE); sampling_mask <- sampling_layer }
  } else { warning("  Coral shapefile path missing/invalid. Sampling from unmasked layer.", call.=FALSE); sampling_mask <- sampling_layer }
} else {
  cat("  Background point sampling not masked to coral reefs.\n")
  sampling_mask <- sampling_layer # Use the unmasked (but potentially cropped) layer
}

# Sample points from the (potentially masked) layer
set.seed(123)
bg_points_for_pca <- terra::spatSample(sampling_mask, size = n_pca_points, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
if (nrow(bg_points_for_pca) < n_pca_points) { warning("Sampled fewer points (", nrow(bg_points_for_pca), ") than requested (", n_pca_points, ") for PCA (likely due to mask).", call.=FALSE) }
if (nrow(bg_points_for_pca) == 0) stop("Failed to sample any valid points for PCA.")
cat("  Sampled", nrow(bg_points_for_pca), "points from the designated sampling area.\n")
rm(sampling_layer, sampling_mask); gc() # Clean up mask layer

# Extract values from the FULL cropped stack at these points
cat("  Extracting values at sampled points for PCA model building...\n")
bg_values_df <- terra::extract(pca_variables_stack_current_cropped, bg_points_for_pca[, c("x","y")], ID = FALSE)
bg_values_df <- na.omit(bg_values_df)

if (nrow(bg_values_df) < config$n_pca_components) { stop("Insufficient non-NA data rows (", nrow(bg_values_df), ") for PCA with ", config$n_pca_components, " components.") }
cat("  Using", nrow(bg_values_df), "complete cases for PCA model.\n")
rm(pca_variables_stack_current_cropped, bg_points_for_pca); gc() # Clean up the full stack

# --- 4. Perform PCA on CURRENT Data ---
# (PCA model building, saving, plotting code remains the same)
cat("--- Performing PCA on scaled current environmental data ---\n")
bg_values_scaled <- scale(bg_values_df, center = TRUE, scale = TRUE)
pca_model <- tryCatch({ stats::prcomp(bg_values_scaled) }, error = function(e) { warning("PCA calculation failed: ", e$message, call. = FALSE); NULL })
if(is.null(pca_model)) stop("PCA failed.")
rm(bg_values_df, bg_values_scaled); gc()

pca_model_save_path <- config$pca_models_rds_path
dir.create(dirname(pca_model_save_path), showWarnings = FALSE, recursive = TRUE)
tryCatch({ saveRDS(pca_model, file = pca_model_save_path); cat("  PCA model object saved to:", pca_model_save_path, "\n")}, error = function(e) { stop("Failed to save PCA model object: ", e$message) })

# PCA Summary and Plots
pca_summary <- summary(pca_model); print(pca_summary)
pca_variance_path <- file.path(config$log_dir_base, "pca_variance_explained_plot.png")
tryCatch({ png(pca_variance_path, width=8, height=6, units="in", res=300); plot(pca_model, type = "l", main = "PCA Variance Explained"); abline(h=0.9, col="red", lty=2); text(x = length(pca_model$sdev), y = 0.95, labels = "90% Variance", col = "red", pos = 1); dev.off(); cat("  PCA variance explained plot saved to:", pca_variance_path, "\n")}, error = function(e) { warning("Failed save PCA variance plot: ", e$message, call.=FALSE)})
pca_biplot_path <- file.path(config$log_dir_base, "pca_biplot_PC1_PC2.png")
tryCatch({ png(pca_biplot_path, width=8, height=8, units="in", res=300); biplot(pca_model, choices = 1:2, scale = 0, cex = 0.7); title("PCA Biplot (PC1 vs PC2)"); dev.off(); cat("  PCA biplot saved to:", pca_biplot_path, "\n")}, error = function(e) { warning("Failed save PCA biplot: ", e$message, call.=FALSE)})

# --- 5. Project PCA onto Rasters for ALL Scenarios (With Cropping) ---
cat("\n--- Projecting PCA model onto rasters for all scenarios ---\n")
pca_raster_paths_list <- list()

for(scenario in config$env_scenarios) {
  cat("  -- Projecting scenario:", scenario, "--\n")
  
  # Load the *full* environmental stack for this specific scenario
  scenario_stack_raw <- load_stack_env_data(scenario, config)
  if(is.null(scenario_stack_raw)){ warning("Could not load stack for scenario '", scenario, "'. Skipping projection.", call.=FALSE); next }
  
  # Apply Indo-Pacific Crop (if configured)
  scenario_stack_cropped <- scenario_stack_raw # Start with raw
  if (config$apply_indo_pacific_crop) {
    cat("    Applying Indo-Pacific crop...\n")
    ip_extent <- terra::ext(config$indo_pacific_bbox)
    scenario_stack_cropped <- tryCatch({ terra::crop(scenario_stack_raw, ip_extent) },
                                       error = function(e) { warning(" Failed crop stack: ", e$message, call. = FALSE); scenario_stack_raw })
    if(terra::ext(scenario_stack_cropped) != terra::ext(scenario_stack_raw)) { cat("    Stack cropped.\n")} else { cat("    Stack extent unchanged after crop (or error occurred).\n")}
  } else { cat("    Skipping Indo-Pacific crop.\n") }
  rm(scenario_stack_raw); gc() # Remove raw stack
  
  # *** Generate the EXPECTED variable names for THIS scenario ***
  expected_vars_for_scenario <- generate_scenario_variable_list(pca_variable_names_original, scenario, config)
  
  # Check if all EXPECTED variables are present in the loaded (and possibly cropped) stack
  if (!all(expected_vars_for_scenario %in% names(scenario_stack_cropped))) {
    missing_projection_vars <- expected_vars_for_scenario[!expected_vars_for_scenario %in% names(scenario_stack_cropped)]
    warning("Expected variables for PCA projection are MISSING in the loaded/cropped stack for scenario '", scenario, "':\n  ",
            paste(missing_projection_vars, collapse=", "), "\nCannot project PCA. Skipping scenario.", call.=FALSE)
    print(paste("Available layers in cropped stack:", paste(names(scenario_stack_cropped), collapse=", ")))
    rm(scenario_stack_cropped); gc(); next
  }
  
  # Subset and reorder layers using the EXPECTED future names
  scenario_stack_subset <- scenario_stack_cropped[[expected_vars_for_scenario]]
  rm(scenario_stack_cropped); gc()
  
  # *** RENAME the subsetted stack layers back to the ORIGINAL names ***
  if (terra::nlyr(scenario_stack_subset) == length(pca_variable_names_original)) {
    names(scenario_stack_subset) <- pca_variable_names_original
    cat("    Renamed scenario stack layers to match original PCA variable names.\n")
  } else {
    warning("Mismatch in layer count after subsetting for scenario '", scenario, "'. Cannot reliably rename for projection. Skipping.", call. = FALSE)
    rm(scenario_stack_subset); gc(); next
  }
  
  # Project using terra::predict
  pca_raster_scenario <- tryCatch({
    terra::predict(scenario_stack_subset, pca_model, index = 1:config$n_pca_components)
  }, error = function(e){
    warning("terra::predict failed for scenario '", scenario, "': ", e$message, call.=FALSE)
    return(NULL)
  })
  
  if (is.null(pca_raster_scenario)) {
    warning("PCA projection failed for scenario '", scenario, "'. Skipping.", call.=FALSE)
    rm(scenario_stack_subset); gc()
    next
  }
  
  # Rename projected layers
  names(pca_raster_scenario) <- paste0("PC", 1:config$n_pca_components)
  
  # Define output path
  pca_raster_filename <- paste0("pca_rasters_", scenario, ".tif")
  pca_raster_save_path <- file.path(config$log_dir_base, pca_raster_filename)
  
  # Save the PCA raster stack
  tryCatch({
    terra::writeRaster(pca_raster_scenario, filename = pca_raster_save_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    cat("    PCA raster stack saved for scenario '", scenario, "' to:", pca_raster_save_path, "\n")
    pca_raster_paths_list[[scenario]] <- pca_raster_save_path # Store the path
  }, error = function(e) {
    warning("Failed to write PCA raster for scenario '", scenario, "': ", e$message, call.=FALSE)
  })
  
  rm(scenario_stack_subset, pca_raster_scenario); gc()
} # End scenario loop

# --- 6. Save the List of PCA Raster Paths ---
pca_paths_rds_save_path <- config$pca_raster_paths_rds_path
if (length(pca_raster_paths_list) > 0) {
  tryCatch({ saveRDS(pca_raster_paths_list, file = pca_paths_rds_save_path)
    cat("--- List of PCA raster paths saved to:", pca_paths_rds_save_path, "---\n") },
    error = function(e){ warning("Failed to save PCA raster paths list: ", e$message, call.=FALSE) })
} else { warning("No PCA rasters generated/saved. Paths list empty.", call.=FALSE) }

cat("--- Script 05b finished. PCA model built and projected onto (potentially cropped) scenarios. BG points masked if configured. ---\n")
#-------------------------------------------------------------------------------