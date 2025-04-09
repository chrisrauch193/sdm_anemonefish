# scripts/05b_preprocess_env_pca_only.R
#-------------------------------------------------------------------------------
# Preprocess Environmental Data using PCA Only (Braun et al. 2023 Method) - v2
#
# 1. Loads ALL relevant current environmental + terrain layers.
# 2. Samples background points from the current stack.
# 3. Extracts values at points & builds ONE PCA model (using stats::prcomp).
# 4. Projects this PCA model onto the CURRENT and ALL FUTURE environmental stacks,
#    handling the necessary renaming of future layers for projection.
# 5. Saves the PCA raster stacks and paths for use in SDM scripts (06a, 06b, 06d).
#-------------------------------------------------------------------------------

cat("--- Running Script 05b: PCA Preprocessing (Paper Method - v2 Corrected Projection) ---\n")

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
pca_variables_stack_current <- load_stack_env_data(current_scenario_name, config)

if(is.null(pca_variables_stack_current) || terra::nlyr(pca_variables_stack_current) < 2) {
  stop("Failed to load sufficient current environmental layers for PCA.")
}
cat("  Loaded current environmental stack with layers:", paste(names(pca_variables_stack_current), collapse=", "), "\n")
# *** Store the ORIGINAL names used to build the PCA model ***
pca_variable_names_original <- names(pca_variables_stack_current)

# --- 3. Sample Background Points & Extract Values from CURRENT Data ---
n_pca_points <- config$background_points_n # Use value from config
cat("--- Sampling", n_pca_points, "random points from current env stack for PCA ---\n")
set.seed(123)
bg_points_for_pca <- terra::spatSample(pca_variables_stack_current[[1]], size = n_pca_points, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
if (nrow(bg_points_for_pca) < n_pca_points) { warning("Sampled fewer points (", nrow(bg_points_for_pca), ") than requested (", n_pca_points, ") for PCA.") }
if (nrow(bg_points_for_pca) == 0) stop("Failed to sample any valid points for PCA.")

cat("  Extracting values at", nrow(bg_points_for_pca), "points for PCA model building...\n")
bg_values_df <- terra::extract(pca_variables_stack_current, bg_points_for_pca[, c("x","y")], ID = FALSE)
bg_values_df <- na.omit(bg_values_df)

if (nrow(bg_values_df) < config$n_pca_components) { stop("Insufficient non-NA data rows (", nrow(bg_values_df), ") for PCA with ", config$n_pca_components, " components.") }
cat("  Using", nrow(bg_values_df), "complete cases for PCA model.\n")
rm(pca_variables_stack_current, bg_points_for_pca); gc() # Clean up current stack

# --- 4. Perform PCA on CURRENT Data ---
cat("--- Performing PCA on scaled current environmental data ---\n")
bg_values_scaled <- scale(bg_values_df, center = TRUE, scale = TRUE)
pca_model <- tryCatch({ stats::prcomp(bg_values_scaled) }, error = function(e) { warning("PCA calculation failed: ", e$message, call. = FALSE); NULL })
if(is.null(pca_model)) stop("PCA failed.")
rm(bg_values_df, bg_values_scaled); gc() # Clean up scaled data

# --- Save PCA Model Object ---
pca_model_save_path <- config$pca_models_rds_path
dir.create(dirname(pca_model_save_path), showWarnings = FALSE, recursive = TRUE)
tryCatch({ saveRDS(pca_model, file = pca_model_save_path); cat("  PCA model object saved to:", pca_model_save_path, "\n")},
         error = function(e) { stop("Failed to save PCA model object: ", e$message) })

# --- PCA Summary and Plots ---
# (Plotting code remains the same)
pca_summary <- summary(pca_model)
print(pca_summary)
pca_variance_path <- file.path(config$log_dir_base, "pca_variance_explained_plot.png") # Use log_dir_base
tryCatch({ png(pca_variance_path, width=8, height=6, units="in", res=300); plot(pca_model, type = "l", main = "PCA Variance Explained"); abline(h=0.9, col="red", lty=2); text(x = length(pca_model$sdev), y = 0.95, labels = "90% Variance", col = "red", pos = 1); dev.off(); cat("  PCA variance explained plot saved to:", pca_variance_path, "\n")}, error = function(e) { warning("Failed save PCA variance plot: ", e$message, call.=FALSE)})
pca_biplot_path <- file.path(config$log_dir_base, "pca_biplot_PC1_PC2.png") # Use log_dir_base
tryCatch({ png(pca_biplot_path, width=8, height=8, units="in", res=300); biplot(pca_model, choices = 1:2, scale = 0, cex = 0.7); title("PCA Biplot (PC1 vs PC2)"); dev.off(); cat("  PCA biplot saved to:", pca_biplot_path, "\n")}, error = function(e) { warning("Failed save PCA biplot: ", e$message, call.=FALSE)})


# --- 5. Project PCA onto Rasters for ALL Scenarios (Corrected Logic) ---
cat("\n--- Projecting PCA model onto rasters for all scenarios ---\n")
pca_raster_paths_list <- list()

for(scenario in config$env_scenarios) {
  cat("  -- Projecting scenario:", scenario, "--\n")
  
  # Load the *full* environmental stack for this specific scenario
  scenario_stack_raw <- load_stack_env_data(scenario, config)
  if(is.null(scenario_stack_raw)){ warning("Could not load stack. Skipping projection.", call.=FALSE); next }
  
  # *** Generate the EXPECTED variable names for THIS scenario ***
  # This uses the helper to convert the original PCA names (baseline) to future names
  expected_vars_for_scenario <- generate_scenario_variable_list(pca_variable_names_original, scenario, config)
  
  # Check if all EXPECTED future variables are present in the loaded stack
  if (!all(expected_vars_for_scenario %in% names(scenario_stack_raw))) {
    missing_projection_vars <- expected_vars_for_scenario[!expected_vars_for_scenario %in% names(scenario_stack_raw)]
    warning("Expected variables for PCA projection are MISSING in the loaded stack for scenario '", scenario, "':\n  ",
            paste(missing_projection_vars, collapse=", "), "\nCannot project PCA. Skipping scenario.", call.=FALSE)
    print(paste("Available layers in stack:", paste(names(scenario_stack_raw), collapse=", ")))
    rm(scenario_stack_raw); gc(); next
  }
  
  # Subset and reorder layers using the EXPECTED future names
  scenario_stack_subset <- scenario_stack_raw[[expected_vars_for_scenario]]
  rm(scenario_stack_raw); gc()
  
  # *** RENAME the subsetted stack layers back to the ORIGINAL names used for PCA model training ***
  # This is crucial for terra::predict to work correctly with the prcomp object
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
  
  # Define output path (saving to log_dir_base as before)
  pca_raster_filename <- paste0("pca_rasters_", scenario, ".tif")
  pca_raster_save_path <- file.path(config$log_dir_base, pca_raster_filename) # Use log_dir_base from config
  
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
pca_paths_rds_save_path <- config$pca_raster_paths_rds_path # Path from config
if (length(pca_raster_paths_list) > 0) {
  tryCatch({ saveRDS(pca_raster_paths_list, file = pca_paths_rds_save_path)
    cat("--- List of PCA raster paths saved to:", pca_paths_rds_save_path, "---\n") },
    error = function(e){ warning("Failed to save PCA raster paths list: ", e$message, call.=FALSE) })
} else { warning("No PCA rasters generated/saved. Paths list empty.", call.=FALSE) }

cat("--- Script 05b finished. PCA model built and projected onto all scenarios. ---\n")
#-------------------------------------------------------------------------------