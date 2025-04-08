# scripts/05b_preprocess_env_pca_only.R
#-------------------------------------------------------------------------------
# Preprocess Environmental Data using PCA Only (Braun et al. 2023 Method)
#
# 1. Loads ALL relevant current environmental + terrain layers.
# 2. Samples background points from the current stack.
# 3. Extracts values at points & builds ONE PCA model (using stats::prcomp).
# 4. Projects this PCA model onto the CURRENT and ALL FUTURE environmental stacks.
# 5. Saves the PCA raster stacks and paths for use in SDM scripts (06a, 06b).
#-------------------------------------------------------------------------------

cat("--- Running Script 05b: PCA Preprocessing (Paper Method) ---\n")

# --- 1. Setup ---
rm(list = ls()) # Clear workspace to ensure clean run
gc()
if (!exists("config")) {
  # Source the configuration file if it doesn't exist
  source("scripts/config.R")
  if (!exists("config")) {
    stop("Failed to load config object from scripts/config.R")
  }
}
pacman::p_load(terra, sf, dplyr, readr, stats, tools, stringr) # Load necessary packages

# Source helpers (only env_processing needed for this script)
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
if (!file.exists(env_helper_path)) stop("Env Helper missing.")
source(env_helper_path)

# --- 2. Define Variables for PCA ---
# According to the paper, PCA is done on the environmental variables.
# We will load *all* available current environmental layers + terrain layers.
# No VIF pre-selection is performed here, following the paper's method section.
cat("--- Loading ALL relevant current environmental layers for PCA ---\n")
current_scenario_name <- "current"
pca_variables_stack_current <- load_stack_env_data(current_scenario_name, config)

if(is.null(pca_variables_stack_current) || terra::nlyr(pca_variables_stack_current) < 2) {
  stop("Failed to load sufficient current environmental layers for PCA.")
}
cat("  Loaded current environmental stack with layers:", paste(names(pca_variables_stack_current), collapse=", "), "\n")
pca_variable_names <- names(pca_variables_stack_current) # Store the names used

# --- 3. Sample Background Points & Extract Values ---
# Paper used 100,000 points; using config$background_points_n for consistency,
# ensure this is large enough (e.g., >= 10000).
n_pca_points <- 10000 # Or use config$background_points_n if preferred and large enough
cat("--- Sampling", n_pca_points, "random points from current env stack for PCA ---\n")

# Sample points from the *first layer* ensuring they are not NA there.
# We'll extract values from the full stack afterwards.
set.seed(123) # for reproducibility
bg_points_for_pca <- terra::spatSample(pca_variables_stack_current[[1]],
                                       size = n_pca_points,
                                       method = "random",
                                       na.rm = TRUE, # Essential
                                       xy = TRUE,
                                       warn = FALSE) # Suppress warning if fewer points found

if (nrow(bg_points_for_pca) < n_pca_points) {
  warning("Sampled fewer points (", nrow(bg_points_for_pca), ") than requested (", n_pca_points, ") for PCA. Check raster NAs.")
}
if (nrow(bg_points_for_pca) == 0) stop("Failed to sample any valid points for PCA.")

cat("  Extracting values at", nrow(bg_points_for_pca), "points for PCA model building...\n")
bg_values_df <- terra::extract(pca_variables_stack_current, bg_points_for_pca[, c("x","y")], ID = FALSE)
bg_values_df <- na.omit(bg_values_df) # Remove rows with NAs in *any* layer

if (nrow(bg_values_df) < config$n_pca_components) {
  stop("Insufficient non-NA data rows (", nrow(bg_values_df), ") remaining after extraction to perform PCA with ", config$n_pca_components, " components.")
}
cat("  Using", nrow(bg_values_df), "complete cases for PCA model.\n")

# --- 4. Perform PCA ---
cat("--- Performing PCA on scaled current environmental data ---\n")
# Center and scale the data
bg_values_scaled <- scale(bg_values_df, center = TRUE, scale = TRUE)

# Run PCA using stats::prcomp
pca_model <- tryCatch({
  stats::prcomp(bg_values_scaled)
}, error = function(e) {
  warning("PCA calculation failed: ", e$message, call. = FALSE)
  return(NULL)
})

if(is.null(pca_model)) stop("PCA failed.")

# --- Save PCA Model Object ---
pca_model_save_path <- config$pca_models_rds_path # Path from config
dir.create(dirname(pca_model_save_path), showWarnings = FALSE, recursive = TRUE)
tryCatch({
  saveRDS(pca_model, file = pca_model_save_path)
  cat("  PCA model object saved to:", pca_model_save_path, "\n")
}, error = function(e) {
  stop("Failed to save PCA model object: ", e$message)
})

# --- PCA Summary and Plots (Optional but Recommended) ---
pca_summary <- summary(pca_model)
print(pca_summary)

# Plot variance explained
pca_variance_path <- file.path(config$log_dir, "pca_variance_explained_plot.png")
tryCatch({
  png(pca_variance_path, width=8, height=6, units="in", res=300)
  plot(pca_model, type = "l", main = "PCA Variance Explained")
  abline(h=0.9, col="red", lty=2) # Add line for 90% variance if needed
  text(x = length(pca_model$sdev), y = 0.95, labels = "90% Variance", col = "red", pos = 1)
  dev.off()
  cat("  PCA variance explained plot saved to:", pca_variance_path, "\n")
}, error = function(e) { warning("Failed to save PCA variance plot: ", e$message, call.=FALSE)})

# Biplot (might be messy with many vars, select first few components)
pca_biplot_path <- file.path(config$log_dir, "pca_biplot_PC1_PC2.png")
tryCatch({
  png(pca_biplot_path, width=8, height=8, units="in", res=300)
  biplot(pca_model, choices = 1:2, scale = 0, cex = 0.7) # Basic biplot
  title("PCA Biplot (PC1 vs PC2)")
  dev.off()
  cat("  PCA biplot saved to:", pca_biplot_path, "\n")
}, error = function(e) { warning("Failed to save PCA biplot: ", e$message, call.=FALSE)})


# --- 5. Project PCA onto Rasters for ALL Scenarios ---
cat("--- Projecting PCA model onto rasters for all scenarios ---\n")
pca_raster_paths_list <- list() # To store paths of the output PCA rasters

for(scenario in config$env_scenarios) {
  cat("  -- Projecting scenario:", scenario, "--\n")
  
  # Load the *full* environmental stack for this scenario
  scenario_stack_raw <- load_stack_env_data(scenario, config)
  if(is.null(scenario_stack_raw)){
    warning("Could not load stack for scenario '", scenario, "'. Skipping projection.", call.=FALSE)
    next
  }
  
  # *** CRITICAL: Ensure layers match those used for PCA ***
  # Check if all variables used to build the PCA are present in this scenario stack
  if (!all(pca_variable_names %in% names(scenario_stack_raw))) {
    missing_projection_vars <- pca_variable_names[!pca_variable_names %in% names(scenario_stack_raw)]
    warning("Variables used for PCA model building are MISSING in the raster stack for scenario '", scenario, "':\n  ",
            paste(missing_projection_vars, collapse=", "), "\nCannot project PCA. Skipping scenario.", call.=FALSE)
    print(paste("Available layers in stack:", paste(names(scenario_stack_raw), collapse=", ")))
    rm(scenario_stack_raw); gc(); next
  }
  
  # Select and reorder layers to exactly match the PCA model input
  scenario_stack_ordered <- scenario_stack_raw[[pca_variable_names]]
  rm(scenario_stack_raw); gc() # Free memory
  
  # Project using terra::predict and the prcomp object
  # The 'predict' method for SpatRaster and prcomp handles scaling automatically
  # based on the 'center' and 'scale' attributes stored in the pca_model object.
  pca_raster_scenario <- tryCatch({
    terra::predict(scenario_stack_ordered, pca_model, index = 1:config$n_pca_components)
  }, error = function(e){
    warning("terra::predict failed for scenario '", scenario, "': ", e$message, call.=FALSE)
    return(NULL)
  })
  
  if (is.null(pca_raster_scenario)) {
    warning("PCA projection failed for scenario '", scenario, "'. Skipping.", call.=FALSE)
    rm(scenario_stack_ordered); gc()
    next
  }
  
  # Rename layers
  names(pca_raster_scenario) <- paste0("PC", 1:config$n_pca_components)
  
  # Define output path
  pca_raster_filename <- paste0("pca_rasters_", scenario, ".tif") # Simplified name
  pca_raster_save_path <- file.path(config$log_dir, pca_raster_filename) # Save in logs dir
  
  # Save the PCA raster stack
  tryCatch({
    terra::writeRaster(pca_raster_scenario, filename = pca_raster_save_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    cat("    PCA raster stack saved for scenario '", scenario, "' to:", pca_raster_save_path, "\n")
    pca_raster_paths_list[[scenario]] <- pca_raster_save_path # Store the path
  }, error = function(e) {
    warning("Failed to write PCA raster for scenario '", scenario, "': ", e$message, call.=FALSE)
  })
  
  rm(scenario_stack_ordered, pca_raster_scenario); gc()
} # End scenario loop

# --- 6. Save the List of PCA Raster Paths ---
pca_paths_rds_save_path <- config$pca_raster_paths_rds_path # Path from config
if (length(pca_raster_paths_list) > 0) {
  tryCatch({
    saveRDS(pca_raster_paths_list, file = pca_paths_rds_save_path)
    cat("--- List of PCA raster paths saved to:", pca_paths_rds_save_path, "---\n")
  }, error = function(e){
    warning("Failed to save the list of PCA raster paths: ", e$message, call.=FALSE)
  })
} else {
  warning("No PCA rasters were successfully generated or saved. The paths list is empty.", call.=FALSE)
}

cat("--- Script 05b finished. PCA model built and projected. ---\n")
#-------------------------------------------------------------------------------