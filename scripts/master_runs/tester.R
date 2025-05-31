

# --- Define Parameters ---
current_scenario_name <- "current"
# Determine the suffix based on config (assuming this is for env-only models)
predictor_suffix <- ifelse(config$use_pca_predictors, "_pca", "_vif")

# --- 0. Create figure directory if it doesn't exist ---
# This uses your config$base_dir to ensure it's relative to your project root
figure_output_dir <- file.path(config$base_dir, "figure_files")
if (!dir.exists(figure_output_dir)) {
  dir.create(figure_output_dir, recursive = TRUE)
  cat("Created directory for figures at:", figure_output_dir, "\n")
}


# --- 1. Load Host Sea Anemone Data ---
cat("--- Loading Host Sea Anemone Predictions (Current Scenario) ---\n")
if (!file.exists(config$anemone_species_list_file)) {
  stop("Anemone species list CSV not found: ", config$anemone_species_list_file)
}
anemone_species_df <- readr::read_csv(config$anemone_species_list_file, show_col_types = FALSE)

host_raster_files <- c()
host_short_names <- c()

for (i in 1:nrow(anemone_species_df)) {
  sp_name_sanitized <- gsub(" ", "_", anemone_species_df$scientificName[i])
  pred_file_path <- construct_prediction_filename(
    species_name_sanitized = sp_name_sanitized,
    scenario_name = current_scenario_name,
    predictor_type_suffix = predictor_suffix, 
    config = config
  )
  
  if (file.exists(pred_file_path)) {
    host_raster_files <- c(host_raster_files, pred_file_path)
    host_short_names <- c(host_short_names, sp_name_sanitized) 
  } else {
    cat("Warning: Host prediction file not found for", sp_name_sanitized, "at", pred_file_path, "\n")
  }
}

host_pred_stack <- NULL
if (length(host_raster_files) > 0) {
  host_pred_stack <- terra::rast(host_raster_files)
  names(host_pred_stack) <- host_short_names
  cat("Loaded", terra::nlyr(host_pred_stack), "host anemone prediction rasters.\n")
  
  host_richness_sum <- sum(host_pred_stack, na.rm = TRUE)
  names(host_richness_sum) <- "HostAnemoneRichness"
  cat("Calculated host anemone summed richness.\n")
  
} else {
  cat("Error: No host anemone prediction rasters found. Cannot proceed with host richness.\n")
  host_richness_sum <- NULL 
}


# --- 2. Load Anemonefish (Environmental-Only Models) Data ---
cat("\n--- Loading Anemonefish (Env-Only) Predictions (Current Scenario) ---\n")
if (!file.exists(config$anemonefish_species_list_file)) {
  stop("Anemonefish species list CSV not found: ", config$anemonefish_species_list_file)
}
anemonefish_species_df <- readr::read_csv(config$anemonefish_species_list_file, show_col_types = FALSE)

fish_raster_files <- c()
fish_short_names <- c()

for (i in 1:nrow(anemonefish_species_df)) {
  sp_name_sanitized <- gsub(" ", "_", anemonefish_species_df$scientificName[i])
  pred_file_path <- construct_prediction_filename(
    species_name_sanitized = sp_name_sanitized,
    scenario_name = current_scenario_name,
    predictor_type_suffix = predictor_suffix, 
    config = config
  )
  
  if (file.exists(pred_file_path)) {
    fish_raster_files <- c(fish_raster_files, pred_file_path)
    fish_short_names <- c(fish_short_names, sp_name_sanitized)
  } else {
    cat("Warning: Anemonefish (env-only) prediction file not found for", sp_name_sanitized, "at", pred_file_path, "\n")
  }
}

fish_pred_stack <- NULL
if (length(fish_raster_files) > 0) {
  fish_pred_stack <- terra::rast(fish_raster_files)
  names(fish_pred_stack) <- fish_short_names
  cat("Loaded", terra::nlyr(fish_pred_stack), "anemonefish (env-only) prediction rasters.\n")
  
  fish_richness_sum <- sum(fish_pred_stack, na.rm = TRUE)
  names(fish_richness_sum) <- "AnemonefishRichness_EnvOnly"
  cat("Calculated anemonefish (env-only) summed richness.\n")
  
} else {
  cat("Error: No anemonefish (env-only) prediction rasters found. Cannot proceed with fish richness.\n")
  fish_richness_sum <- NULL 
}


# --- 3. Crop to Indo-Pacific Extent (if rasters were successfully loaded) ---
cat("\n--- Cropping Rasters to Indo-Pacific Extent ---\n")
if (config$apply_indo_pacific_crop) {
  ip_extent <- terra::ext(config$indo_pacific_bbox)
  cat("Using Indo-Pacific Bounding Box for cropping:", 
      paste(config$indo_pacific_bbox, collapse=", "), "\n")
  
  if (!is.null(host_pred_stack)) {
    host_pred_stack_cropped <- tryCatch({
      terra::crop(host_pred_stack, ip_extent)
    }, error = function(e) {
      cat("Warning: Failed to crop host_pred_stack:", e$message, "\n"); host_pred_stack
    })
  } else { host_pred_stack_cropped <- NULL }
  
  if (!is.null(host_richness_sum)) {
    host_richness_sum_cropped <- tryCatch({
      terra::crop(host_richness_sum, ip_extent)
    }, error = function(e) {
      cat("Warning: Failed to crop host_richness_sum:", e$message, "\n"); host_richness_sum
    })
  } else { host_richness_sum_cropped <- NULL }
  
  if (!is.null(fish_pred_stack)) {
    fish_pred_stack_cropped <- tryCatch({
      terra::crop(fish_pred_stack, ip_extent)
    }, error = function(e) {
      cat("Warning: Failed to crop fish_pred_stack:", e$message, "\n"); fish_pred_stack
    })
  } else { fish_pred_stack_cropped <- NULL }
  
  if (!is.null(fish_richness_sum)) {
    fish_richness_sum_cropped <- tryCatch({
      terra::crop(fish_richness_sum, ip_extent)
    }, error = function(e) {
      cat("Warning: Failed to crop fish_richness_sum:", e$message, "\n"); fish_richness_sum
    })
  } else { fish_richness_sum_cropped <- NULL }
  
  cat("Cropping complete.\n")
} else {
  cat("Skipping Indo-Pacific cropping based on config.\n")
  host_pred_stack_cropped <- host_pred_stack
  host_richness_sum_cropped <- host_richness_sum
  fish_pred_stack_cropped <- fish_pred_stack
  fish_richness_sum_cropped <- fish_richness_sum
}

# --- 4. Plot Cropped Richness Maps AND SAVE THEM ---
cat("\n--- Plotting and Saving Cropped Richness Maps ---\n")

# Define plot parameters to make them look nicer for saving
# You might want to adjust col, main titles, etc.
plot_params <- list(
  col = rev(terrain.colors(255)), # Example color palette
  plg = list(loc = "right", title = "Richness", cex = 0.8), # Legend parameters
  pax = list(cex.axis = 0.8) # Axis parameters
)


if (!is.null(host_richness_sum_cropped)) {
  # Define filename for host richness plot
  host_plot_filename <- file.path(figure_output_dir, "host_richness_current_cropped.png")
  
  # Open PNG device
  png(filename = host_plot_filename, width = 800, height = 600, units = "px", res = 100)
  
  # Plot host richness (use the plot_params list)
  plot(host_richness_sum_cropped, 
       main = "Summed Host Anemone Richness (Current, Cropped)", 
       col = plot_params$col, 
       plg = plot_params$plg, 
       pax = plot_params$pax)
  
  # Close PNG device
  dev.off()
  cat("Saved host richness plot to:", host_plot_filename, "\n")
  
  # Also display in RMarkdown output if knitting
  plot(host_richness_sum_cropped, 
       main = "Summed Host Anemone Richness (Current, Cropped)", 
       col = plot_params$col, 
       plg = plot_params$plg, 
       pax = plot_params$pax)
}

if (!is.null(fish_richness_sum_cropped)) {
  # Define filename for fish richness plot
  fish_plot_filename <- file.path(figure_output_dir, "fish_richness_current_env_only_cropped.png")
  
  # Open PNG device
  png(filename = fish_plot_filename, width = 800, height = 600, units = "px", res = 100)
  
  # Plot fish richness (use the plot_params list)
  plot(fish_richness_sum_cropped, 
       main = "Summed Anemonefish Richness (Current, Env-Only, Cropped)",
       col = plot_params$col, 
       plg = plot_params$plg, 
       pax = plot_params$pax)
  
  # Close PNG device
  dev.off()
  cat("Saved anemonefish richness plot to:", fish_plot_filename, "\n")
  
  # Also display in RMarkdown output if knitting
  plot(fish_richness_sum_cropped, 
       main = "Summed Anemonefish Richness (Current, Env-Only, Cropped)",
       col = plot_params$col, 
       plg = plot_params$plg, 
       pax = plot_params$pax)
}

# Optional: Plot cropped individual stacks if needed for visual check
if (!is.null(host_pred_stack_cropped)) {
  plot(host_pred_stack_cropped)
}
if (!is.null(fish_pred_stack_cropped)) {
  plot(fish_pred_stack_cropped)
}

cat("\n--- First section of results processing finished. ---\n")
# The objects host_pred_stack_cropped, host_richness_sum_cropped, 
# fish_pred_stack_cropped, and fish_richness_sum_cropped are now available.