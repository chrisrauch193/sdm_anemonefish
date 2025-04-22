# --- Plotting Section ---
cat("\n--- Generating Comparison Plot ---\n")

# --- 1. Load Libraries ---
library(ggplot2)
library(sf)
library(terra)
library(tidyterra) # For plotting SpatRasters with ggplot
library(rnaturalearth) # For a base map
library(viridis) # For nice color scales

# --- 2. Prepare Data for Plotting ---

# Base World Map
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

# Convert points to sf objects (assuming they are data frames with lon/lat)
# Ensure correct CRS (should match the raster, likely EPSG:4326)
raster_crs <- terra::crs(tuning_predictor_stack_species) # Get CRS from the species raster

occ_sf_plot <- sf::st_as_sf(occs_coords_df, coords = c("longitude", "latitude"), crs = raster_crs)

# Check if the global background points exist before trying to plot
global_bg_exists <- exists("background_points_df") && !is.null(background_points_df) && nrow(background_points_df) > 0
if (global_bg_exists) {
  bg_sf_global_plot <- sf::st_as_sf(background_points_df, coords = c("longitude", "latitude"), crs = raster_crs)
} else {
  cat("INFO: Global background points ('background_points_df') not found or empty. Skipping from plot.\n")
}

bg_sf_species_plot <- sf::st_as_sf(background_points_df_species, coords = c("longitude", "latitude"), crs = raster_crs)

# Convert calibration vector polygon to sf object
calibration_sf_plot <- sf::st_as_sf(species_calibration_vect)

# --- 3. Define Plot Extent ---
# Get extent from the species-specific raster + a small buffer
plot_extent_vec <- terra::ext(tuning_predictor_stack_species)
xmin_plot <- as.numeric(plot_extent_vec$xmin) - 1 # Add 1-degree buffer
xmax_plot <- as.numeric(plot_extent_vec$xmax) + 1
ymin_plot <- as.numeric(plot_extent_vec$ymin) - 1
ymax_plot <- as.numeric(plot_extent_vec$ymax) + 1

# --- 4. Create the Plot ---
comparison_plot <- ggplot() +
  # Base map layer (light grey)
  geom_sf(data = world_sf, fill = "grey80", color = "white", size = 0.1) +
  
  # Species-specific predictor background (e.g., PC1)
  # Using tidyterra for easy plotting
  geom_spatraster(data = tuning_predictor_stack_species[[1]]) +
  scale_fill_viridis_c(option = "plasma", na.value = NA, name = names(tuning_predictor_stack_species)[1]) + # Use viridis colors
  
  # Calibration Area Polygon Outline (e.g., blue outline)
  geom_sf(data = calibration_sf_plot, fill = NA, color = "blue", linewidth = 0.8) +
  
  # Plot Species-Specific Background Points (e.g., small black dots)
  geom_sf(data = bg_sf_species_plot, color = "black", size = 0.2, shape = ".", alpha = 0.5) +
  
  # Plot Global Background Points IF THEY EXIST (e.g., small grey dots)
  {
    if (global_bg_exists)
      geom_sf(data = bg_sf_global_plot, color = "grey50", size = 0.1, shape = ".", alpha = 0.3)
  } +
  
  # Plot Occurrence Points (e.g., larger red points)
  geom_sf(data = occ_sf_plot, color = "red", size = 1.5, shape = 16, alpha = 0.8) +
  
  # Set Coordinate System and Limits
  coord_sf(crs = raster_crs, # Use the raster's CRS
           xlim = c(xmin_plot, xmax_plot),
           ylim = c(ymin_plot, ymax_plot),
           expand = FALSE) + # Prevent ggplot from adding extra space
  
  # Labels and Theme
  labs(
    title = paste("SDM Extent Comparison:", species_name),
    subtitle = "Occurrences (Red), Species-Specific BG (Black), Global BG (Grey, if shown)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Light blue ocean
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "right"
  )

# --- 5. Print the Plot ---
print(comparison_plot)

# --- 6. Save the Plot (Optional) ---
plot_save_filename <- file.path(config$log_dir_base, paste0("bg_comparison_", species_name_sanitized, predictor_type_suffix, ".png"))
tryCatch({
  ggsave(plot_save_filename, plot = comparison_plot, width = 10, height = 7, dpi = 300, bg = "white")
  cat("INFO: Comparison plot saved to:", plot_save_filename, "\n")
}, error = function(e) {
  cat("WARN: Failed to save comparison plot:", e$message, "\n")
})

cat("--- Plotting Section Finished ---\n")