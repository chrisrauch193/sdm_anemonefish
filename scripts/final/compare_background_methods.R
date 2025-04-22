# scripts/final/compare_background_methods_sdmtune.R
#-------------------------------------------------------------------------------
# Script to compare background sampling strategies for ONE species using SDMtune.
# Methods Compared:
# 1. Global Background Sampling (potentially masked by coral, Indo-Pacific extent)
# 2. Species-Specific Alpha Hull Extent Background Sampling
# 3. Species-Specific OBIS/MPAEU Ecoregion/Depth Extent Background Sampling
#-------------------------------------------------------------------------------
rm(list=ls()); gc() # Clean workspace

cat("--- Running SDMtune Background Comparison Script ---\n")

# --- 1. Setup ---
cat("--- Loading Config and Packages ---\n")
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
config_path_options <- c(file.path(script_dir, "..", "config.R"), file.path("scripts", "config.R"), "config.R")
config_path <- NULL; for(path_opt in config_path_options) if(file.exists(path_opt)) { config_path <- path_opt; break }
if (is.null(config_path)) stop("FATAL: config.R not found.")
source(config_path); if (!exists("config")) stop("FATAL: config list missing.")

# Load required packages
pacman::p_load(terra, sf, dplyr, readr, tools, stringr, SDMtune, ggplot2, tidyterra,
               rnaturalearth, viridis, concaveman, # Added concaveman
               log4r, future, furrr, progressr) # Added logging/parallel just in case helpers use them

cat("--- Sourcing Helper Functions ---\n")
helper_paths <- c(
  file.path(config$helpers_dir, "env_processing_helpers.R"),
  file.path(config$helpers_dir, "sdm_modeling_helpers.R") # Ensure generate_sdm_background_species_extent & generate_sdm_background_obis_extent are here
)
missing_helpers <- helper_paths[!file.exists(helper_paths)]
if(length(missing_helpers) > 0) stop("Missing helper(s): ", paste(missing_helpers, collapse=", "))
invisible(sapply(helper_paths, source))

# --- 2. Define Target Species & Scenario ---
group_name <- "anemone" # Or "anemonefish" - CHOOSE THE GROUP
species_list_file <- config$anemone_species_list_file # Adjust if needed
occurrence_dir <- config$anemone_occurrence_dir # Adjust if needed

use_pca <- config$use_pca_predictors # Ensure this matches the predictors you want to test
tuning_scenario <- "current" # Scenario to use for tuning/comparison

cat("--- Target Group:", group_name, "---\n")
cat("--- Using Predictors:", ifelse(use_pca, "PCA", "VIF"), "---\n")
cat("--- Tuning Scenario:", tuning_scenario, "---\n")

species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
species_row <- species_df[1, ] # <<< --- SELECT THE TEST SPECIES ROW --- >>>
species_name <- species_row$scientificName
species_name_sanitized <- gsub(" ", "_", species_name)
species_aphia_id <- species_row$AphiaID
cat("--- Target Species:", species_name, "(AphiaID:", species_aphia_id, ") ---\n")

# --- 3. Load Global Predictor Stack ---
cat("--- Loading GLOBAL Predictor Stack for Scenario:", tuning_scenario, "---\n")
global_predictor_stack <- NULL
if (use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path; if (!file.exists(pca_paths_rds)) stop("Global PCA paths RDS file not found.")
  predictor_paths_list <- readRDS(pca_paths_rds)
  global_predictor_path <- predictor_paths_list[[tuning_scenario]]
  if (is.null(global_predictor_path) || !file.exists(global_predictor_path)) stop("GLOBAL PCA stack path missing for tuning.")
  global_predictor_stack <- tryCatch(terra::rast(global_predictor_path), error = function(e) stop("Failed load GLOBAL PCA stack: ", e$message))
} else {
  # Load VIF selected variables if not using PCA
  selected_vars_vif <- if(group_name == "anemone") config$final_vars_vif_anemone else config$final_vars_vif_anemonefish_env
  if(is.null(selected_vars_vif)) stop("VIF variable list not defined in config for ", group_name)
  scenario_vif_vars <- generate_scenario_variable_list(selected_vars_vif, tuning_scenario, config)
  global_predictor_stack <- load_selected_env_data(tuning_scenario, scenario_vif_vars, config) # Load only selected VIF vars
}
if (is.null(global_predictor_stack)) stop("Failed to load GLOBAL predictor stack.")
if (terra::crs(global_predictor_stack) == "") stop("GLOBAL predictor stack missing CRS.")
cat("  GLOBAL stack loaded. Layers:", paste(names(global_predictor_stack), collapse=", "), "\n")

# --- 4. Load & Prep Occurrences ---
cat("--- Loading and Cleaning Occurrences ---\n")
config_for_occ_load <- config; config_for_occ_load$predictor_stack_for_thinning <- global_predictor_stack
occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL)
if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) stop("Insufficient occurrences.")
occs_coords_df <- as.data.frame(occ_data_result$coords); colnames(occs_coords_df) <- c("longitude", "latitude")
occs_sf_clean <- sf::st_as_sf(occs_coords_df, coords=c("longitude","latitude"), crs=config$occurrence_crs) # CRS from config
# Ensure occs_sf_clean has the correct CRS matching the predictor stack
if(sf::st_crs(occs_sf_clean) != sf::st_crs(terra::crs(global_predictor_stack))){
  occs_sf_clean <- sf::st_transform(occs_sf_clean, crs=sf::st_crs(terra::crs(global_predictor_stack)))
  cat("  Transformed occurrence CRS to match predictors.\n")
}
cat("  Occurrence count:", nrow(occs_coords_df), "\n")

# --- 5. Method 1: Global Background & Tuning ---
cat("\n--- Running Method 1: GLOBAL Background Sampling & SDMtune ---\n")
bg_global_df <- generate_sdm_background(global_predictor_stack, config$background_points_n, config, logger = NULL, seed = 123)
if (is.null(bg_global_df)) stop("Failed global background generation.")
colnames(bg_global_df) <- c("longitude", "latitude") # Ensure consistent naming for SDMtune
cat("  Generated", nrow(bg_global_df), "global background points.\n")

cat("  Running SDMtune (Global Background)...\n")
tuning_output_global <- run_sdm_tuning_scv(occs_coords_df, global_predictor_stack, bg_global_df, config, logger=NULL, species_name=species_name)
if (is.null(tuning_output_global)) warning("SDMtune failed for Global Background.")

# --- 6. Method 2: Species-Specific Alpha Hull Background & Tuning ---
cat("\n--- Running Method 2: SPECIES-SPECIFIC Alpha Hull Background Sampling & SDMtune ---\n")
bg_species_result_hull <- generate_sdm_background_species_extent( # ASSUMES THIS HELPER IS DEFINED
  occs_sf = occs_sf_clean, # Pass the sf object with correct CRS
  global_predictor_stack = global_predictor_stack,
  config = config,
  logger = NULL,
  seed_offset = 1 # Use different offset for seed
)

tuning_output_species_hull <- NULL # Initialize
predictor_stack_species_hull <- NULL
species_calibration_vect_hull <- NULL
bg_species_df_hull <- NULL

if(is.null(bg_species_result_hull)) {
  cat("WARN: Failed species-specific (Alpha Hull) background generation helper.\n")
} else {
  bg_species_df_hull <- bg_species_result_hull$background_points
  predictor_stack_species_hull <- bg_species_result_hull$species_specific_stack
  species_calibration_vect_hull <- bg_species_result_hull$calibration_polygon # This is the alpha hull SpatVector
  cat("  Generated", nrow(bg_species_df_hull), "species-specific (Alpha Hull) background points.\n")
  cat("  Alpha Hull specific stack created.\n")
  
  cat("  Running SDMtune (Species Alpha Hull Background)...\n")
  tuning_output_species_hull <- run_sdm_tuning_scv(
    occs_coords_df, # Original coordinates DF
    predictor_stack_species_hull, # The HULL-masked stack
    bg_species_df_hull, # The HULL background points
    config,
    logger=NULL,
    species_name=species_name
  )
  if (is.null(tuning_output_species_hull)) warning("SDMtune failed for Species Alpha Hull Background.")
}

# --- 7. Method 3: OBIS Ecoregion/Depth Background & Tuning ---
cat("\n--- Running Method 3: OBIS Ecoregion/Depth Background Sampling & SDMtune ---\n")
bg_obis_result <- generate_sdm_background_obis_extent( # ASSUMES THIS HELPER IS DEFINED
  occs_sf = occs_sf_clean, # Use the same cleaned/projected occurrences
  global_predictor_stack = global_predictor_stack, # Pass the global stack
  config = config, # Pass the config list
  logger = NULL, # Pass logger if you have one initialized
  species_log_file = NULL, # Pass species log file if needed
  seed_offset = 2 # Use different offset
)

tuning_output_obis <- NULL # Initialize
predictor_stack_obis <- NULL
obis_calibration_vect <- NULL
bg_obis_df <- NULL

if(is.null(bg_obis_result)) {
  cat("WARN: Failed OBIS-specific background generation helper.\n")
} else {
  bg_obis_df <- bg_obis_result$background_points
  predictor_stack_obis <- bg_obis_result$species_specific_stack
  obis_calibration_vect <- bg_obis_result$calibration_polygon # This is the ecoregion SpatVector
  cat("  Generated", nrow(bg_obis_df), "OBIS-extent background points.\n")
  cat("  OBIS-specific stack created.\n")
  
  cat("  Running SDMtune (OBIS Background)...\n")
  tuning_output_obis <- run_sdm_tuning_scv(
    occs_coords_df, # Original coordinates DF
    predictor_stack_obis, # The OBIS-masked stack
    bg_obis_df, # The OBIS background points
    config,
    logger=NULL,
    species_name=species_name
  )
  if (is.null(tuning_output_obis)) warning("SDMtune failed for OBIS Background.")
}


# --- 8. Compare Results ---
cat("\n--- Comparing SDMtune Results ---\n")
comparison_list <- list()

# Global
if (!is.null(tuning_output_global) && inherits(tuning_output_global, "SDMtune")) {
  res_g <- tuning_output_global@results
  best_g_idx <- if(tolower(config$sdm_evaluation_metric) == "aicc") which.min(res_g$AICc) else which.max(res_g[[paste0("test_", toupper(config$sdm_evaluation_metric))]])
  if(length(best_g_idx) == 0 || is.na(best_g_idx)) best_g_idx <- 1 # Fallback
  comparison_list$GlobalBG <- res_g[best_g_idx, ]
} else { cat("WARN: Global tuning results missing or invalid.\n") }

# Species Hull
if (!is.null(tuning_output_species_hull) && inherits(tuning_output_species_hull, "SDMtune")) {
  res_s_hull <- tuning_output_species_hull@results
  best_s_hull_idx <- if(tolower(config$sdm_evaluation_metric) == "aicc") which.min(res_s_hull$AICc) else which.max(res_s_hull[[paste0("test_", toupper(config$sdm_evaluation_metric))]])
  if(length(best_s_hull_idx) == 0 || is.na(best_s_hull_idx)) best_s_hull_idx <- 1 # Fallback
  comparison_list$SpeciesHullBG <- res_s_hull[best_s_hull_idx, ]
} else { cat("WARN: Species-specific (Alpha Hull) tuning results missing or invalid.\n") }

# OBIS Extent
if (!is.null(tuning_output_obis) && inherits(tuning_output_obis, "SDMtune")) {
  res_o <- tuning_output_obis@results
  best_o_idx <- if(tolower(config$sdm_evaluation_metric) == "aicc") which.min(res_o$AICc) else which.max(res_o[[paste0("test_", toupper(config$sdm_evaluation_metric))]])
  if(length(best_o_idx) == 0 || is.na(best_o_idx)) best_o_idx <- 1 # Fallback
  comparison_list$OBISBG <- res_o[best_o_idx, ]
} else { cat("WARN: OBIS-specific tuning results missing or invalid.\n") }

# Combine and Print Comparison
if (length(comparison_list) > 0) {
  comparison_df <- dplyr::bind_rows(comparison_list, .id = "Method")
  metrics_to_show <- intersect(c("Method", "reg", "fc", "test_AUC", "test_TSS", "AICc"), names(comparison_df))
  cat("Comparison of Best Models (selected by", config$sdm_evaluation_metric, "):\n")
  print(comparison_df[, metrics_to_show], row.names = FALSE, digits = 4)
} else { cat("Could not generate comparison table.\n") }


# --- 9. Visualization ---
cat("\n--- Generating Comparison Plot ---\n")

# --- 9a. Prepare Data for Plotting ---
world_sf <- tryCatch(ne_countries(scale = "medium", returnclass = "sf"), error = function(e) { cat("Failed load rnaturalearth map."); NULL })
plot_crs_terra <- terra::crs(global_predictor_stack) # Use stack's CRS
plot_crs_sf <- sf::st_crs(plot_crs_terra) # Corresponding sf CRS

# Points
occ_sf_plot <- sf::st_as_sf(occs_coords_df, coords = c("longitude", "latitude"), crs = plot_crs_sf)
bg_sf_global_plot <- sf::st_as_sf(bg_global_df, coords = c("longitude", "latitude"), crs = plot_crs_sf) # Adjusted col names
bg_sf_species_hull_plot <- if(!is.null(bg_species_df_hull)) sf::st_as_sf(bg_species_df_hull, coords = c("longitude", "latitude"), crs = plot_crs_sf) else NULL # Adjusted col names
bg_sf_obis_plot <- if(!is.null(bg_obis_df)) sf::st_as_sf(bg_obis_df, coords = c("longitude", "latitude"), crs = plot_crs_sf) else NULL # Adjusted col names

# Polygons
alpha_hull_sf_plot <- if (!is.null(species_calibration_vect_hull)) sf::st_as_sf(species_calibration_vect_hull) else NULL
obis_calibration_sf_plot <- if (!is.null(obis_calibration_vect)) sf::st_as_sf(obis_calibration_vect) else NULL

# Ensure polygons have the correct CRS for plotting
if(!is.null(alpha_hull_sf_plot) && sf::st_crs(alpha_hull_sf_plot) != plot_crs_sf) alpha_hull_sf_plot <- sf::st_transform(alpha_hull_sf_plot, crs=plot_crs_sf)
if(!is.null(obis_calibration_sf_plot) && sf::st_crs(obis_calibration_sf_plot) != plot_crs_sf) obis_calibration_sf_plot <- sf::st_transform(obis_calibration_sf_plot, crs=plot_crs_sf)


# --- 9b. Define Plot Extent ---
# Create a list of valid polygons for extent calculation
extent_polys <- list(alpha_hull_sf_plot, obis_calibration_sf_plot)
extent_polys <- extent_polys[!sapply(extent_polys, is.null)] # Remove NULLs

if(length(extent_polys) > 0) {
  # Calculate combined bounding box of valid polygons
  combined_bbox <- sf::st_bbox(do.call(rbind, extent_polys))
} else {
  # Fallback to occurrences if no valid polygons exist
  combined_bbox <- sf::st_bbox(occ_sf_plot)
}

# Add a buffer (e.g., 5 degrees) to the combined bounding box
plot_bbox_buffered <- sf::st_bbox(sf::st_buffer(sf::st_as_sfc(combined_bbox), dist = 5))

xmin_plot <- plot_bbox_buffered$xmin; xmax_plot <- plot_bbox_buffered$xmax
ymin_plot <- plot_bbox_buffered$ymin; ymax_plot <- plot_bbox_buffered$ymax


# --- 9c. Create Plot ---
plot_raster <- global_predictor_stack[[1]] # Use first layer of GLOBAL stack for background color

comparison_plot <- ggplot()
if(!is.null(world_sf)) { comparison_plot <- comparison_plot + geom_sf(data = world_sf, fill = "grey80", color = "white", size = 0.1) }
comparison_plot <- comparison_plot +
  # Background raster (using the global one cropped implicitly by coord_sf)
  geom_spatraster(data = plot_raster) +
  scale_fill_viridis_c(option = "plasma", na.value = NA, name = names(plot_raster)[1]) +
  # Polygons
  { if (!is.null(alpha_hull_sf_plot)) geom_sf(data = alpha_hull_sf_plot, fill = NA, color = "blue", linewidth = 0.8, aes(linetype = "Alpha Hull Extent")) } +
  { if (!is.null(obis_calibration_sf_plot)) geom_sf(data = obis_calibration_sf_plot, fill = NA, color = "darkgreen", linewidth = 0.8, aes(linetype = "OBIS Extent")) } +
  # Background Points
  geom_sf(data = bg_sf_global_plot, aes(shape = "Global BG", size = "Global BG", alpha = "Global BG", color = "Global BG")) +
  { if (!is.null(bg_sf_species_hull_plot)) geom_sf(data = bg_sf_species_hull_plot, aes(shape = "Alpha Hull BG", size = "Alpha Hull BG", alpha = "Alpha Hull BG", color = "Alpha Hull BG")) } +
  { if (!is.null(bg_sf_obis_plot)) geom_sf(data = bg_sf_obis_plot, aes(shape = "OBIS BG", size = "OBIS BG", alpha = "OBIS BG", color = "OBIS BG")) } +
  # Occurrence Points
  geom_sf(data = occ_sf_plot, aes(shape = "Occurrence", size = "Occurrence", alpha = "Occurrence", color = "Occurrence")) +
  # Coordinate System and Limits using calculated buffered bbox
  coord_sf(crs = plot_crs_sf, # Use the defined plot CRS
           xlim = c(xmin_plot, xmax_plot),
           ylim = c(ymin_plot, ymax_plot),
           expand = FALSE) +
  # Manual Scales
  scale_shape_manual(name = "Data Type", values = c("Global BG" = 1, "Alpha Hull BG" = 16, "OBIS BG" = 16, "Occurrence" = 17)) +
  scale_size_manual(name = "Data Type", values = c("Global BG" = 0.8, "Alpha Hull BG" = 0.5, "OBIS BG" = 0.5, "Occurrence" = 2.0)) +
  scale_alpha_manual(name = "Data Type", values = c("Global BG" = 0.4, "Alpha Hull BG" = 0.6, "OBIS BG" = 0.6, "Occurrence" = 0.8)) +
  scale_color_manual(name = "Data Type", values = c("Global BG" = "grey50", "Alpha Hull BG" = "black", "OBIS BG" = "darkorange", "Occurrence" = "red")) + # Changed OBIS BG color
  scale_linetype_manual(name = "Extent", values = c("Alpha Hull Extent" = "solid", "OBIS Extent" = "dashed")) +
  guides(
    shape = guide_legend(title = "Data Type", order = 1, override.aes = list(alpha = 0.8)),
    color = guide_legend(title = "Data Type", order = 1, override.aes = list(alpha = 0.8)),
    size = "none", alpha = "none",
    linetype = guide_legend(title = "Extent", order = 2)
  ) +
  labs(
    title = paste("Background Sampling Comparison:", species_name),
    subtitle = "Occs(Red), GlobalBG(Grey), HullBG(Black), OBISBG(Orange). Polygons show extents.",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "right",
    legend.box = "vertical"
  )

# --- 9d. Print and Save Plot ---
print(comparison_plot)
plot_save_dir <- "plots" # Save in a 'plots' subdirectory relative to script
dir.create(plot_save_dir, showWarnings = FALSE)
plot_save_filename <- file.path(plot_save_dir, paste0("bg_comparison_all_", species_name_sanitized, "_sdmtune.png"))
ggsave(plot_save_filename, plot = comparison_plot, width = 11, height = 7, dpi = 300, bg = "white")
cat(paste("\nComparison plot saved to:", plot_save_filename, "\n"))

cat("\n--- compare_background_methods_sdmtune.R finished. ---\n")
#-------------------------------------------------------------------------------