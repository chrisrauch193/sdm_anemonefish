# # scripts/final/compare_background_methods_sdmtune.R
# #-------------------------------------------------------------------------------
# # Script to compare background sampling strategies for ONE species using SDMtune.
# # Methods Compared:
# # 1. Global Background Sampling (potentially masked by coral, Indo-Pacific extent)
# # 2. Species-Specific Alpha Hull Extent Background Sampling
# # 4. Species-Specific OBIS Ecoregion/Depth Extent Background Sampling (Exact Replication Attempt)
# # Generates SEPARATE plots for each method with CONSISTENT point styling.
# # Version: v4 - Consistent Plot Styling.
# #-------------------------------------------------------------------------------
# rm(list=ls()); gc() # Clean workspace
# 
# cat("--- Running SDMtune Background Comparison Script (Separate Plots v4 - Consistent Styling) ---\n")
# 
# # --- 1. Setup ---
# cat("--- Loading Config and Packages ---\n")
# script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
# # Adjust path finding for config.R based on your execution context
# config_path_options <- c(
#   file.path(script_dir, "config.R"),              # If run from scripts/final
#   file.path(script_dir, "..", "config.R"),          # If run from scripts/
#   file.path(script_dir, "..", "..", "config.R"),     # If run from base directory
#   file.path(script_dir, "..", "..", "scripts", "config.R") # Common structure
# )
# config_path <- NULL
# for(path_opt in config_path_options) {
#   if(file.exists(path_opt)) {
#     config_path <- path_opt
#     cat("Found config at:", config_path, "\n")
#     break
#   }
# }
# if (is.null(config_path)) stop("FATAL: config.R not found in expected locations.")
# source(config_path); if (!exists("config")) stop("FATAL: config list missing after sourcing.")
# 
# pacman::p_load(terra, sf, dplyr, readr, tools, stringr, SDMtune, ggplot2, tidyterra,
#                rnaturalearth, viridis, concaveman, # Added concaveman
#                log4r, future, furrr, progressr) # Added logging/parallel just in case helpers use them
# 
# cat("--- Sourcing Helper Functions ---\n")
# helper_paths <- c(
#   file.path(config$helpers_dir, "env_processing_helpers.R"),
#   file.path(config$helpers_dir, "sdm_modeling_helpers.R") # Ensure generate_sdm_background_..., run_sdm_tuning_scv etc. are here
# )
# missing_helpers <- helper_paths[!file.exists(helper_paths)]
# if(length(missing_helpers) > 0) stop("Missing helper(s): ", paste(missing_helpers, collapse=", "))
# invisible(sapply(helper_paths, source))
# cat("Helpers sourced.\n")
# 
# # --- 2. Define Target Species & Scenario ---
# group_name <- "anemone" # Or "anemonefish" - CHOOSE THE GROUP
# species_list_file <- config$anemone_species_list_file # Adjust if needed
# occurrence_dir <- config$anemone_occurrence_dir # Adjust if needed
# 
# use_pca <- config$use_pca_predictors # Ensure this matches the predictors you want to test
# tuning_scenario <- "current" # Scenario to use for tuning/comparison
# 
# cat("--- Target Group:", group_name, "---\n")
# cat("--- Using Predictors:", ifelse(use_pca, "PCA", "VIF"), "---\n")
# cat("--- Tuning Scenario:", tuning_scenario, "---\n")
# 
# species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
# species_row <- species_df[1, ] # <<< --- SELECT THE TEST SPECIES ROW --- >>>
# species_name <- species_row$scientificName
# species_name_sanitized <- gsub(" ", "_", species_name)
# species_aphia_id <- species_row$AphiaID
# cat("--- Target Species:", species_name, "(AphiaID:", species_aphia_id, ") ---\n")
# 
# # --- 3. Load Global Predictor Stack ---
# cat("--- Loading GLOBAL Predictor Stack for Scenario:", tuning_scenario, "---\n")
# global_predictor_stack <- NULL
# if (use_pca) {
#   pca_paths_rds <- config$pca_raster_paths_rds_path; if (!file.exists(pca_paths_rds)) stop("Global PCA paths RDS file not found.")
#   predictor_paths_list <- readRDS(pca_paths_rds)
#   global_predictor_path <- predictor_paths_list[[tuning_scenario]]
#   if (is.null(global_predictor_path) || !file.exists(global_predictor_path)) stop("GLOBAL PCA stack path missing for tuning.")
#   global_predictor_stack <- tryCatch(terra::rast(global_predictor_path), error = function(e) stop("Failed load GLOBAL PCA stack: ", e$message))
# } else {
#   selected_vars_vif <- if(group_name == "anemone") config$final_vars_vif_anemone else config$final_vars_vif_anemonefish_env
#   if(is.null(selected_vars_vif)) stop("VIF variable list not defined in config for ", group_name)
#   scenario_vif_vars <- generate_scenario_variable_list(selected_vars_vif, tuning_scenario, config)
#   global_predictor_stack <- load_selected_env_data(tuning_scenario, scenario_vif_vars, config)
# }
# if (is.null(global_predictor_stack)) stop("Failed to load GLOBAL predictor stack.")
# if (terra::crs(global_predictor_stack) == "") stop("GLOBAL predictor stack missing CRS.")
# cat("  GLOBAL stack loaded. Layers:", paste(names(global_predictor_stack), collapse=", "), "\n")
# 
# # --- 4. Load & Prep Occurrences ---
# cat("--- Loading and Cleaning Occurrences ---\n")
# config_for_occ_load <- config; config_for_occ_load$predictor_stack_for_thinning <- global_predictor_stack
# occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL)
# if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) stop("Insufficient occurrences.")
# occs_coords_df <- as.data.frame(occ_data_result$coords); colnames(occs_coords_df) <- c("longitude", "latitude")
# occs_sf_clean <- sf::st_as_sf(occs_coords_df, coords=c("longitude","latitude"), crs=config$occurrence_crs) # CRS from config
# # Ensure occs_sf_clean has the correct CRS matching the predictor stack
# if(sf::st_crs(occs_sf_clean) != sf::st_crs(terra::crs(global_predictor_stack))){
#   occs_sf_clean <- sf::st_transform(occs_sf_clean, crs=sf::st_crs(terra::crs(global_predictor_stack)))
#   cat("  Transformed occurrence CRS to match predictors.\n")
# }
# cat("  Occurrence count:", nrow(occs_coords_df), "\n")
# 
# 
# # --- 5. Method 1: Global Background & Tuning ---
# cat("\n--- Running Method 1: GLOBAL Background Sampling & SDMtune ---\n")
# bg_global_df <- generate_sdm_background(global_predictor_stack, config$background_points_n, config, logger = NULL, seed = 123)
# if (is.null(bg_global_df)) stop("Failed global background generation.")
# colnames(bg_global_df) <- c("longitude", "latitude")
# cat("  Generated", nrow(bg_global_df), "global background points.\n")
# 
# cat("  Running SDMtune (Global Background)...\n")
# tuning_output_global <- run_sdm_tuning_scv(occs_coords_df, global_predictor_stack, bg_global_df, config, logger=NULL, species_name=species_name)
# if (is.null(tuning_output_global)) warning("SDMtune failed for Global Background.")
# 
# 
# # --- 6. Method 2: Species-Specific Alpha Hull Background & Tuning ---
# cat("\n--- Running Method 2: SPECIES-SPECIFIC Alpha Hull Background Sampling & SDMtune ---\n")
# bg_species_result_hull <- generate_sdm_background_species_extent( # Make sure this helper is defined
#   occs_sf = occs_sf_clean,
#   global_predictor_stack = global_predictor_stack,
#   config = config,
#   logger = NULL,
#   seed_offset = 1
# )
# 
# tuning_output_species_hull <- NULL # Initialize
# predictor_stack_species_hull <- NULL
# species_calibration_vect_hull <- NULL
# bg_species_df_hull <- NULL
# 
# if(is.null(bg_species_result_hull)) {
#   cat("WARN: Failed species-specific (Alpha Hull) background generation helper.\n")
# } else {
#   bg_species_df_hull <- bg_species_result_hull$background_points
#   predictor_stack_species_hull <- bg_species_result_hull$species_specific_stack
#   species_calibration_vect_hull <- bg_species_result_hull$calibration_polygon
#   colnames(bg_species_df_hull) <- c("longitude", "latitude") # Ensure consistent naming
#   cat("  Generated", nrow(bg_species_df_hull), "species-specific (Alpha Hull) background points.\n")
#   cat("  Alpha Hull specific stack created.\n")
#   
#   cat("  Running SDMtune (Species Alpha Hull Background)...\n")
#   tuning_output_species_hull <- run_sdm_tuning_scv(
#     occs_coords_df,
#     predictor_stack_species_hull,
#     bg_species_df_hull,
#     config,
#     logger=NULL,
#     species_name=species_name
#   )
#   if (is.null(tuning_output_species_hull)) warning("SDMtune failed for Species Alpha Hull Background.")
# }

# --- 7. Method 4: OBIS Ecoregion/Depth Background (Exact Replication Attempt) ---
# (Renumbered from Method 3 in previous response as Method 3 was commented out)
cat("\n--- Running Method 4: OBIS Ecoregion/Depth (Exact) Background Sampling & SDMtune ---\n")
bg_obis_exact_result <- generate_sdm_background_obis_part2( # Ensure this helper exists
  occs_sf = occs_sf_clean,
  global_predictor_stack = global_predictor_stack,
  config = config,
  seed_offset = 3
)

tuning_output_obis_exact <- NULL # Initialize
predictor_stack_obis_exact <- NULL
obis_calibration_vect_exact <- NULL # This will be the polygon *before* depth masking
bg_obis_exact_df <- NULL

if(is.null(bg_obis_exact_result)) {
  cat("WARN: Failed OBIS-Exact background generation helper.\n")
} else {
  bg_obis_exact_df <- bg_obis_exact_result$background_points
  predictor_stack_obis_exact <- bg_obis_exact_result$species_specific_stack # Stack after masking
  obis_calibration_vect_exact <- bg_obis_exact_result$calibration_polygon # Polygon before masking
  colnames(bg_obis_exact_df) <- c("longitude", "latitude") # Ensure consistent naming
  cat("  Generated", nrow(bg_obis_exact_df), "OBIS-Exact extent background points.\n")
  cat("  OBIS-Exact specific stack created.\n")
  
  cat("  Running SDMtune (OBIS-Exact Background)...\n")
  tuning_output_obis_exact <- run_sdm_tuning_scv(
    occs_coords_df,
    predictor_stack_obis_exact, # Use the OBIS-Exact masked stack
    bg_obis_exact_df,
    config,
    logger=NULL,
    species_name=species_name
  )
  if (is.null(tuning_output_obis_exact)) warning("SDMtune failed for OBIS-Exact Background.")
}


# --- 8. Compare Results ---
# (Keep this section exactly as in the previous version, including saving the comparison table) ...
cat("\n--- Comparing SDMtune Results ---\n")
comparison_list <- list()
if (!is.null(tuning_output_global) && inherits(tuning_output_global, "SDMtune")) { res_g <- tuning_output_global@results; metric_col <- paste0("test_", toupper(config$sdm_evaluation_metric)); best_g_idx <- if(tolower(config$sdm_evaluation_metric) == "aicc") which.min(res_g$AICc) else which.max(res_g[[metric_col]]); if(length(best_g_idx) == 0 || is.na(best_g_idx)) best_g_idx <- 1; comparison_list$GlobalBG <- res_g[best_g_idx, ] } else { cat("WARN: Global tuning results missing or invalid.\n") }
if (!is.null(tuning_output_species_hull) && inherits(tuning_output_species_hull, "SDMtune")) { res_s_hull <- tuning_output_species_hull@results; metric_col <- paste0("test_", toupper(config$sdm_evaluation_metric)); best_s_hull_idx <- if(tolower(config$sdm_evaluation_metric) == "aicc") which.min(res_s_hull$AICc) else which.max(res_s_hull[[metric_col]]); if(length(best_s_hull_idx) == 0 || is.na(best_s_hull_idx)) best_s_hull_idx <- 1; comparison_list$SpeciesHullBG <- res_s_hull[best_s_hull_idx, ] } else { cat("WARN: Species-specific (Alpha Hull) tuning results missing or invalid.\n") }
if (!is.null(tuning_output_obis_exact) && inherits(tuning_output_obis_exact, "SDMtune")) { res_o_exact <- tuning_output_obis_exact@results; metric_col <- paste0("test_", toupper(config$sdm_evaluation_metric)); best_o_exact_idx <- if(tolower(config$sdm_evaluation_metric) == "aicc") which.min(res_o_exact$AICc) else which.max(res_o_exact[[metric_col]]); if(length(best_o_exact_idx) == 0 || is.na(best_o_exact_idx)) best_o_exact_idx <- 1; comparison_list$OBISExactBG <- res_o_exact[best_o_exact_idx, ] } else { cat("WARN: OBIS-Exact tuning results missing or invalid.\n") }
if (length(comparison_list) > 0) { comparison_df <- dplyr::bind_rows(comparison_list, .id = "Method"); metrics_to_show <- intersect(c("Method", "reg", "fc", "test_AUC", "test_TSS", "AICc"), names(comparison_df)); cat("Comparison of Best Models (selected by", config$sdm_evaluation_metric, "):\n"); print(comparison_df[, metrics_to_show], row.names = FALSE, digits = 4); comp_save_dir <- "plots"; dir.create(comp_save_dir, showWarnings = FALSE); comp_save_filename <- file.path(comp_save_dir, paste0("bg_comparison_table_", species_name_sanitized, "_sdmtune.csv")); tryCatch(write.csv(comparison_df, comp_save_filename, row.names = FALSE), error=function(e){cat("WARN: Failed to save comparison table:", e$message, "\n")}); cat(paste("Comparison table saved to:", comp_save_filename, "\n")) } else { cat("Could not generate comparison table.\n") }


# --- 9. Visualization ---
cat("\n--- Generating Separate Comparison Plots ---\n")

# --- 9a. Common Plotting Elements ---
world_sf <- tryCatch(ne_countries(scale = "medium", returnclass = "sf"), error = function(e) { cat("Failed load rnaturalearth map."); NULL })
plot_crs_terra <- terra::crs(global_predictor_stack)
plot_crs_sf <- sf::st_crs(plot_crs_terra)
occ_sf_plot <- sf::st_as_sf(occs_coords_df, coords = c("longitude", "latitude"), crs = plot_crs_sf)
plot_raster_global <- global_predictor_stack[[1]] # Use first layer for background color

# --- 9b. Define Overall Plot Extent ---
# (Keep extent calculation based on polygons and points as before)
all_geoms_for_extent_list <- list(occ_sf_plot)
if(!is.null(species_calibration_vect_hull)) all_geoms_for_extent_list[[length(all_geoms_for_extent_list)+1]] <- sf::st_as_sf(species_calibration_vect_hull)
if(!is.null(obis_calibration_vect_exact)) all_geoms_for_extent_list[[length(all_geoms_for_extent_list)+1]] <- sf::st_as_sf(obis_calibration_vect_exact)
if(!is.null(bg_global_df)) all_geoms_for_extent_list[[length(all_geoms_for_extent_list)+1]] <- sf::st_as_sf(bg_global_df, coords = c("longitude", "latitude"), crs = plot_crs_sf)
if(!is.null(bg_species_df_hull)) all_geoms_for_extent_list[[length(all_geoms_for_extent_list)+1]] <- sf::st_as_sf(bg_species_df_hull, coords = c("longitude", "latitude"), crs = plot_crs_sf)
if(!is.null(bg_obis_exact_df)) all_geoms_for_extent_list[[length(all_geoms_for_extent_list)+1]] <- sf::st_as_sf(bg_obis_exact_df, coords = c("longitude", "latitude"), crs = plot_crs_sf)
all_geoms_for_extent_list <- all_geoms_for_extent_list[!sapply(all_geoms_for_extent_list, is.null)]
if(length(all_geoms_for_extent_list) > 0) { geometries_list <- lapply(all_geoms_for_extent_list, sf::st_geometry); combined_geometries <- do.call(c, geometries_list); combined_bbox_all <- sf::st_bbox(combined_geometries) } else { combined_bbox_all <- sf::st_bbox(occ_sf_plot) }
if(all(is.finite(combined_bbox_all))) { plot_bbox_buffered_all <- sf::st_bbox(sf::st_buffer(sf::st_as_sfc(combined_bbox_all), dist=5)) } else { cat("WARN: Combined bounding box is invalid. Using occurrence bbox + buffer for plot extent.\n"); plot_bbox_buffered_all <- sf::st_bbox(sf::st_buffer(occ_sf_plot, dist=5)) }
xmin_plot <- plot_bbox_buffered_all$xmin; xmax_plot <- plot_bbox_buffered_all$xmax
ymin_plot <- plot_bbox_buffered_all$ymin; ymax_plot <- plot_bbox_buffered_all$ymax
cat("DEBUG: Plot Limits - x:", xmin_plot, xmax_plot, "y:", ymin_plot, ymax_plot, "\n")

# --- 9c. Define Consistent Aesthetics ---
occ_color <- "red"; occ_shape <- 17; occ_size <- 1.8; occ_alpha <- 0.8
bg_color <- "black"; bg_shape <- 1; bg_size <- 0.6; bg_alpha <- 0.5

# --- 9d. Base Plot Function ---
create_base_map <- function(world_data, raster_data, xlims, ylims, plot_crs) {
  p <- ggplot()
  if(!is.null(world_data)) { p <- p + geom_sf(data = world_data, fill = "grey80", color = "white", linewidth = 0.1) } # Use linewidth
  p <- p +
    geom_spatraster(data = raster_data) +
    scale_fill_viridis_c(option = "plasma", na.value = NA, name = names(raster_data)[1]) +
    coord_sf(crs = plot_crs, xlim = xlims, ylim = ylims, expand = FALSE) +
    theme_minimal(base_size = 10) +
    theme( panel.background = element_rect(fill = "aliceblue", color = NA), panel.grid.major = element_line(color = "grey90", linewidth = 0.2), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 9), legend.position = "right" )
  return(p)
}

# --- 9e. Plot for Method 1: Global Background ---
cat("  Generating Plot 1: Global Background\n")
bg_sf_global_plot <- if(!is.null(bg_global_df)) sf::st_as_sf(bg_global_df, coords = c("longitude", "latitude"), crs = plot_crs_sf) else NULL

plot_global <- create_base_map(world_sf, plot_raster_global, c(xmin_plot, xmax_plot), c(ymin_plot, ymax_plot), plot_crs_sf) +
  { if (!is.null(bg_sf_global_plot)) geom_sf(data = bg_sf_global_plot, color = bg_color, shape = bg_shape, size = bg_size, alpha = bg_alpha) } + # Fixed aesthetics
  geom_sf(data = occ_sf_plot, color = occ_color, shape = occ_shape, size = occ_size, alpha = occ_alpha) + # Fixed aesthetics
  labs(title = paste("Background Comparison:", species_name), subtitle = "Method 1: Global Background Points", x = "Longitude", y = "Latitude") +
  theme(legend.position = "right") # Only show raster legend

print(plot_global)
plot_save_dir <- "plots"; dir.create(plot_save_dir, showWarnings = FALSE)
plot_save_filename_global <- file.path(plot_save_dir, paste0("bg_plot_global_", species_name_sanitized, "_sdmtune.png"))
ggsave(plot_save_filename_global, plot = plot_global, width = 11, height = 7, dpi = 300, bg = "white")
cat(paste("  Plot 1 saved to:", plot_save_filename_global, "\n"))

# --- 9f. Plot for Method 2: Alpha Hull Background ---
if (!is.null(predictor_stack_species_hull)) {
  cat("  Generating Plot 2: Alpha Hull Background\n")
  bg_sf_species_hull_plot <- if(!is.null(bg_species_df_hull)) sf::st_as_sf(bg_species_df_hull, coords = c("longitude", "latitude"), crs = plot_crs_sf) else NULL
  alpha_hull_sf_plot <- if (!is.null(species_calibration_vect_hull)) sf::st_as_sf(species_calibration_vect_hull) else NULL
  if(!is.null(alpha_hull_sf_plot) && sf::st_crs(alpha_hull_sf_plot) != plot_crs_sf) alpha_hull_sf_plot <- sf::st_transform(alpha_hull_sf_plot, crs=plot_crs_sf)
  
  plot_hull <- create_base_map(world_sf, predictor_stack_species_hull[[1]], c(xmin_plot, xmax_plot), c(ymin_plot, ymax_plot), plot_crs_sf) +
    { if (!is.null(alpha_hull_sf_plot)) geom_sf(data = alpha_hull_sf_plot, fill = NA, color = "blue", linewidth = 0.8, aes(linetype = "Calibration Extent")) } +
    { if (!is.null(bg_sf_species_hull_plot)) geom_sf(data = bg_sf_species_hull_plot, color = bg_color, shape = bg_shape, size = bg_size, alpha = bg_alpha) } + # Fixed aesthetics
    geom_sf(data = occ_sf_plot, color = occ_color, shape = occ_shape, size = occ_size, alpha = occ_alpha) + # Fixed aesthetics
    scale_linetype_manual(name = "Extent", values = c("Calibration Extent" = "solid")) +
    guides(linetype = guide_legend(title = "Extent")) +
    labs(title = paste("Background Comparison:", species_name), subtitle = "Method 2: Alpha Hull Background Points & Extent", x = "Longitude", y = "Latitude") +
    theme(legend.position = "right")
  
  print(plot_hull)
  plot_save_filename_hull <- file.path(plot_save_dir, paste0("bg_plot_hull_", species_name_sanitized, "_sdmtune.png"))
  ggsave(plot_save_filename_hull, plot = plot_hull, width = 11, height = 7, dpi = 300, bg = "white")
  cat(paste("  Plot 2 saved to:", plot_save_filename_hull, "\n"))
} else { cat("WARN: Skipping Alpha Hull plot due to missing data.\n") }


# --- 9g. Plot for Method 4: OBIS Exact Background ---
if (!is.null(predictor_stack_obis_exact)) {
  cat("  Generating Plot 4: OBIS Exact Background\n")
  bg_sf_obis_exact_plot <- if(!is.null(bg_obis_exact_df)) sf::st_as_sf(bg_obis_exact_df, coords = c("longitude", "latitude"), crs = plot_crs_sf) else NULL
  obis_calibration_sf_plot <- if (!is.null(obis_calibration_vect_exact)) sf::st_as_sf(obis_calibration_vect_exact) else NULL
  if(!is.null(obis_calibration_sf_plot) && sf::st_crs(obis_calibration_sf_plot) != plot_crs_sf) obis_calibration_sf_plot <- sf::st_transform(obis_calibration_sf_plot, crs=plot_crs_sf)
  
  plot_obis_exact <- create_base_map(world_sf, predictor_stack_obis_exact[[1]], c(xmin_plot, xmax_plot), c(ymin_plot, ymax_plot), plot_crs_sf) +
    { if (!is.null(obis_calibration_sf_plot)) geom_sf(data = obis_calibration_sf_plot, fill = NA, color = "green", linewidth = 0.8, aes(linetype = "Calibration Extent")) } +
    { if (!is.null(bg_sf_obis_exact_plot)) geom_sf(data = bg_sf_obis_exact_plot, color = bg_color, shape = bg_shape, size = bg_size, alpha = bg_alpha) } + # Fixed aesthetics
    geom_sf(data = occ_sf_plot, color = occ_color, shape = occ_shape, size = occ_size, alpha = occ_alpha) + # Fixed aesthetics
    scale_linetype_manual(name = "Extent", values = c("Calibration Extent" = "dashed")) +
    guides(linetype = guide_legend(title = "Extent")) +
    labs(title = paste("Background Comparison:", species_name), subtitle = "Method 4: OBIS Ecoregion/Depth Background Points & Extent", x = "Longitude", y = "Latitude") +
    theme(legend.position = "right")
  
  print(plot_obis_exact)
  plot_save_filename_obis_exact <- file.path(plot_save_dir, paste0("bg_plot_obis_exact_", species_name_sanitized, "_sdmtune.png"))
  ggsave(plot_save_filename_obis_exact, plot = plot_obis_exact, width = 11, height = 7, dpi = 300, bg = "white")
  cat(paste("  Plot 4 saved to:", plot_save_filename_obis_exact, "\n"))
} else { cat("WARN: Skipping OBIS Exact plot due to missing data.\n") }


cat("\n--- compare_background_methods_sdmtune.R finished. ---\n")
#-------------------------------------------------------------------------------