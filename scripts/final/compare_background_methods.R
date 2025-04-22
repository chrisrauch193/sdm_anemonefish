# scripts/final/compare_background_methods_sdmtune.R
#-------------------------------------------------------------------------------
# Bare-bones script to compare background sampling strategies for ONE species
# using SDMtune workflow.
# 1. Global Background Sampling (using existing generate_sdm_background helper)
# 2. Species-Specific Background Sampling (Alpha Hull + new helper)
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

pacman::p_load(terra, sf, dplyr, readr, tools, stringr, SDMtune, ggplot2, tidyterra, rnaturalearth, viridis, concaveman) # Added concaveman

cat("--- Sourcing Helper Functions ---\n")
helper_paths <- c( file.path(config$helpers_dir, "env_processing_helpers.R"), file.path(config$helpers_dir, "sdm_modeling_helpers.R") )
missing_helpers <- helper_paths[!file.exists(helper_paths)]; if(length(missing_helpers) > 0) stop("Missing helper(s): ", paste(missing_helpers, collapse=", "))
invisible(sapply(helper_paths, source))

# --- 2. Define Target Species & Scenario ---
group_name <- "anemone"
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors
tuning_scenario <- "current"
cat("--- Target Group:", group_name, "---\n")
cat("--- Using Predictors:", ifelse(use_pca, "PCA", "VIF"), "---\n")

species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
species_row <- species_df[1, ] # <<< --- Test species --- >>>
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
} else { stop("VIF logic not implemented in this script.") }
if (is.null(global_predictor_stack)) stop("Failed to load GLOBAL predictor stack.")
if (terra::crs(global_predictor_stack) == "") stop("GLOBAL predictor stack missing CRS.")
cat("  GLOBAL stack loaded.\n")

# --- 4. Load & Prep Occurrences ---
cat("--- Loading and Cleaning Occurrences ---\n")
config_for_occ_load <- config; config_for_occ_load$predictor_stack_for_thinning <- global_predictor_stack
occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL)
if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) stop("Insufficient occurrences.")
occs_coords_df <- as.data.frame(occ_data_result$coords); colnames(occs_coords_df) <- c("longitude", "latitude")
occs_sf_clean <- sf::st_as_sf(occs_coords_df, coords=c("longitude","latitude"), crs=config$occurrence_crs) # Needed for new helper
cat("  Occurrence count:", nrow(occs_coords_df), "\n")

# --- 5. Method 1: Global Background & Tuning ---
cat("\n--- Running Method 1: GLOBAL Background Sampling & SDMtune ---\n")
bg_global_df <- generate_sdm_background(global_predictor_stack, config$background_points_n, config, logger = NULL, seed = 123)
if (is.null(bg_global_df)) stop("Failed global background generation.")
colnames(bg_global_df) <- c("longitude", "latitude")
cat("  Generated", nrow(bg_global_df), "global background points.\n")

cat("  Running SDMtune (Global Background)...\n")
tuning_output_global <- run_sdm_tuning_scv(occs_coords_df, global_predictor_stack, bg_global_df, config, logger=NULL, species_name=species_name)
if (is.null(tuning_output_global)) warning("SDMtune failed for Global Background.")

# --- 6. Method 2: Species-Specific Background & Tuning ---
cat("\n--- Running Method 2: SPECIES-SPECIFIC Background Sampling & SDMtune ---\n")
bg_species_result <- generate_sdm_background_species_extent(
  occs_sf = occs_sf_clean, # Pass the sf object
  global_predictor_stack = global_predictor_stack,
  config = config,
  logger = NULL,
  seed_offset = 1 # Use different offset for seed
)
if(is.null(bg_species_result)) stop("Failed species-specific background generation helper.")
bg_species_df <- bg_species_result$background_points
predictor_stack_species <- bg_species_result$species_specific_stack
species_calibration_vect <- bg_species_result$calibration_polygon
cat("  Generated", nrow(bg_species_df), "species-specific background points.\n")
cat("  Species-specific stack created and potentially depth-filtered.\n")

cat("  Running SDMtune (Species Background)...\n")
tuning_output_species <- run_sdm_tuning_scv(occs_coords_df, predictor_stack_species, bg_species_df, config, logger=NULL, species_name=species_name)
if (is.null(tuning_output_species)) warning("SDMtune failed for Species Background.")


# --- 7. Compare Results ---
cat("\n--- Comparing SDMtune Results ---\n")
comparison_list <- list()
if (!is.null(tuning_output_global) && inherits(tuning_output_global, "SDMtune")) {
  res_g <- tuning_output_global@results
  best_g_idx <- if(config$sdm_evaluation_metric == "AICc") which.min(res_g$AICc) else which.max(res_g[[paste0("test_", toupper(config$sdm_evaluation_metric))]])
  if(length(best_g_idx) == 0) best_g_idx <- 1 # Fallback
  comparison_list$GlobalBG <- res_g[best_g_idx, ]
} else { cat("WARN: Global tuning results missing or invalid.\n") }

if (!is.null(tuning_output_species) && inherits(tuning_output_species, "SDMtune")) {
  res_s <- tuning_output_species@results
  best_s_idx <- if(config$sdm_evaluation_metric == "AICc") which.min(res_s$AICc) else which.max(res_s[[paste0("test_", toupper(config$sdm_evaluation_metric))]])
  if(length(best_s_idx) == 0) best_s_idx <- 1 # Fallback
  comparison_list$SpeciesBG <- res_s[best_s_idx, ]
} else { cat("WARN: Species-specific tuning results missing or invalid.\n") }

if (length(comparison_list) > 0) {
  comparison_df <- dplyr::bind_rows(comparison_list, .id = "Method")
  metrics_to_show <- intersect(c("Method", "reg", "fc", "test_AUC", "test_TSS", "AICc"), names(comparison_df)) # Adjust metrics if needed
  cat("Comparison of Best Models (selected by", config$sdm_evaluation_metric, "):\n")
  print(comparison_df[, metrics_to_show], row.names = FALSE, digits = 4)
} else { cat("Could not generate comparison table.\n") }




# --- 8. Visualization ---
cat("\n--- Generating Comparison Plot ---\n")
cat("Generating Comparison Plot")

# --- 8a. Load Libraries (if not already loaded) ---
pacman::p_load(ggplot2, sf, terra, tidyterra, rnaturalearth, viridis)

# --- 8b. Prepare Data ---
# Base World Map
world_sf <- tryCatch(ne_countries(scale = "medium", returnclass = "sf"), error = function(e) { cat("Failed load rnaturalearth map."); NULL })

# Define the CRS from a valid object (e.g., the global stack)
plot_crs <- "4326" # terra::crs(global_predictor_stack)

# Convert points to sf objects, specifying correct column names
occ_sf_plot <- sf::st_as_sf(occs_coords_df, coords = c("longitude", "latitude"), crs = terra::crs(global_predictor_stack))

# **** CORRECTED BACKGROUND SF CONVERSION ****
# Use coords = c("x", "y") for background points generated by generate_sdm_background
bg_sf_global_plot <- sf::st_as_sf(bg_global_df, coords = c("longitude", "latitude"), crs = terra::crs(global_predictor_stack))
bg_sf_species_plot <- sf::st_as_sf(bg_species_df, coords = c("x", "y"), crs = terra::crs(global_predictor_stack))
# **** END CORRECTION ****

calibration_sf_plot <- if (!is.null(species_calibration_vect)) sf::st_as_sf(species_calibration_vect) else NULL

# --- 8c. Define Plot Extent ---
if (!is.null(calibration_sf_plot)) {
  # Ensure polygon CRS matches plot CRS before buffering
  # if(sf::st_crs(calibration_sf_plot) != sf::st_crs(plot_crs)) {
  #   calibration_sf_plot <- sf::st_transform(calibration_sf_plot, crs = plot_crs)
  # }
  # Project to estimate buffer distance accurately (using Mollweide centered on data)
  mean_lon_plot <- mean(st_coordinates(occ_sf_plot)[,1])
  proj_crs_plot <- paste0("+proj=moll +lon_0=", round(mean_lon_plot), " +datum=WGS84")
  tryCatch({
    calibration_sf_proj <- sf::st_transform(calibration_sf_plot, crs=proj_crs_plot)
    buffered_proj <- sf::st_buffer(calibration_sf_proj, dist = 200000) # 200km buffer example in meters
    plot_bbox_sf <- sf::st_transform(buffered_proj, crs = 4326)
    plot_bbox_coords <- sf::st_bbox(plot_bbox_sf)
  }, error = function(e){
    cat("Buffering polygon for plot extent failed, using occurrence bbox + 5 deg.");
    plot_bbox_coords <- sf::st_bbox(sf::st_buffer(occ_sf_plot, dist=5)) # Fallback: 5 degrees
  })
  
} else { # Fallback if polygon failed
  plot_bbox_coords <- sf::st_bbox(sf::st_buffer(occ_sf_plot, dist=5)) # 5-degree buffer around occurrences
}
xmin_plot <- plot_bbox_coords$xmin; xmax_plot <- plot_bbox_coords$xmax
ymin_plot <- plot_bbox_coords$ymin; ymax_plot <- plot_bbox_coords$ymax

# --- 8d. Create Plot ---
comparison_plot <- ggplot()
if(!is.null(world_sf)) { comparison_plot <- comparison_plot + geom_sf(data = world_sf, fill = "grey80", color = "white", size = 0.1) }
comparison_plot <- comparison_plot +
  # Use the SPECIES-SPECIFIC masked predictor stack for background visualization
  geom_spatraster(data = predictor_stack_species[[1]]) +
  scale_fill_viridis_c(option = "plasma", na.value = NA, name = names(predictor_stack_species)[1]) +
  { if (!is.null(calibration_sf_plot)) geom_sf(data = calibration_sf_plot, fill = NA, color = "blue", linewidth = 0.8, aes(linetype = "Calibration Extent")) } +
  geom_sf(data = bg_sf_global_plot, aes(shape = "Global BG", size = "Global BG", alpha = "Global BG"), color = "grey50") +
  geom_sf(data = bg_sf_species_plot, aes(shape = "Species BG", size = "Species BG", alpha = "Species BG"), color = "black") +
  geom_sf(data = occ_sf_plot, aes(shape = "Occurrence", size = "Occurrence", alpha = "Occurrence"), color = "red") +
  # Set Coordinate System and Limits using calculated bbox
  coord_sf(crs = 4326, # Use the defined plot CRS
           xlim = c(xmin_plot, xmax_plot),
           ylim = c(ymin_plot, ymax_plot),
           expand = FALSE) +
  # Manual scales for legend
  scale_shape_manual(name = "Data Type", values = c("Global BG" = 1, "Species BG" = 16, "Occurrence" = 16)) +
  scale_size_manual(name = "Data Type", values = c("Global BG" = 0.8, "Species BG" = 0.5, "Occurrence" = 2.0)) +
  scale_alpha_manual(name = "Data Type", values = c("Global BG" = 0.4, "Species BG" = 0.6, "Occurrence" = 0.8)) +
  scale_linetype_manual(name = "Extent", values = c("Calibration Extent" = "solid")) +
  guides( shape = guide_legend(title = "Data Type", order = 1, override.aes = list(alpha = 0.8)), size = "none", alpha = "none", linetype = guide_legend(title = "Extent", order = 2) ) +
  labs( title = paste("Background Sampling Comparison:", species_name), subtitle = "Species BG sampled within blue polygon & env data extent (+ depth filter if applied)", x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 12) +
  theme( panel.background = element_rect(fill = "aliceblue", color = NA), panel.grid.major = element_line(color = "grey90", linewidth = 0.2), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 9), legend.position = "right", legend.box = "vertical" )


# comparison_plot
 
# --- 8e. Print and Save Plot ---
print(comparison_plot)
log_dir_comp <- "plots"
plot_save_filename <- file.path(log_dir_comp, paste0("bg_comparison_sdmtune_", species_name_sanitized, ".png"))
ggsave(plot_save_filename, plot = comparison_plot, width = 10, height = 7, dpi = 300, bg = "white")
cat(paste("Comparison plot saved to:", plot_save_filename))

cat("\n--- compare_background_methods.R finished. ---\n")
#-------------------------------------------------------------------------------