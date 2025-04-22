# scripts/final/basic_run2.R

# --- Preamble: Cleanup ---
# source("scripts/final/cleanup.R", echo = TRUE) # Use relative path if needed

cat("--- Running basic_run2.R (with Alpha Hull Calibration Extent) ---\n")

# --- 1. Setup: Load Config FIRST ---
if (file.exists("../config.R")) {
  source("../config.R")
} else if (file.exists("scripts/config.R")) {
  source("scripts/config.R")
} else {
  stop("FATAL: Configuration file 'scripts/config.R' not found.")
}
if (!exists("config") || !is.list(config)) {
  stop("FATAL: 'config' list object not found or invalid after sourcing config.R")
}


# --- 2. Load Required Packages ---
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
# Added sf and concaveman
pacman::p_load(terra, sf, dplyr, readr, tools, stringr, log4r, future, furrr, progressr, ENMeval, predicts, concaveman)


# --- 3. Source Helper Functions ---
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers_enmeval.R"))


# --- 4. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path, log_level = config$log_level, append = config$log_append, log_to_console = config$log_to_console, console_level = config$log_console_level)
log4r::info(logger, "--- Starting basic_run2.R (with Alpha Hull Calibration Extent) ---")


# --- 5. Define Group Specifics & Predictor Type ---
group_name <- "anemone" # Example: testing with anemones first
config$group_name <- group_name # Ensure config has this for helpers
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors # Assuming PCA is used based on 05b
predictor_type_suffix <- ifelse(use_pca, "_enmeval_pca", "_enmeval_vif")
log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---"))
# ... (log other ENMeval settings as before) ...


# --- 6. Load Predictor Information (Global PCA Paths) ---
predictor_paths_or_list <- NULL
if (use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (is.null(pca_paths_rds) || !file.exists(pca_paths_rds)) { log4r::fatal(logger, "PCA paths RDS missing/invalid."); stop("PCA paths RDS invalid.") }
  predictor_paths_or_list <- tryCatch(readRDS(pca_paths_rds), error = function(e) { log4r::fatal(logger, paste("Failed load PCA paths RDS:", e$message)); NULL })
  if (!is.list(predictor_paths_or_list) || length(predictor_paths_or_list) == 0) { log4r::fatal(logger, "PCA paths list empty/invalid."); stop("PCA paths list invalid.") }
  log4r::info(logger, paste("Loaded GLOBAL PCA raster paths for scenarios:", paste(names(predictor_paths_or_list), collapse = ", ")))
} else {
  # Add VIF logic here if needed for the test, otherwise assume PCA
  log4r::fatal(logger, "This test script currently assumes use_pca_predictors = TRUE.")
  stop("VIF logic not implemented in this basic test script.")
}

# --- 7. Create Intermediate Output Dirs ---
# ... (keep as is) ...
base_intermediate_model_path <- config$models_dir_intermediate
base_intermediate_results_path <- config$results_dir_intermediate
base_species_log_path <- config$species_log_dir
intermediate_models_dir <- file.path(base_intermediate_model_path, paste0(group_name, predictor_type_suffix))
intermediate_results_dir <- file.path(base_intermediate_results_path, paste0(group_name, predictor_type_suffix))
tryCatch({ dir.create(intermediate_models_dir, recursive = TRUE, showWarnings = FALSE); dir.create(intermediate_results_dir, recursive = TRUE, showWarnings = FALSE); dir.create(base_species_log_path, recursive = TRUE, showWarnings = FALSE); log4r::debug(logger, paste("Intermediate dirs created/checked:", group_name, predictor_type_suffix))
}, error = function(e) { log4r::fatal(logger, paste("Failed create intermediate dirs:", e$message)); stop("Directory creation failed.") })


# --- 8. Load Species List & Select Test Species ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE) }, error = function(e) { log4r::fatal(logger, paste("Failed load species list:", e$message)); stop("Species list failed.") })
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))

# --- Select a species to test (e.g., the first one) ---
species_row <- species_df[1, ]
species_name <- species_row$scientificName
species_name_sanitized <- gsub(" ", "_", species_name)
species_aphia_id <- species_row$AphiaID

# --- Setup Species Log ---
species_log_file <- file.path(config$species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_basic_run2_detail.log"))
slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[", level, "]"), paste0("[", species_name, "]"), paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
slog("INFO", paste0("--- Starting basic ENMeval processing (", predictor_type_suffix, ") with Alpha Hull Extent ---"))

# --- Load Global Predictors for Tuning ---
tuning_scenario <- "current"
slog("INFO", "Loading GLOBAL predictor stack for tuning scenario:", tuning_scenario)
tuning_predictor_stack_global <- NULL
if (use_pca) {
  tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
  if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) { msg <- paste0("FATAL: GLOBAL PCA stack path for tuning scenario '", tuning_scenario, "' not found."); slog("ERROR", msg); stop(msg) }
  tuning_predictor_stack_global <- tryCatch(terra::rast(tuning_predictor_path), error = function(e) { slog("ERROR", "Failed to load GLOBAL tuning PCA stack:", e$message); NULL })
} else {
  # VIF logic would go here if needed
  msg <- "VIF predictor loading not implemented in this basic script."; slog("ERROR", msg); stop(msg)
}
if (is.null(tuning_predictor_stack_global)) { msg <- "FATAL: Failed to load GLOBAL predictor stack for tuning."; slog("ERROR", msg); stop(msg) }
if (terra::crs(tuning_predictor_stack_global) == "") { msg <- "FATAL: GLOBAL Tuning predictor stack has no CRS."; slog("ERROR", msg); stop(msg) }
slog("INFO", "GLOBAL Tuning predictor stack loaded. Names:", paste(names(tuning_predictor_stack_global), collapse = ", "))


# --- Load Occurrences ---
config_for_occ_load <- config
config_for_occ_load$predictor_stack_for_thinning <- tuning_predictor_stack_global # Use global stack for initial thinning extent
slog("INFO", "Loading/cleaning/thinning occurrences.")
occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger = NULL, species_log_file = species_log_file)
if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) {
  msg <- paste0("FATAL: Insufficient valid/thinned occurrences (found ", occ_data_result$count %||% 0, "). Required: ", config$min_occurrences_sdm); slog("ERROR", msg); stop(msg)
}
occs_coords_df <- as.data.frame(occ_data_result$coords); colnames(occs_coords_df) <- c("longitude", "latitude")
occurrence_count_after_thinning <- occ_data_result$count
slog("INFO", "Occurrence count after clean/thin:", occurrence_count_after_thinning)

# --- Create Species-Specific Calibration Extent (Buffered Alpha Hull) ---
slog("INFO", "Creating species-specific calibration extent (Buffered Alpha Hull)...")
species_calibration_vect <- NULL
tryCatch({
  occ_sf_species <- sf::st_as_sf(occs_coords_df, coords = c("longitude", "latitude"), crs = config$occurrence_crs)
  
  # 1. Project for accurate buffering (Choose a suitable projection!)
  # Example: Mollweide centered near Indo-Pacific. Replace if needed.
  target_crs_proj <- "+proj=moll +lon_0=125 +datum=WGS84"
  slog("DEBUG", paste("Projecting occurrences to:", target_crs_proj))
  occ_sf_proj <- sf::st_transform(occ_sf_species, crs = target_crs_proj)
  
  # 2. Calculate Alpha Hull (or fallback)
  concavity_value <- 2.5 # Example value
  buffer_km <- 150     # Example buffer distance
  slog("DEBUG", paste("Calculating hull/polygon with concavity =", concavity_value))
  
  if (nrow(occ_sf_proj) >= 5) {
    ahull_poly_proj <- concaveman::concaveman(occ_sf_proj, concavity = concavity_value)
  } else if (nrow(occ_sf_proj) >= 3) {
    slog("WARN", "Fewer than 5 points, using Convex Hull instead of Alpha Hull.")
    ahull_poly_proj <- sf::st_convex_hull(sf::st_union(occ_sf_proj))
  } else {
    slog("WARN", "Fewer than 3 points, using simple circular buffer.")
    ahull_poly_proj <- sf::st_union(sf::st_buffer(occ_sf_proj, dist = buffer_km * 1000))
  }
  slog("DEBUG", "Raw polygon calculated in projected CRS.")
  
  # 3. Buffer the Polygon
  slog("DEBUG", paste("Buffering polygon by", buffer_km, "km."))
  ahull_buffered_proj <- sf::st_buffer(ahull_poly_proj, dist = buffer_km * 1000) # Buffer takes meters
  
  # 4. Project back to Raster CRS
  raster_crs <- sf::st_crs(terra::crs(tuning_predictor_stack_global))
  slog("DEBUG", "Projecting buffered polygon back to raster CRS:", raster_crs$input)
  species_calibration_poly <- sf::st_transform(ahull_buffered_proj, crs = raster_crs)
  
  # 5. Convert to SpatVector
  species_calibration_vect <- terra::vect(species_calibration_poly)
  slog("INFO", "Species-specific calibration vector created.")
  
}, error = function(e) {
  msg <- paste("FATAL: Failed to create calibration polygon:", e$message); slog("ERROR", msg); stop(msg)
})
if (is.null(species_calibration_vect)) stop("Polygon creation failed.")

# --- Mask Global Tuning Stack with Species Polygon ---
slog("INFO", "Masking GLOBAL tuning predictors with species polygon...")
tuning_predictor_stack_species <- tryCatch({
  cropped_stack <- terra::crop(tuning_predictor_stack_global, species_calibration_vect, snap = "near")
  terra::mask(cropped_stack, species_calibration_vect)
}, error = function(e) {
  slog("ERROR", paste("Failed to mask tuning stack:", e$message)); NULL
})
if (is.null(tuning_predictor_stack_species)) stop("Masking tuning stack failed.")
min_max_check <- terra::global(tuning_predictor_stack_species[[1]], c("min", "max"), na.rm = TRUE)
if (is.na(min_max_check$min) && is.na(min_max_check$max)) {
  msg <- "FATAL: No valid cells remain in tuning stack after masking."; slog("ERROR", msg); stop(msg)
}
slog("INFO", "GLOBAL tuning stack successfully masked to species extent.")

# --- Apply Depth Filter (Optional) ---
if (config$apply_depth_filter) {
  slog("INFO", "Applying depth filter...")
  depth_layer_name <- config$terrain_variables_final[grepl("bathymetry", config$terrain_variables_final, ignore.case = TRUE)][1]
  
  bath_path <- file.path(config$terrain_folder, paste0(depth_layer_name, ".tif"))
  if (file.exists(bath_path)) {
    bath_layer_global <- tryCatch(terra::rast(bath_path), error = function(e) { slog("WARN", "Failed load bathy layer:", e$message); NULL })
    if (!is.null(bath_layer_global)) {
      
      # --- START CORRECTED DEPTH MASKING BLOCK ---
      if (terra::crs(bath_layer_global) != terra::crs(tuning_predictor_stack_species)) {
        slog("WARN", "Projecting GLOBAL bathymetry layer to match species stack CRS.")
        bath_layer_global <- tryCatch(terra::project(bath_layer_global, tuning_predictor_stack_species), error = function(e) { slog("ERROR", "Failed projecting global bathy:", e$message); NULL })
      }
      
      if (!is.null(bath_layer_global)) {
        bath_layer_species <- tryCatch(
          terra::mask(
            terra::crop(bath_layer_global, species_calibration_vect, snap = "near"),
            species_calibration_vect
          ), error = function(e) { slog("ERROR", paste("Failed crop/mask bathymetry:", e$message)); NULL }
        )
        
        if (!is.null(bath_layer_species)) {
          depth_mask_logical <- bath_layer_species >= config$depth_min & bath_layer_species <= config$depth_max
          
          slog("DEBUG", "Aligning depth mask raster with species predictor stack...")
          # Use resample to force alignment
          depth_mask_aligned <- tryCatch(
            terra::resample(depth_mask_logical, tuning_predictor_stack_species, method = "near"),
            error = function(e) { slog("ERROR", paste("Failed aligning depth mask:", e$message)); NULL }
          )
          
          if (!is.null(depth_mask_aligned)) {
            tuning_predictor_stack_species_before_depth <- tuning_predictor_stack_species # Keep a copy for check
            # Mask where condition is FALSE (i.e., outside depth range)
            tuning_predictor_stack_species <- terra::mask(tuning_predictor_stack_species, depth_mask_aligned, maskvalues = FALSE, updatevalue = NA)
            
            slog("INFO", "Depth filter applied using aligned mask.")
            
            min_max_check_depth <- terra::global(tuning_predictor_stack_species[[1]], c("min", "max"), na.rm = TRUE) # Check the correct object
            if (is.na(min_max_check_depth$min) && is.na(min_max_check_depth$max)) {
              msg <- "FATAL: No valid cells remain in tuning stack after depth filter."; slog("ERROR", msg); stop(msg)
            } else {
              cells_before <- terra::global(!is.na(tuning_predictor_stack_species_before_depth[[1]]), "sum", na.rm = TRUE)$sum
              cells_after <- terra::global(!is.na(tuning_predictor_stack_species[[1]]), "sum", na.rm = TRUE)$sum
              slog("DEBUG", paste("Cells before depth mask:", cells_before, " | Cells after depth mask:", cells_after))
            }
            rm(tuning_predictor_stack_species_before_depth) # Clean up
          } else { slog("WARN", "Skipping depth filter application due to alignment error.") }
          rm(depth_mask_logical, depth_mask_aligned) # Clean up
        } else { slog("WARN", "Skipping depth filter application because species-specific bathymetry failed.") }
        rm(bath_layer_species) # Clean up
      } else { slog("WARN", "Skipping depth filter application because global bathymetry projection failed.") }
      # --- END CORRECTED DEPTH MASKING BLOCK ---
      
    } else { slog("WARN", "Failed loading bathymetry, skipping depth filter.") }
    if (exists("bath_layer_global")) rm(bath_layer_global) # Clean up
  } else { slog("WARN", "Bathymetry layer not found at expected path, skipping depth filter:", bath_path) }
} else {
  slog("INFO", "Depth filtering disabled in config.")
}


# --- Generate Background Points from Species-Specific Extent ---
slog("INFO", "Generating background points from SPECIES-SPECIFIC masked extent...")
background_points_df_species <- generate_sdm_background(
  predictor_stack = tuning_predictor_stack_species, # USE THE MASKED & potentially depth-filtered STACK
  n_background = config$background_points_n,
  config = config,
  logger = NULL,
  species_log_file = species_log_file,
  seed = species_aphia_id
)
if (is.null(background_points_df_species)) { msg <- paste0("FATAL: Failed background point generation within species extent."); slog("ERROR", msg); stop(msg) }
colnames(background_points_df_species) <- c("longitude", "latitude")
slog("INFO", "Background points generated from species extent:", nrow(background_points_df_species))


# --- Run ENMevaluate with Species-Specific Data ---
# # Prepare arguments list, using the species-specific stacks and background points
# enmeval_args <- list(
#   occs = occs_coords_df,
#   envs = tuning_predictor_stack_species, # Use species-specific stack
#   bg = background_points_df_species,    # Use species-specific background
#   tune.args = config$enmeval_tuning_settings,
#   algorithm = config$enmeval_algorithms,
#   partitions = config$enmeval_partitions,
#   other.settings = list(pred.type = config$enmeval_pred_type, abs.auc.diff = FALSE, validation.bg = "partition"),
#   clamp = config$enmeval_clamp,
#   parallel = config$use_parallel && config$enmeval_num_cores > 1,
#   numCores = if (config$use_parallel && config$enmeval_num_cores > 1) config$enmeval_num_cores else 1,
#   quiet = FALSE # Turn quiet off for this test run to see messages
# )



# prep data for custom block sampling
pts <- rbind(occs_coords_df, background_points_df_species)
pts$occ <- c(rep(1, nrow(occs_coords_df)), rep(0, nrow(background_points_df_species)))

# investigate spatial autocorrelation in the landscape to choose a suitable size for spatial blocks
# folds.cor <- cv_spatial_autocor(r = tuning_predictor_stack, x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326), column = 'occ', num_sample = 10000, plot = T, progress = T)

# generate folds
# scv <- cv_spatial(x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326), column = 'occ', r = tuning_predictor_stack, k = 5, hexagon = T, flat_top = F, size = folds.cor$range,
#                   selection = 'random', iteration = 50, progress = T, report = T, plot = T, raster_colors = terrain.colors(10, rev = T))
scv <- cv_spatial(x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326), column = 'occ', r = tuning_predictor_stack_species, k = 5, hexagon = T, flat_top = F, size = 20000,
                  selection = 'random', iteration = 50, progress = T, report = T, plot = T, raster_colors = terrain.colors(10, rev = T))

# plot folds
cv_plot(cv = scv, x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326))

# separate occs and bg folds
fold_ids <- data.frame(fold_ids = scv$folds_ids)
folds <-  cbind(pts, fold_ids)
head(folds)

occs.folds <- folds %>% dplyr::filter(occ == 1) %>% dplyr::select('fold_ids')
bg.folds <- folds %>% dplyr::filter(occ == 0) %>% dplyr::select('fold_ids')

# export folds
saveRDS(as.vector(occs.folds)$fold_ids, 'outputs/folds/occs_fold.rds')
saveRDS(as.vector(bg.folds)$fold_ids, 'outputs/folds/bg_fold.rds')





enmeval_results <- ENMevaluate(taxon.name = 'Radianthus magnifica',
                               occs = occs_coords_df,
                               envs = tuning_predictor_stack_species,
                               bg = background_points_df_species,
                               tune.args = list(rm = seq(0.5, 5, by = 0.5),
                                                fc = c('L', 'LQ', 'H', 'LQH', 'LQHP', 'LQHPT')),
                               partitions = 'user',
                               user.grp = list(occs.grp = as.vector(occs.folds)$fold_ids,
                                               bg.grp = as.vector(bg.folds)$fold_ids),
                               algorithm = 'maxnet',
                               doClamp = T,
                               parallel = config$use_parallel && config$num_cores > 1,
                               numCores = if(config$use_parallel && config$num_cores > 1) config$num_cores else 1,
                               updateProgress = T)


i.mods <- enmeval_results

i.tune.res <- eval.results(i.mods)
i.find.opt <- i.tune.res %>%
  dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
  dplyr::filter(cbi.val.avg == max(cbi.val.avg)) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg)) %>%
  print()