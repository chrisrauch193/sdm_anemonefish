# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow v2)
# - Saving logic moved to dedicated helper functions.
# - Path construction uses new config structure for target output.
# - Added Variable Importance helper.
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr, blockCV, ggplot2) # Added readr

#' Load, Clean, and Get Coordinates for Individual Species Occurrences
#' (Returns list with data and count)
#'
#' @param species_aphia_id AphiaID of the target species.
#' @param occurrence_dir Directory containing individual species CSV files.
#' @param config Project configuration list (needs `occurrence_crs`, `thinning_method`, `min_occurrences_sdm`, and optionally `predictor_stack_for_thinning` if thinning is 'cell').
#' @param logger A log4r logger object (can be NULL for less verbose output).
#' @param species_log_file Optional path to a species-specific log file for detailed messages.
#' @return A list containing cleaned coordinates `coords` (matrix) and the final `count`, or NULL on error.
load_clean_individual_occ_coords <- function(species_aphia_id, occurrence_dir, config, logger, species_log_file = NULL) {
  
  # Local logging function for this helper
  hlog <- function(level, ...) {
    msg <- paste(Sys.time(), paste0("[",level,"]"), "[OccLoadHelper]", paste0(..., collapse = " "))
    if (!is.null(species_log_file)) { # Prefer species log if provided
      cat(msg, "\n", file = species_log_file, append = TRUE)
    } else if (!is.null(logger)) { # Fallback to main logger
      log_level_func <- switch(level, DEBUG=log4r::debug, INFO=log4r::info, WARN=log4r::warn, log4r::error)
      log_level_func(logger, msg) # Log to main logger
    } else {
      # cat(msg, "\n") # Fallback to console if no loggers provided
    }
  }
  
  hlog("DEBUG", paste("Loading occurrences for AphiaID:", species_aphia_id))
  occ_file <- file.path(occurrence_dir, paste0(species_aphia_id, ".csv"))
  if (!file.exists(occ_file)) {
    hlog("WARN", paste("Occurrence file not found:", basename(occ_file)))
    return(list(coords = NULL, count = 0))
  }
  
  tryCatch({
    occ_df <- readr::read_csv(occ_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
    
    if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occ_df))) {
      hlog("WARN", paste("Missing coordinate columns in file:", basename(occ_file)))
      return(list(coords = NULL, count = 0))
    }
    
    occ_clean <- occ_df %>%
      dplyr::mutate(
        decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
        decimalLatitude = suppressWarnings(as.numeric(decimalLatitude))
      ) %>%
      dplyr::filter(
        !is.na(decimalLongitude), !is.na(decimalLatitude),
        decimalLongitude >= -180, decimalLongitude <= 180,
        decimalLatitude >= -90, decimalLatitude <= 90
      )
    
    count_after_clean <- nrow(occ_clean)
    hlog("DEBUG", paste("  Retained", count_after_clean, "records after coordinate cleaning."))
    
    if (count_after_clean == 0) {
      hlog("WARN", "No valid coordinates after cleaning.")
      return(list(coords = NULL, count = 0))
    }
    
    # --- Spatial Thinning (Optional) ---
    if (!is.null(config$thinning_method) && config$thinning_method == "cell" && !is.null(config$predictor_stack_for_thinning)) {
      hlog("DEBUG", "  Applying cell-based thinning...")
      predictor_stack_thin <- config$predictor_stack_for_thinning
      if (is.null(predictor_stack_thin) || !inherits(predictor_stack_thin, "SpatRaster") || terra::nlyr(predictor_stack_thin) == 0) {
        hlog("WARN", " predictor_stack_for_thinning invalid/missing. Skipping thinning.")
        occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      } else {
        occs_sf <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs)
        target_crs <- terra::crs(predictor_stack_thin)
        if (sf::st_crs(occs_sf) != sf::st_crs(target_crs)) {
          occs_sf <- sf::st_transform(occs_sf, crs = target_crs)
        }
        occs_spatvector <- terra::vect(occs_sf)
        occs_cells <- tryCatch(terra::extract(predictor_stack_thin[[1]], occs_spatvector, cells = TRUE), error=function(e){hlog("WARN",paste("Error extracting cells:",e$message));NULL})
        
        if(is.null(occs_cells)){
          hlog("WARN", " Cell extraction failed. Skipping thinning.")
          occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
        } else {
          valid_cells_indices <- which(!is.na(occs_cells$cell))
          if(length(valid_cells_indices) > 0) {
            unique_cell_indices_in_valid <- which(!duplicated(occs_cells$cell[valid_cells_indices]))
            original_indices_to_keep <- valid_cells_indices[unique_cell_indices_in_valid]
            occs_thinned_sf <- occs_sf[original_indices_to_keep, ]
            occ_thinned_coords <- sf::st_coordinates(occs_thinned_sf)
            colnames(occ_thinned_coords) <- c("decimalLongitude", "decimalLatitude") # Ensure names
          } else {
            hlog("WARN", " No valid cells found for occurrences on thinning raster. Skipping thinning.")
            occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
          }
        }
      }
      count_after_thin <- nrow(occ_thinned_coords)
      hlog("DEBUG", paste("  Retained", count_after_thin, "records after thinning."))
      if (count_after_thin == 0) {hlog("WARN", "No records left after thinning."); return(list(coords = NULL, count = 0))}
    } else {
      hlog("DEBUG", "  Skipping spatial thinning or method not 'cell'.")
      occ_thinned_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
      count_after_thin <- nrow(occ_thinned_coords)
    }
    
    return(list(coords = occ_thinned_coords, count = count_after_thin))
    
  }, error = function(e) {
    hlog("ERROR", paste("Error processing occurrences for AphiaID", species_aphia_id, ":", e$message))
    return(list(coords = NULL, count = 0))
  })
}


#' Generate Background Points within the Raster Extent (v2 - with Coral Masking)
#' Can optionally mask sampling to coral reef areas defined in config.
#' @param predictor_stack SpatRaster stack (used for extent/masking).
#' @param n_background Number of points to attempt generating.
#' @param config Project configuration list (needs `mask_background_points_to_coral`, `apply_coral_mask`, `coral_shapefile`). #<< ADDED config
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @param seed Optional random seed for reproducibility.
#' @return A data frame of background point coordinates (x, y), or NULL on error.
generate_sdm_background <- function(predictor_stack, n_background, config, logger, species_log_file = NULL, seed = 123) { #<< ADDED config
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[BgGenHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  hlog("DEBUG", paste("Generating up to", n_background, "background points...")) # Changed log level
  if (is.null(predictor_stack) || !inherits(predictor_stack, "SpatRaster") || terra::nlyr(predictor_stack) == 0) { hlog("ERROR", "Predictor stack is required."); return(NULL) }
  if (!is.null(seed)) set.seed(seed)
  
  # --- Determine the layer to sample from ---
  sampling_layer <- predictor_stack[[1]] # Start with the first layer
  sampling_mask <- NULL # Initialize
  
  # Apply Coral Mask (if configured in config)
  if (config$mask_background_points_to_coral && config$apply_coral_mask) {
    hlog("DEBUG", "  Attempting coral reef mask for background point sampling...")
    if (!is.null(config$coral_shapefile) && file.exists(config$coral_shapefile)) {
      coral_areas_sf <- tryCatch({ sf::st_read(config$coral_shapefile, quiet = TRUE) }, error = function(e) {hlog("WARN",paste("   Failed load coral shapefile:", e$message)); NULL})
      if (!is.null(coral_areas_sf)) {
        coral_areas_vect <- tryCatch({ terra::vect(coral_areas_sf) }, error = function(e) {hlog("WARN",paste("   Failed convert coral sf to vect:", e$message)); NULL})
        if (!is.null(coral_areas_vect)) {
          if(terra::crs(coral_areas_vect) != terra::crs(sampling_layer)){
            hlog("DEBUG", "    Projecting coral shapefile CRS...")
            coral_areas_vect <- tryCatch(terra::project(coral_areas_vect, terra::crs(sampling_layer)), error = function(e){hlog("WARN",paste("     Failed project coral shapefile:", e$message)); NULL})
          }
          if(!is.null(coral_areas_vect)){
            sampling_mask <- tryCatch(terra::mask(sampling_layer, coral_areas_vect), error=function(e){hlog("WARN",paste("     Failed mask sampling layer:", e$message)); NULL})
            if(!is.null(sampling_mask)) { hlog("DEBUG", "  Coral reef mask applied for sampling.") }
            else { hlog("WARN", "  Failed to create mask. Sampling from unmasked layer.") }
          } else { hlog("WARN", "  CRS projection failed. Sampling from unmasked layer.") }
        } else { hlog("WARN", "  sf to vect conversion failed. Sampling from unmasked layer.") }
      } else { hlog("WARN", "  Coral shapefile loading failed. Sampling from unmasked layer.") }
    } else { hlog("WARN", "  Coral shapefile path missing/invalid. Sampling from unmasked layer.") }
  } else { hlog("DEBUG", "  Background point sampling not masked to coral reefs.") }
  
  # Use the masked layer if created, otherwise the original first layer
  layer_to_sample_from <- if(!is.null(sampling_mask)) sampling_mask else sampling_layer
  
  # --- Sample Points ---
  tryCatch({
    bg_points <- terra::spatSample(layer_to_sample_from, size = n_background, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
    if (nrow(bg_points) < n_background) { hlog("WARN", paste("Could only sample", nrow(bg_points), "background points (requested", n_background, ", likely due to mask/raster extent).")) }
    if (nrow(bg_points) == 0) { hlog("ERROR", "Failed to generate ANY background points from the sampling area."); return(NULL) }
    hlog("INFO", paste("Generated", nrow(bg_points), "background points.")) # Changed log level
    return(as.data.frame(bg_points[, c("x", "y")]))
  }, error = function(e) { hlog("ERROR", paste("Error generating background points:", e$message)); return(NULL) })
}


#' Generate Background Points using OBIS/MPAEU Ecoregion/Depth Method (Exact Replication v3)
#'
#' Replicates the OBIS/MPAEU Part 2 spatial extent definition AND background point sampling.
#' 1. Finds intersecting/adjacent ecoregions based on occurrences.
#' 2. Optionally buffers the ecoregion polygon.
#' 3. Optionally filters by depth range.
#' 4. Masks global predictors using the final spatial extent.
#' 5. Calculates background point number (quad_n) based on OBIS logic.
#' 6. Samples background points within the final masked extent.
#'
#' @param occs_sf sf object with cleaned, thinned occurrence points (needs CRS).
#' @param global_predictor_stack SpatRaster stack covering the global study area (needs CRS).
#' @param config Project configuration list. Requires relevant paths, flags, and
#'        `background_points_n` (used as the base target number, mimicking `quad_samp`).
#' @param seed_offset Integer offset for background sampling seed.
#' @return A list containing `$background_points` (data frame of coordinates),
#'         `$species_specific_stack` (SpatRaster, masked to final extent),
#'         `$calibration_polygon` (SpatVector, pre-depth mask extent),
#'         `$quad_n_calculated` (Numeric, the number of background points sampled),
#'         or NULL on error.
generate_sdm_background_obis_test <- function(occs_sf, global_predictor_stack, config, seed_offset = 4) {
  
  # Use cat for direct output in this comparison script
  cat("INFO: Starting OBIS/MPAEU (Part 2 Replica v3) background generation method...\n")
  
  # --- Input Checks ---
  if (!inherits(occs_sf, "sf") || nrow(occs_sf) < 3) { cat("ERROR: Valid sf occurrence object with >= 3 points required.\n"); return(NULL) }
  if (is.na(sf::st_crs(occs_sf))) { cat("ERROR: Occurrence sf object needs a CRS.\n"); return(NULL) }
  if (!inherits(global_predictor_stack, "SpatRaster") || terra::nlyr(global_predictor_stack) < 1) { cat("ERROR: Valid global SpatRaster required.\n"); return(NULL) }
  if (terra::crs(global_predictor_stack) == "") { cat("ERROR: Global predictor stack needs a CRS.\n"); return(NULL) }
  if (!file.exists(config$ecoregion_shapefile)) { cat("ERROR: Ecoregion shapefile not found:", config$ecoregion_shapefile, "\n"); return(NULL) }
  if (config$limit_by_depth_obis && !file.exists(config$bathymetry_file)) { cat("ERROR: Bathymetry file required for depth filtering not found:", config$bathymetry_file, "\n"); return(NULL) }
  
  species_aphia_id_internal <- tryCatch(occs_sf$AphiaID[1], error = function(e) { cat("WARN: AphiaID missing from occs_sf, using random seed.\n"); sample.int(1e6, 1) })
  
  # --- 1. Ecoregion Extent Definition ---
  # ... (Keep the ecoregion loading, intersection, adjacency, buffering logic from the previous version) ...
  # --- Load Ecoregions ---
  cat("DEBUG: Loading ecoregions...\n")
  ecoregions <- tryCatch(terra::vect(config$ecoregion_shapefile), error = function(e) { cat("ERROR: Failed load ecoregions:", e$message, "\n"); NULL })
  if (is.null(ecoregions)) return(NULL)
  occs_vect <- tryCatch(terra::vect(occs_sf), error = function(e) { cat("ERROR: Failed convert occs sf to vect:", e$message, "\n"); NULL })
  if(is.null(occs_vect)) return(NULL)
  if (terra::crs(occs_vect) != terra::crs(ecoregions)) {
    cat("DEBUG: Projecting occurrences to match ecoregions CRS...\n")
    occs_vect <- tryCatch(terra::project(occs_vect, terra::crs(ecoregions)), error = function(e) { cat("ERROR: Failed project occurrences:", e$message, "\n"); NULL })
    if(is.null(occs_vect)) return(NULL)
  }
  # --- Identify Intersecting & Adjacent Ecoregions (OBIS Step) ---
  cat("DEBUG: Identifying intersecting ecoregions...\n")
  intersect_idx <- tryCatch(terra::is.related(ecoregions, occs_vect, "intersects"), error=function(e){cat("WARN: is.related failed:",e$message,"\n");NULL})
  if (is.null(intersect_idx) || !any(intersect_idx)) { cat("WARN: No occurrences intersect with ecoregions. Cannot define extent this way.\n"); return(NULL) }
  ecoreg_occ <- ecoregions[intersect_idx, ]
  cat("DEBUG: Finding adjacent ecoregions using buffered points...\n")
  sf::sf_use_s2(FALSE)
  occs_buffered_sf <- tryCatch(sf::st_buffer(sf::st_as_sf(occs_vect), dist = config$poly_buffer_obis), error = function(e) { cat("ERROR: Failed buffer points:", e$message, "\n"); NULL })
  if(is.null(occs_buffered_sf)) { sf::sf_use_s2(TRUE); return(NULL) }
  adjacent_idx <- tryCatch(terra::is.related(ecoregions, terra::vect(occs_buffered_sf), "intersects"), error=function(e){cat("WARN: is.related failed for adjacent:",e$message,"\n");NULL})
  if(is.null(adjacent_idx)) adjacent_idx <- intersect_idx # Fallback
  sf::sf_use_s2(TRUE)
  final_indices <- unique(c(which(intersect_idx), which(adjacent_idx)))
  ecoreg_sel <- ecoregions[final_indices, ]
  cat("DEBUG: Buffering final selected ecoregion polygon...\n")
  sf::sf_use_s2(FALSE)
  ecoreg_sel_buffered_sf <- tryCatch(sf::st_buffer(sf::st_as_sf(ecoreg_sel), dist = config$poly_buffer_final), error = function(e) { cat("ERROR: Failed buffer final polygon:", e$message, "\n"); NULL })
  if(is.null(ecoreg_sel_buffered_sf)) { sf::sf_use_s2(TRUE); return(NULL) }
  sf::sf_use_s2(TRUE)
  raster_crs_terra <- terra::crs(global_predictor_stack)
  obis_calibration_poly_sf <- sf::st_transform(ecoreg_sel_buffered_sf, crs = sf::st_crs(raster_crs_terra))
  obis_calibration_vect <- terra::vect(obis_calibration_poly_sf)
  cat("INFO: Ecoregion-based calibration polygon created (before depth filter).\n")
  
  # --- 2. Depth Filtering (Optional) ---
  # ... (Keep the depth filtering logic exactly as corrected in the previous response) ...
  depth_mask <- NULL
  limit_by_depth_flag <- config$limit_by_depth_obis
  if (limit_by_depth_flag) {
    cat("INFO: Applying OBIS depth filter...\n")
    bath_global <- tryCatch(terra::rast(config$bathymetry_file), error = function(e) { cat("ERROR: Failed load bathymetry:", e$message, "\n"); NULL })
    if (is.null(bath_global)) { limit_by_depth_flag <- FALSE; cat("WARN: Disabling depth filter because bathymetry loading failed.\n") }
    else {
      if (config$apply_indo_pacific_crop) { cat("DEBUG: Applying Indo-Pacific crop to global bathymetry layer...\n"); ip_extent <- terra::ext(config$indo_pacific_bbox); bath_global <- tryCatch(terra::crop(bath_global, ip_extent), error=function(e){ cat("WARN: Failed cropping global bathy:", e$message, "\n"); NULL}); if(is.null(bath_global)){ cat("WARN: Proceeding without depth filter due to bathy crop error.\n"); limit_by_depth_flag <- FALSE } else { cat("DEBUG: Global bathymetry cropped to Indo-Pacific extent.\n")} }
      if(limit_by_depth_flag){ if (terra::crs(bath_global) != raster_crs_terra) { cat("DEBUG: Projecting global bathymetry...\n"); bath_global <- tryCatch(terra::project(bath_global, raster_crs_terra), error = function(e) { cat("ERROR: Failed project bathymetry:", e$message, "\n"); NULL }); if (is.null(bath_global)) { limit_by_depth_flag <- FALSE; cat("WARN: Disabling depth filter because projection failed.\n") } } }
      if (limit_by_depth_flag) {
        bath_species_extent <- tryCatch(terra::crop(bath_global, obis_calibration_vect, snap="near"), error = function(e) { cat("WARN: Failed cropping bathymetry for depth filter:", e$message, "\n"); NULL })
        if (is.null(bath_species_extent)) { cat("WARN: Proceeding without depth filter due to cropping error.\n"); limit_by_depth_flag <- FALSE; }
        else {
          occs_vect_proj <- tryCatch(terra::project(occs_vect, raster_crs_terra), error=function(e)NULL)
          if(is.null(occs_vect_proj)){ cat("WARN: Could not project occurrences to raster CRS for depth extraction. Skipping depth filter.\n"); limit_by_depth_flag <- FALSE; }
          else {
            bath_pts_vals <- tryCatch(terra::extract(bath_species_extent, occs_vect_proj, ID = FALSE), error=function(e){cat("WARN: Depth extraction failed:", e$message, "\n");NULL})
            if(is.null(bath_pts_vals) || nrow(bath_pts_vals) == 0 || all(is.na(bath_pts_vals[[1]]))) { cat("WARN: No valid depth values extracted for occurrences within ecoregion extent. Skipping depth filter.\n"); limit_by_depth_flag <- FALSE; }
            else {
              depth_range_pts <- range(bath_pts_vals[[1]], na.rm = TRUE); min_depth <- depth_range_pts[1] - config$depth_buffer_obis; max_depth <- depth_range_pts[2] + config$depth_buffer_obis; max_depth <- min(0, max_depth); cat("DEBUG:   Depth range filter:", min_depth, "to", max_depth, "m\n")
              depth_mask <- bath_species_extent; depth_mask[depth_mask < min_depth | depth_mask > max_depth] <- NA; depth_mask[!is.na(depth_mask)] <- 1
              max_val_check <- terra::global(depth_mask, "max", na.rm = TRUE)$max
              if (is.na(max_val_check) || max_val_check == 0) { cat("WARN: Depth filter masked all cells within ecoregion extent. Skipping depth filter.\n"); depth_mask <- NULL; limit_by_depth_flag <- FALSE; }
              else { cat("INFO: Depth filter mask created based on ecoregion extent and occurrence depths.\n") }
            }
          }
          if(exists("bath_species_extent_masked")) rm(bath_species_extent_masked); gc()
          if(exists("bath_pts_vals")) rm(bath_pts_vals); gc()
        }
        if(exists("bath_species_extent")) rm(bath_species_extent); gc()
      }
      if(exists("bath_global")) rm(bath_global); gc()
    }
  } else { cat("INFO: Depth filter disabled by config.\n") }
  
  # --- 3. Mask Environmental Layers ---
  cat("INFO: Masking global predictors with final species extent (Ecoregions +/- Depth)...\n")
  predictor_stack_obis_part2 <- global_predictor_stack
  
  if(!is.null(depth_mask) && limit_by_depth_flag){
    predictor_stack_obis_part2 <- tryCatch(terra::mask(terra::crop(predictor_stack_obis_part2, depth_mask, snap="near"), depth_mask),
                                           error = function(e){cat("ERROR: Failed depth masking predictor stack:", e$message, "\n"); NULL})
    if(is.null(predictor_stack_obis_part2)){ rm(depth_mask); gc(); return(NULL)}
    rm(depth_mask); gc()
    cat("DEBUG: Predictor stack masked by depth.\n")
  }
  
  predictor_stack_obis_part2 <- tryCatch(terra::mask(terra::crop(predictor_stack_obis_part2, obis_calibration_vect, snap="near"), obis_calibration_vect),
                                         error = function(e){cat("ERROR: Final ecoregion masking failed:", e$message, "\n"); NULL})
  if(is.null(predictor_stack_obis_part2)) return(NULL)
  
  min_max_check_final <- terra::global(predictor_stack_obis_part2[[1]], c("min", "max"), na.rm = TRUE)
  if (is.na(min_max_check_final$min) && is.na(min_max_check_final$max)) { cat("ERROR:   No valid cells after final masking.\n"); return(NULL) }
  cat("INFO:   Final species-specific stack created.\n")
  
  # --- 4. Calculate Quadrature Number (OBIS logic - mimicking .cm_calc_quad) ---
  cat("DEBUG: Calculating number of background points (quad_n) using OBIS logic...\n")
  
  # Use the target number from config as the base 'quad_samp' value
  # Their code uses a default of 50000 or 1% of available cells if quad_samp is between 0-1.
  # We'll use config$background_points_n as the target, but apply their adjustment logic.
  quad_samp_target <- config$background_points_n # Base target number
  n_presences <- nrow(occs_sf) # Number of presence points (after cleaning/thinning)
  
  # Count valid cells in the final masked stack
  env_size_t <- tryCatch({
    valid_cells <- terra::global(predictor_stack_obis_part2[[1]], fun="notNA") # Count non-NA cells
    valid_cells$notNA
  }, error=function(e){ cat("WARN: Could not count valid cells for quad_n calc.\n"); NA})
  
  if(is.na(env_size_t) || env_size_t == 0){ cat("ERROR: No valid cells in final masked stack to sample background from.\n"); return(NULL) }
  cat("DEBUG:   Number of valid cells in final masked stack:", env_size_t, "\n")
  
  # OBIS logic: if n_presences >= target background number, potentially double target,
  # BUT ensure total points (pres + bg) doesn't exceed available cells.
  quad_n_final <- quad_samp_target # Start with the config target
  
  if (n_presences >= quad_samp_target) {
    est_tsize <- n_presences * 2 # Potential doubling
    if (est_tsize > env_size_t) {
      # If doubled presences exceed available cells, cap bg points
      # And potentially reduce presences (OBIS does this, we won't here for simplicity in comparison)
      quad_n_final <- env_size_t - n_presences # Max possible background points
      if (quad_n_final < 0) quad_n_final <- 10 # Ensure at least some background points
      cat("WARN:   Doubled presences exceed available cells. Capping background points to:", quad_n_final, "(Total cells:", env_size_t, ", Presences:", n_presences, ")\n")
      # Note: OBIS code would also downsample presences here, we skip this.
    } else {
      # If doubling fits, use double the presences as the background target
      quad_n_final <- n_presences * 2
      cat("DEBUG:   Using doubled presences as background target:", quad_n_final, "\n")
    }
  }
  
  # Ensure quad_n doesn't exceed available cells
  quad_n_final <- min(quad_n_final, env_size_t)
  cat("INFO:   Final calculated number of background points (quad_n):", quad_n_final, "\n")
  
  # --- 5. Sample Background Points ---
  cat("INFO: Generating background points using OBIS Part2 logic & stack...\n")
  set.seed(config$background_seed %||% (111 + seed_offset)) # Use consistent seed logic if available
  
  # Sample directly from the final masked stack
  bg_obis_part2_df <- tryCatch({
    terra::spatSample(predictor_stack_obis_part2,
                      size = quad_n_final, # Use the calculated number
                      method = "random",
                      na.rm = TRUE,
                      xy = TRUE,
                      warn = FALSE)
  }, error = function(e) { cat("ERROR: spatSample failed for background points:", e$message, "\n"); NULL })
  
  if (is.null(bg_obis_part2_df)) { cat("ERROR:   Background point generation failed.\n"); return(NULL) }
  if (nrow(bg_obis_part2_df) < quad_n_final) { cat("WARN:   Could only sample", nrow(bg_obis_part2_df), "background points (requested", quad_n_final, ").\n") }
  if (nrow(bg_obis_part2_df) == 0) { cat("ERROR:   Failed to generate ANY background points.\n"); return(NULL) }
  
  # Keep only coordinates
  bg_coords_df <- as.data.frame(bg_obis_part2_df[, c("x", "y")])
  colnames(bg_coords_df) <- c("longitude", "latitude") # Match expected output format
  cat("INFO:   Generated", nrow(bg_coords_df), "OBIS Part2 background points.\n")
  
  # --- 6. Return Results ---
  cat("INFO: Returning background points, species stack, and polygon.\n")
  return(list(
    background_points = bg_coords_df,
    species_specific_stack = predictor_stack_obis_part2,
    calibration_polygon = obis_calibration_vect, # Polygon *before* depth masking
    quad_n_calculated = quad_n_final # Return the actual number sampled/calculated
  ))
}


generate_sdm_background_obis <- function(occs_sf, global_predictor_stack, config, logger, species_log_file = NULL, seed_offset = 4) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[GenerateSDMBackgroundOBISPart2]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  # --- Input Checks ---
  if (!inherits(occs_sf, "sf") || nrow(occs_sf) < 3) { hlog("ERROR", "Valid sf occurrence object with >= 3 points required.\n"); return(NULL) }
  if (is.na(sf::st_crs(occs_sf))) { hlog("ERROR", "Occurrence sf object needs a CRS.\n"); return(NULL) }
  if (!inherits(global_predictor_stack, "SpatRaster") || terra::nlyr(global_predictor_stack) < 1) { hlog("ERROR", "Valid global SpatRaster required.\n"); return(NULL) }
  if (terra::crs(global_predictor_stack) == "") { hlog("ERROR", "Global predictor stack needs a CRS.\n"); return(NULL) }
  if (!file.exists(config$ecoregion_shapefile)) { hlog("ERROR", "Ecoregion shapefile not found:", config$ecoregion_shapefile, "\n"); return(NULL) }
  if (config$limit_by_depth_obis && !file.exists(config$bathymetry_file)) { hlog("ERROR", "Bathymetry file required for depth filtering not found:", config$bathymetry_file, "\n"); return(NULL) }
  
  species_aphia_id_internal <- tryCatch(occs_sf$AphiaID[1], error = function(e) { hlog("WARN", "AphiaID missing from occs_sf, using random seed.\n"); sample.int(1e6, 1) })
  
  # --- 1. Ecoregion Extent Definition ---
  # ... (Keep the ecoregion loading, intersection, adjacency, buffering logic from the previous version) ...
  # --- Load Ecoregions ---
  hlog("DEBUG", "Loading ecoregions...\n")
  ecoregions <- tryCatch(terra::vect(config$ecoregion_shapefile), error = function(e) { hlog("ERROR", "Failed load ecoregions:", e$message, "\n"); NULL })
  if (is.null(ecoregions)) return(NULL)
  occs_vect <- tryCatch(terra::vect(occs_sf), error = function(e) { hlog("ERROR", "Failed convert occs sf to vect:", e$message, "\n"); NULL })
  if(is.null(occs_vect)) return(NULL)
  if (terra::crs(occs_vect) != terra::crs(ecoregions)) {
    hlog("DEBUG", "Projecting occurrences to match ecoregions CRS...\n")
    occs_vect <- tryCatch(terra::project(occs_vect, terra::crs(ecoregions)), error = function(e) { hlog("ERROR", "Failed project occurrences:", e$message, "\n"); NULL })
    if(is.null(occs_vect)) return(NULL)
  }
  # --- Identify Intersecting & Adjacent Ecoregions (OBIS Step) ---
  hlog("DEBUG", "Identifying intersecting ecoregions...\n")
  intersect_idx <- tryCatch(terra::is.related(ecoregions, occs_vect, "intersects"), error=function(e){hlog("WARN", "is.related failed:",e$message,"\n");NULL})
  if (is.null(intersect_idx) || !any(intersect_idx)) { hlog("WARN", "No occurrences intersect with ecoregions. Cannot define extent this way.\n"); return(NULL) }
  ecoreg_occ <- ecoregions[intersect_idx, ]
  hlog("DEBUG", "Finding adjacent ecoregions using buffered points...\n")
  sf::sf_use_s2(FALSE)
  occs_buffered_sf <- tryCatch(sf::st_buffer(sf::st_as_sf(occs_vect), dist = config$poly_buffer_obis), error = function(e) { hlog("ERROR", "Failed buffer points:", e$message, "\n"); NULL })
  if(is.null(occs_buffered_sf)) { sf::sf_use_s2(TRUE); return(NULL) }
  adjacent_idx <- tryCatch(terra::is.related(ecoregions, terra::vect(occs_buffered_sf), "intersects"), error=function(e){hlog("WARN", "is.related failed for adjacent:",e$message,"\n");NULL})
  if(is.null(adjacent_idx)) adjacent_idx <- intersect_idx # Fallback
  sf::sf_use_s2(TRUE)
  final_indices <- unique(c(which(intersect_idx), which(adjacent_idx)))
  ecoreg_sel <- ecoregions[final_indices, ]
  hlog("DEBUG", "Buffering final selected ecoregion polygon...\n")
  sf::sf_use_s2(FALSE)
  ecoreg_sel_buffered_sf <- tryCatch(sf::st_buffer(sf::st_as_sf(ecoreg_sel), dist = config$poly_buffer_final), error = function(e) { hlog("ERROR", "Failed buffer final polygon:", e$message, "\n"); NULL })
  if(is.null(ecoreg_sel_buffered_sf)) { sf::sf_use_s2(TRUE); return(NULL) }
  sf::sf_use_s2(TRUE)
  raster_crs_terra <- terra::crs(global_predictor_stack)
  obis_calibration_poly_sf <- sf::st_transform(ecoreg_sel_buffered_sf, crs = sf::st_crs(raster_crs_terra))
  obis_calibration_vect <- terra::vect(obis_calibration_poly_sf)
  hlog("INFO", "Ecoregion-based calibration polygon created (before depth filter).\n")
  
  # --- 2. Depth Filtering (Optional) ---
  # ... (Keep the depth filtering logic exactly as corrected in the previous response) ...
  depth_mask <- NULL
  limit_by_depth_flag <- config$limit_by_depth_obis
  if (limit_by_depth_flag) {
    hlog("INFO", "Applying OBIS depth filter...\n")
    bath_global <- tryCatch(terra::rast(config$bathymetry_file), error = function(e) { hlog("ERROR", "Failed load bathymetry:", e$message, "\n"); NULL })
    if (is.null(bath_global)) { limit_by_depth_flag <- FALSE; hlog("WARN", "Disabling depth filter because bathymetry loading failed.\n") }
    else {
      if (config$apply_indo_pacific_crop) { hlog("DEBUG", "Applying Indo-Pacific crop to global bathymetry layer...\n"); ip_extent <- terra::ext(config$indo_pacific_bbox); bath_global <- tryCatch(terra::crop(bath_global, ip_extent), error=function(e){ hlog("WARN", "Failed cropping global bathy:", e$message, "\n"); NULL}); if(is.null(bath_global)){ hlog("WARN", "Proceeding without depth filter due to bathy crop error.\n"); limit_by_depth_flag <- FALSE } else { hlog("DEBUG", "Global bathymetry cropped to Indo-Pacific extent.\n")} }
      if(limit_by_depth_flag){ if (terra::crs(bath_global) != raster_crs_terra) { hlog("DEBUG", "Projecting global bathymetry...\n"); bath_global <- tryCatch(terra::project(bath_global, raster_crs_terra), error = function(e) { hlog("ERROR", "Failed project bathymetry:", e$message, "\n"); NULL }); if (is.null(bath_global)) { limit_by_depth_flag <- FALSE; hlog("WARN", "Disabling depth filter because projection failed.\n") } } }
      if (limit_by_depth_flag) {
        bath_species_extent <- tryCatch(terra::crop(bath_global, obis_calibration_vect, snap="near"), error = function(e) { hlog("WARN", "Failed cropping bathymetry for depth filter:", e$message, "\n"); NULL })
        if (is.null(bath_species_extent)) { hlog("WARN", "Proceeding without depth filter due to cropping error.\n"); limit_by_depth_flag <- FALSE; }
        else {
          occs_vect_proj <- tryCatch(terra::project(occs_vect, raster_crs_terra), error=function(e)NULL)
          if(is.null(occs_vect_proj)){ hlog("WARN", "Could not project occurrences to raster CRS for depth extraction. Skipping depth filter.\n"); limit_by_depth_flag <- FALSE; }
          else {
            bath_pts_vals <- tryCatch(terra::extract(bath_species_extent, occs_vect_proj, ID = FALSE), error=function(e){hlog("WARN", "Depth extraction failed:", e$message, "\n");NULL})
            if(is.null(bath_pts_vals) || nrow(bath_pts_vals) == 0 || all(is.na(bath_pts_vals[[1]]))) { hlog("WARN", "No valid depth values extracted for occurrences within ecoregion extent. Skipping depth filter.\n"); limit_by_depth_flag <- FALSE; }
            else {
              depth_range_pts <- range(bath_pts_vals[[1]], na.rm = TRUE); min_depth <- depth_range_pts[1] - config$depth_buffer_obis; max_depth <- depth_range_pts[2] + config$depth_buffer_obis; max_depth <- min(0, max_depth); hlog("DEBUG", "  Depth range filter:", min_depth, "to", max_depth, "m\n")
              depth_mask <- bath_species_extent; depth_mask[depth_mask < min_depth | depth_mask > max_depth] <- NA; depth_mask[!is.na(depth_mask)] <- 1
              max_val_check <- terra::global(depth_mask, "max", na.rm = TRUE)$max
              if (is.na(max_val_check) || max_val_check == 0) { hlog("WARN", "Depth filter masked all cells within ecoregion extent. Skipping depth filter.\n"); depth_mask <- NULL; limit_by_depth_flag <- FALSE; }
              else { hlog("INFO", "Depth filter mask created based on ecoregion extent and occurrence depths.\n") }
            }
          }
          if(exists("bath_species_extent_masked")) rm(bath_species_extent_masked); gc()
          if(exists("bath_pts_vals")) rm(bath_pts_vals); gc()
        }
        if(exists("bath_species_extent")) rm(bath_species_extent); gc()
      }
      if(exists("bath_global")) rm(bath_global); gc()
    }
  } else { hlog("INFO", "Depth filter disabled by config.\n") }
  
  # --- 3. Mask Environmental Layers ---
  hlog("INFO", "Masking global predictors with final species extent (Ecoregions +/- Depth)...\n")
  predictor_stack_obis_part2 <- global_predictor_stack
  
  if(!is.null(depth_mask) && limit_by_depth_flag){
    predictor_stack_obis_part2 <- tryCatch(terra::mask(terra::crop(predictor_stack_obis_part2, depth_mask, snap="near"), depth_mask),
                                           error = function(e){hlog("ERROR", "Failed depth masking predictor stack:", e$message, "\n"); NULL})
    if(is.null(predictor_stack_obis_part2)){ rm(depth_mask); gc(); return(NULL)}
    rm(depth_mask); gc()
    hlog("DEBUG", "Predictor stack masked by depth.\n")
  }
  
  predictor_stack_obis_part2 <- tryCatch(terra::mask(terra::crop(predictor_stack_obis_part2, obis_calibration_vect, snap="near"), obis_calibration_vect),
                                         error = function(e){hlog("ERROR", "Final ecoregion masking failed:", e$message, "\n"); NULL})
  if(is.null(predictor_stack_obis_part2)) return(NULL)
  
  min_max_check_final <- terra::global(predictor_stack_obis_part2[[1]], c("min", "max"), na.rm = TRUE)
  if (is.na(min_max_check_final$min) && is.na(min_max_check_final$max)) { hlog("ERROR", "  No valid cells after final masking.\n"); return(NULL) }
  hlog("INFO", "  Final species-specific stack created.\n")
  
  # --- 4. Calculate Quadrature Number (OBIS logic - mimicking .cm_calc_quad) ---
  hlog("DEBUG", "Calculating number of background points (quad_n) using OBIS logic...\n")
  
  # Use the target number from config as the base 'quad_samp' value
  # Their code uses a default of 50000 or 1% of available cells if quad_samp is between 0-1.
  # We'll use config$background_points_n as the target, but apply their adjustment logic.
  quad_samp_target <- config$background_points_n # Base target number
  n_presences <- nrow(occs_sf) # Number of presence points (after cleaning/thinning)
  
  # Count valid cells in the final masked stack
  env_size_t <- tryCatch({
    valid_cells <- terra::global(predictor_stack_obis_part2[[1]], fun="notNA") # Count non-NA cells
    valid_cells$notNA
  }, error=function(e){ hlog("WARN", "Could not count valid cells for quad_n calc.\n"); NA})
  
  if(is.na(env_size_t) || env_size_t == 0){ hlog("ERROR", "No valid cells in final masked stack to sample background from.\n"); return(NULL) }
  hlog("DEBUG", "  Number of valid cells in final masked stack:", env_size_t, "\n")
  
  # OBIS logic: if n_presences >= target background number, potentially double target,
  # BUT ensure total points (pres + bg) doesn't exceed available cells.
  quad_n_final <- quad_samp_target # Start with the config target
  
  if (n_presences >= quad_samp_target) {
    est_tsize <- n_presences * 2 # Potential doubling
    if (est_tsize > env_size_t) {
      # If doubled presences exceed available cells, cap bg points
      # And potentially reduce presences (OBIS does this, we won't here for simplicity in comparison)
      quad_n_final <- env_size_t - n_presences # Max possible background points
      if (quad_n_final < 0) quad_n_final <- 10 # Ensure at least some background points
      hlog("WARN", "  Doubled presences exceed available cells. Capping background points to:", quad_n_final, "(Total cells:", env_size_t, ", Presences:", n_presences, ")\n")
      # Note: OBIS code would also downsample presences here, we skip this.
    } else {
      # If doubling fits, use double the presences as the background target
      quad_n_final <- n_presences * 2
      hlog("DEBUG", "  Using doubled presences as background target:", quad_n_final, "\n")
    }
  }
  
  # Ensure quad_n doesn't exceed available cells
  quad_n_final <- min(quad_n_final, env_size_t)
  hlog("INFO", "  Final calculated number of background points (quad_n):", quad_n_final, "\n")
  
  # --- 5. Sample Background Points ---
  hlog("INFO", "Generating background points using OBIS Part2 logic & stack...\n")
  set.seed(config$background_seed %||% (111 + seed_offset)) # Use consistent seed logic if available
  
  # Sample directly from the final masked stack
  bg_obis_part2_df <- tryCatch({
    terra::spatSample(predictor_stack_obis_part2,
                      size = quad_n_final, # Use the calculated number
                      method = "random",
                      na.rm = TRUE,
                      xy = TRUE,
                      warn = FALSE)
  }, error = function(e) { hlog("ERROR", "spatSample failed for background points:", e$message, "\n"); NULL })
  
  if (is.null(bg_obis_part2_df)) { hlog("ERROR", "  Background point generation failed.\n"); return(NULL) }
  if (nrow(bg_obis_part2_df) < quad_n_final) { hlog("WARN", "  Could only sample", nrow(bg_obis_part2_df), "background points (requested", quad_n_final, ").\n") }
  if (nrow(bg_obis_part2_df) == 0) { hlog("ERROR", "  Failed to generate ANY background points.\n"); return(NULL) }
  
  # Keep only coordinates
  bg_coords_df <- as.data.frame(bg_obis_part2_df[, c("x", "y")])
  
  hlog("INFO", "  Generated", nrow(bg_coords_df), "OBIS Part2 background points.\n")
  
  return(bg_coords_df)
}







#' Create Spatial CV Folds using blockCV (Simplified iobis/mpaeu_sdm Logic)
#'
#' Calculates block size via autocorrelation (if requested) and generates
#' spatial_grid folds using the 'size' argument. Optionally handles fixed range.
#' Directly returns the blockCV folds object ready for SDMtune::train.
#'
#' @param full_swd_data SWD object (output from SDMtune::prepareSWD).
#' @param predictor_stack SpatRaster stack (used for CRS).
#' @param config Project configuration list.
#' @param logger log4r logger object.
#' @param species_log_file Path to species-specific log file.
#' @return A blockCV folds object suitable for SDMtune::train, or NULL on error.
create_spatial_cv_folds_simplified <- function(full_swd_data, predictor_stack, config, logger, species_log_file) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[CreateCVFoldsSimp]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  # --- Extract Config Parameters ---
  cv_method      <- config$sdm_spatial_cv_type_to_use
  
  nfolds         <- config$sdm_n_folds         %||% 5
  auto_range     <- config$blockcv_auto_range  %||% TRUE
  range          <- config$blockcv_range_m     %||% 200000
  use_hexagon    <- config$blockcv_hexagon     %||% FALSE
  selection_type <- config$blockcv_selection   %||% "systematic"
  n_iterate      <- config$blockcv_n_iterate   %||% 300
  lat_blocks     <- config$blockcv_lat_blocks  %||% 50
  
  hlog("DEBUG", paste("Attempting to create", nfolds, "spatial folds (blockCV simplified logic)..."))
  
  # --- Prepare sf object from SWD data ---
  swd_coords_df <- as.data.frame(full_swd_data@coords); colnames(swd_coords_df) <- c("X", "Y")
  swd_coords_df$pa <- full_swd_data@pa
  target_crs_str <- terra::crs(predictor_stack, proj=TRUE); if (is.null(target_crs_str)|target_crs_str=="") target_crs_str <- "EPSG:4326"
  swd_sf <- tryCatch({ sf::st_as_sf(swd_coords_df, coords = c("X", "Y"), crs = target_crs_str) },
                     error = function(e){ hlog("ERROR", paste("Failed create sf from SWD:", e$message)); NULL })
  if(is.null(swd_sf)) return(NULL)
  hlog("DEBUG", paste("Created sf object with", nrow(swd_sf), "rows."))
  
  if (auto_range) {
    hlog("INFO", "Using auto range")
    sf::sf_use_s2(FALSE)
    auto_cor <- blockCV::cv_spatial_autocor(x = sf::st_make_valid(swd_sf), column = "pa")
    range = auto_cor$range
    hlog("INFO", paste("Auto range result: ", range))
  }
  
  
  if (cv_method == "spatial_grid") {
    hlog("INFO", "Generating grid blocks.")
    hlog("INFO", paste("Final range to use: ", range))
    spatial_folds <- blockCV::cv_spatial(r = predictor_stack, x = swd_sf, column = "pa", iteration = n_iterate, size = range,
                                         hexagon = use_hexagon, k = nfolds, progress = F, report = T, plot = F, selection = selection_type)
    blockCV::cv_plot(cv = spatial_folds, x = swd_sf, r = predictor_stack) + geom_sf(data = swd_sf, 
                                                                                    alpha = .5)
  }
  
  if (cv_method == "spatial_lat") {
    hlog("INFO", "Generating latitudinal blocks.")
    spatial_folds <- blockCV::cv_spatial(r = predictor_stack, x = swd_sf, column = "pa", iteration = n_iterate, rows_cols = c(lat_blocks, 1),
                                         hexagon = use_hexagon, k = nfolds, progress = F, report = T, plot = F, extend = 0.5)
    blockCV::cv_plot(cv = spatial_folds, x = swd_sf) + geom_sf(data = swd_sf, 
                                                               alpha = .5)
  }
  
  if (cv_method == "random") {
    hlog("INFO", "Generating random samples.")
    spatial_folds <- SDMtune::randomFolds(size = range, k = nfolds, only_presence = TRUE)
    # spatial_folds <- blockCV::cv_spatial(swd_sf, column = "pa", iteration = n_iterate, size = range,
    #                                   hexagon = use_hexagon, k = nfolds, progress = F, report = T, plot = F, selection = selection_type)
  }
  
  if(is.null(spatial_folds)) { hlog("ERROR", "blockCV::cv_spatial returned NULL."); return(NULL) }
  
  hlog("DEBUG", paste("Spatial folds object created using blockCV (k=", nfolds, ")."))
  
  return(spatial_folds) # Return the blockCV object directly
}


#' Tune SDM Hyperparameters using SDMtune gridSearch with SPATIAL k-fold CV
#' (Wrapper calling simplified blockCV helper)
run_sdm_tuning_scv <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species", species_log_file = NULL) { # Name kept as requested
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SCV_TuneHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  k_folds_used <- config$sdm_n_folds
  hlog("INFO", paste("Tuning hyperparameters for", species_name, "using SPATIAL", k_folds_used, "-fold CV (blockCV simplified)..."))
  
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df)) { hlog("ERROR", "Invalid inputs for tuning."); return(NULL) }
  
  # --- 1. Prepare SWD object (handles NA removal) ---
  full_swd_data <- tryCatch({
    SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE)
  }, error = function(e) { hlog("ERROR", paste("Failed prepareSWD:", e$message)); return(NULL) })
  if (is.null(full_swd_data)) return(NULL)
  hlog("DEBUG", paste("SWD object prepared. Rows:", nrow(full_swd_data@coords)))
  
  # --- 2. Create Spatial Folds using the new simplified helper ---
  spatial_folds_blockcv <- create_spatial_cv_folds_simplified(
    full_swd_data,
    predictor_stack, # Pass stack for CRS info
    config,
    logger,
    species_log_file
  )
  if(is.null(spatial_folds_blockcv)) { hlog("ERROR", "Failed to create spatial folds."); return(NULL) }
  
  # --- 3. Train initial CV model ---
  hlog("DEBUG", "Training initial CV base model for gridSearch...")
  initial_cv_model <- tryCatch({
    SDMtune::train(method = config$sdm_method, data = full_swd_data, folds = spatial_folds_blockcv, progress = FALSE)
  }, error = function(e) {
    if (grepl("glmnet failed", e$message, ignore.case = TRUE)) {
      hlog("ERROR", paste("Failed train CV base model:", e$message))
      hlog("ERROR", "  ---> glmnet error occurred. Data issues within folds possible or block size might be problematic. Check data or adjust blockcv_auto_range/blockcv_range_m.")
    } else { hlog("ERROR", paste("Failed train CV base model:", e$message)) }
    return(NULL)
  })
  if (is.null(initial_cv_model) || !inherits(initial_cv_model, "SDMmodelCV")) { hlog("ERROR", "Failed create valid SDMmodelCV object."); return(NULL)}
  
  # --- 4. Hyperparameter Tuning ---
  # ... (Rest of the function remains the same) ...
  hyper_grid <- config$sdm_tune_grid; tuning_results <- NULL
  tryCatch({
    hlog("DEBUG", "Running SDMtune::gridSearch with spatial k-fold CV...")
    show_progress <- FALSE
    if (!is.null(logger)) { tryCatch({ current_log_level_num <- log4r::log_level(config$log_level %||% "INFO"); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num }, error=function(e){}) }
    
    tuning_results <- SDMtune::gridSearch(model = initial_cv_model, hypers = hyper_grid, metric = config$sdm_evaluation_metric, save_models = TRUE, progress = show_progress, interactive = FALSE)
    hlog("DEBUG", "SDMtune::gridSearch completed.")
    
    res_df <- tuning_results@results
    if(is.null(res_df) || nrow(res_df) == 0) { hlog("WARN", "No results found in tuning object."); return(NULL) }
    
    # Result selection logic
    metric_base_upper <- toupper(config$sdm_evaluation_metric)
    target_metric_col <- if (tolower(config$sdm_evaluation_metric) == "aicc") "AICc" else paste0("test_", metric_base_upper)
    if (!target_metric_col %in% colnames(res_df)) { hlog("ERROR", paste("CV metric '", target_metric_col, "' not found. Available:", paste(colnames(res_df), collapse=", "))); print(head(res_df)); return(NULL) }
    valid_metric_indices <- which(!is.na(res_df[[target_metric_col]]))
    if(length(valid_metric_indices) == 0) { hlog("WARN", paste("All values for metric '", target_metric_col, "' are NA.")); return(NULL) }
    select_fun <- if (target_metric_col == "AICc") which.min else which.max
    best_row_relative_index <- select_fun(res_df[[target_metric_col]][valid_metric_indices])
    best_row_index <- valid_metric_indices[best_row_relative_index]
    if(length(best_row_index) == 0 || best_row_index > nrow(res_df)) { hlog("WARN", "Could not determine best model index, using first valid row."); best_row_index <- valid_metric_indices[1]}
    best_hypers_df <- res_df[best_row_index, names(hyper_grid), drop = FALSE]
    hlog("INFO", paste("  Best hypers (Mean CV", target_metric_col,"=", round(res_df[[target_metric_col]][best_row_index], 4),"): ", paste(names(best_hypers_df), best_hypers_df[1,], collapse=", ")))
    attr(tuning_results, "best_hypers") <- best_hypers_df
    return(tuning_results)
    
  }, error = function(e) { hlog("ERROR", paste("SDMtune::gridSearch failed:", e$message)); return(NULL) })
}




#' Train Final SDM Model
#' Returns the trained model object.
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to a species-specific log file.
#' @return An `SDMmodel` object, or NULL on error.
train_final_sdm <- function(occs_coords, predictor_stack, background_df, best_hypers, config, logger, species_name = "species", species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[TrainHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  hlog("INFO", paste("Training final model for", species_name, "on full dataset..."))
  if (is.null(occs_coords) || nrow(occs_coords) < config$min_occurrences_sdm || is.null(predictor_stack) || is.null(background_df) || is.null(best_hypers)) { hlog("ERROR", "Invalid inputs."); return(NULL) }
  
  full_swd_data <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_df, env = predictor_stack, verbose = FALSE)}, error = function(e) { hlog("ERROR", paste("Failed prepareSWD:", e$message)); return(NULL)})
  if (is.null(full_swd_data)) return(NULL)
  
  train_args <- list(method = config$sdm_method, data = full_swd_data, progress = FALSE)
  for (h_name in names(best_hypers)) { h_value <- best_hypers[[h_name]]; if (h_name == "reg") h_value <- as.numeric(h_value); train_args[[h_name]] <- h_value }
  hlog("DEBUG", paste("  Using final hyperparameters:", paste(names(train_args)[-(1:2)], sapply(train_args[-(1:2)], as.character), collapse=", ")))
  
  final_model <- tryCatch({ do.call(SDMtune::train, train_args)}, error = function(e) { hlog("ERROR", paste("Failed train final model:", e$message)); return(NULL)})
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Final model training returned invalid object."); return(NULL) }
  
  hlog("INFO", paste("  Final model trained successfully for", species_name))
  return(final_model)
}


#' Construct Target Prediction Filename
#' Builds the expected filename and path for a prediction raster based on target structure.
#' Used for checking if a prediction file already exists.
#'
#' @param species_name_sanitized Sanitized species name (e.g., "Genus_species").
#' @param scenario_name The name of the scenario (e.g., "current", "ssp119_2050").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_combined_pca").
#' @param config The configuration list.
#' @return Character string with the full path to the expected prediction file.
construct_prediction_filename <- function(species_name_sanitized, scenario_name, predictor_type_suffix, config) {
  
  # Construct base filename part (e.g., mean_pred_SPECIES)
  # Using sanitized name for now.
  base_filename <- paste0("mean_pred_", species_name_sanitized)
  
  target_dir <- NULL
  target_filename_stem <- NULL # Filename without extension
  
  if (scenario_name == "current") {
    target_dir <- config$target_predictions_current_dir
    # Filename format: mean_pred_SPECIES[_MODELTYPE?]
    target_filename_stem <- paste0(base_filename, predictor_type_suffix)
  } else {
    # Future scenario
    ssp_dir_name <- config$ssp_scenario_map[[scenario_name]]
    if (is.null(ssp_dir_name)) {
      warning("No SSP directory mapping for scenario: ", scenario_name, call. = FALSE)
      return(NULL) # Return NULL if mapping is missing
    }
    target_dir <- file.path(config$target_predictions_future_base, ssp_dir_name)
    
    # Extract time tag for filename
    time_tag_clean <- ifelse(grepl("2050$", scenario_name), "dec50", ifelse(grepl("2100$", scenario_name), "dec100", "UNKNOWN"))
    if(time_tag_clean == "UNKNOWN") {
      warning("Cannot extract time tag from scenario: ", scenario_name, call. = FALSE)
      return(NULL) # Return NULL if time tag is unknown
    }
    
    # Filename format: mean_pred_SPECIES_SSP_TIME[_MODELTYPE?]
    target_filename_stem <- paste0(base_filename, "_", ssp_dir_name, "_", time_tag_clean, predictor_type_suffix)
  }
  
  if(is.null(target_dir) || is.null(target_filename_stem)) {
    warning("Could not construct target path components.", call. = FALSE)
    return(NULL)
  }
  
  # Add extension
  target_filename <- paste0(target_filename_stem, ".tif")
  
  return(file.path(target_dir, target_filename))
}


#' Predict SDM Suitability (Returns Raster or Error Message)
#' @return A SpatRaster object with the prediction, or a character string with the error message on failure.
predict_sdm_suitability <- function(final_sdm_model, predictor_stack, config, logger, species_log_file = NULL, output_type = "cloglog") {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[PredictHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  hlog("DEBUG", paste("Attempting prediction (type:", output_type, ")..."))
  if (is.null(final_sdm_model) || !inherits(final_sdm_model, "SDMmodel")) { msg <- "Invalid final_sdm_model."; hlog("ERROR", msg); return(msg) }
  if (is.null(predictor_stack)) { msg <- "Predictor stack required."; hlog("ERROR", msg); return(msg) }
  
  prediction_result <- NULL
  tryCatch({
    prediction_result <- SDMtune::predict(object = final_sdm_model, data = predictor_stack, type = output_type, clamp = TRUE)
    if (is.null(prediction_result) || !inherits(prediction_result, "SpatRaster") || terra::nlyr(prediction_result) == 0) { msg <- "Prediction returned NULL or empty raster."; hlog("WARN", msg); return(msg) }
    names(prediction_result) <- "suitability"; hlog("DEBUG", "Prediction raster generated."); return(prediction_result)
  }, error = function(e) { err_msg <- paste("SDMtune::predict failed:", e$message); hlog("ERROR", err_msg); return(err_msg) })
}

#' Save Tuning Results (RDS and CSV)
#' Saves the full tuning object and the results table separately.
#' Constructs paths based on target structure in config.
#'
#' @param tuning_output The object returned by `run_sdm_tuning_kfold` (contains `@results`).
#' @param species_name_sanitized Sanitized species name (e.g., "Genus_species").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_combined_pca").
#' @param config The configuration list (needs target_results_base, model_output_subdir_map).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_tuning_results <- function(tuning_output, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveTuneHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  if (is.null(tuning_output) || !inherits(tuning_output, "SDMtune")) { hlog("ERROR", "Invalid tuning object provided."); return(FALSE) }
  
  # Determine target subdirectory
  subdir_name <- config$model_output_subdir_map[[predictor_type_suffix]]
  if (is.null(subdir_name)) { hlog("ERROR", paste("No output subdirectory mapping found for suffix:", predictor_type_suffix)); return(FALSE) }
  target_subdir <- file.path(config$target_results_base, subdir_name)
  dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
  
  # Construct filenames matching the target repo style (using original species short codes if available, else sanitized name)
  # NOTE: Using sanitized name for now, post-analysis script will need to adapt.
  target_base_name <- paste0("CV_Results_", species_name_sanitized)
  
  # 1. Save full tuning object (RDS - for potential internal reuse)
  rds_file <- file.path(target_subdir, paste0(target_base_name, "_tuning_object.rds")) # Different name to avoid clash
  tryCatch({ saveRDS(tuning_output, file = rds_file); hlog("DEBUG", "Full tuning object saved (RDS):", basename(rds_file)) },
           error = function(e) { hlog("ERROR", paste("Failed save tuning RDS:", e$message)) }) # Log error but continue
  
  # 2. Save results table (CSV - for analysis mirroring target)
  csv_file <- file.path(target_subdir, paste0(target_base_name, ".csv"))
  results_df <- tuning_output@results
  if (is.null(results_df) || nrow(results_df) == 0) { hlog("WARN", "No results table found in tuning object."); return(TRUE) } # Return TRUE as RDS might have saved
  
  tryCatch({ readr::write_csv(results_df, file = csv_file); hlog("DEBUG", "Tuning results table saved (CSV):", basename(csv_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save tuning CSV:", e$message)); FALSE })
}


#' Save Final SDM Model Object
#' Saves the trained SDMmodel object as an RDS file.
#' NOTE: This saves to an *intermediate* location (`config$models_dir`), not the target analysis structure.
#'
#' @param final_model The trained `SDMmodel` object.
#' @param species_name_sanitized Sanitized species name.
#' @param predictor_type_suffix Suffix indicating model type.
#' @param config The configuration list (needs `models_dir`).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_final_model <- function(final_model, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveModelHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Invalid model object provided."); return(FALSE) }
  
  # Create model type subdirectory if needed
  model_subdir <- file.path(config$models_dir, paste0(basename(config$anemone_occurrence_dir), predictor_type_suffix)) # Example based on group occurrence dir name
  dir.create(model_subdir, recursive = TRUE, showWarnings = FALSE)
  
  model_file <- file.path(model_subdir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  tryCatch({ saveRDS(final_model, file = model_file); hlog("DEBUG", "Final model object saved (RDS):", basename(model_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save final model RDS:", e$message)); FALSE })
}


#' Save SDM Prediction Raster
#' Saves the prediction SpatRaster to the target output directory structure.
#' Handles current vs future scenario paths.
#'
#' @param prediction_raster The predicted SpatRaster object.
#' @param species_name_sanitized Sanitized species name.
#' @param scenario_name The name of the scenario (e.g., "current", "ssp119_2050").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_combined_pca").
#' @param config The configuration list (needs target prediction dirs, ssp_scenario_map).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_sdm_prediction <- function(prediction_raster, species_name_sanitized, scenario_name, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SavePredHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}#cat(msg, "\n")} }
  
  if (is.null(prediction_raster) || !inherits(prediction_raster, "SpatRaster")) { hlog("ERROR", "Invalid prediction raster provided."); return(FALSE) }
  
  # Determine target directory and filename structure
  target_dir <- NULL
  target_filename <- NULL
  
  # Construct base filename part (e.g., mean_pred_SPECIES)
  # NOTE: Using sanitized name. Adjust post-analysis or add mapping if target uses codes.
  base_filename <- paste0("mean_pred_", species_name_sanitized) # Using SDMtune mean convention implicitly
  
  if (scenario_name == "current") {
    target_dir <- config$target_predictions_current_dir
    # Filename format: mean_pred_SPECIES[_MODELTYPE?].tif
    target_filename <- paste0(base_filename, predictor_type_suffix, ".tif") # Add suffix to distinguish model types
  } else {
    # Future scenario
    ssp_dir_name <- config$ssp_scenario_map[[scenario_name]]
    if (is.null(ssp_dir_name)) { hlog("ERROR", paste("No SSP directory mapping found for scenario:", scenario_name)); return(FALSE) }
    target_dir <- file.path(config$target_predictions_future_base, ssp_dir_name)
    
    # Extract time tag for filename
    time_tag_clean <- ifelse(grepl("2050$", scenario_name), "dec50", ifelse(grepl("2100$", scenario_name), "dec100", "UNKNOWN"))
    if(time_tag_clean == "UNKNOWN") { hlog("ERROR", paste("Cannot extract time tag from scenario:", scenario_name)); return(FALSE) }
    
    # Filename format: mean_pred_SPECIES_SSP_TIME[_MODELTYPE?].tif
    target_filename <- paste0(base_filename, "_", ssp_dir_name, "_", time_tag_clean, predictor_type_suffix, ".tif") # Add suffix
  }
  
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  pred_file_path <- file.path(target_dir, target_filename)
  
  tryCatch({
    terra::writeRaster(prediction_raster, filename = pred_file_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    hlog("DEBUG", paste("Prediction raster saved:", basename(pred_file_path), "to", target_dir))
    TRUE
  }, error = function(e) {
    hlog("ERROR", paste("Failed save prediction raster:", e$message))
    FALSE
  })
}


# scripts/helpers/sdm_modeling_helpers.R

# ... (keep other functions above) ...

#' Calculate and Save Variable Importance (v9 - Handles MaxNet Directly)
#' Calculates permutation importance. For MaxNet models, it extracts the
#' importance calculated during training. For other models, it uses SDMtune::varImp.
#' Saves results to the target structure.
#'
#' @param final_model An `SDMmodel` object (output from `train_final_sdm`).
#' @param training_swd An `SWD` object (potentially needed for varImp on non-MaxNet).
#' @param species_name_sanitized Sanitized species name.
#' @param group_name Group name ("anemone" or "anemonefish").
#' @param predictor_type_suffix Suffix indicating model type. Used in filename.
#' @param config The configuration list.
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
#' @export
calculate_and_save_vi <- function(final_model, training_swd, # training_swd kept for potential use by other methods
                                  species_name_sanitized, group_name, predictor_type_suffix,
                                  config, logger, species_log_file = NULL) {
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[VarImpHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  # --- Input validation ---
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) {
    hlog("ERROR", "Invalid SDMmodel object provided for VI.")
    return(FALSE)
  }
  # Training SWD validation might only be strictly needed if varImp is called below
  # if (is.null(training_swd) || !inherits(training_swd, "SWD")) {
  #   hlog("ERROR", "Training SWD object required for varImp calculation.")
  #   return(FALSE)
  # }
  
  hlog("INFO", "Calculating/Extracting variable importance...")
  vi_results_df <- NULL
  
  tryCatch({
    model_method <- final_model@method # Get the method used (e.g., "Maxnet")
    
    if (model_method == "Maxnet") {
      hlog("DEBUG", "Maxnet model detected. Extracting pre-calculated permutation importance.")
      # Access the raw maxnet model object stored by SDMtune
      raw_maxnet_model <- final_model@model
      if (!is.null(raw_maxnet_model) && !is.null(raw_maxnet_model$variable.importance)) {
        # The importance is stored as a named numeric vector
        importance_vector <- raw_maxnet_model$variable.importance
        # Convert to a standard data frame
        vi_results_df <- data.frame(
          Variable = names(importance_vector),
          Importance = as.numeric(importance_vector),
          stringsAsFactors = FALSE
        )
        vi_results_df <- vi_results_df[order(-vi_results_df$Importance), ] # Order descending
        hlog("DEBUG", "Successfully extracted MaxNet variable importance.")
      } else {
        hlog("WARN", "Could not find pre-calculated variable importance in the MaxNet model object.")
      }
    } else {
      # Fallback to SDMtune::varImp for other methods (e.g., RF, BRT) - may need testing
      hlog("DEBUG", paste("Model method is", model_method, ". Attempting SDMtune::varImp..."))
      if (is.null(training_swd) || !inherits(training_swd, "SWD")) {
        hlog("ERROR", "Training SWD object required for varImp calculation with non-MaxNet models.")
        return(FALSE)
      }
      show_progress <- FALSE
      if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::log_level(config$log_level %||% "INFO"); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE}) }
      
      vi_results <- SDMtune::varImp(
        model = final_model,
        permut = 10,
        progress = show_progress,
        test = training_swd # Provide test data for permutation
      )
      if (!is.null(vi_results) && nrow(vi_results) > 0) {
        vi_results_df <- vi_results # Already a data frame
        hlog("DEBUG", "SDMtune::varImp successful for non-MaxNet model.")
      } else {
        hlog("WARN", "SDMtune::varImp returned empty results for non-MaxNet model.")
      }
    }
    
    # --- Saving Logic ---
    if (is.null(vi_results_df)) {
      hlog("ERROR", "Variable importance could not be obtained.")
      return(FALSE)
    }
    
    vi_subdir_name <- paste0("vi_", group_name)
    # Use the correct target base directory from config
    target_subdir <- file.path(config$target_vi_base, vi_subdir_name) # Use config$target_vi_base
    dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
    vi_filename <- paste0("vi_", species_name_sanitized, predictor_type_suffix, ".csv")
    vi_file_path <- file.path(target_subdir, vi_filename)
    
    readr::write_csv(vi_results_df, vi_file_path)
    hlog("INFO", paste("Variable importance saved to:", vi_file_path))
    return(TRUE)
    
  }, error = function(e) {
    hlog("ERROR", paste("Variable importance processing/saving failed:", e$message))
    # print(rlang::trace_back()) # Uncomment for deeper debugging
    return(FALSE)
  })
}

#-------------------------------------------------------------------------