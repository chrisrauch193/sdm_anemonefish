# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow v2)
# - Saving logic moved to dedicated helper functions.
# - Path construction uses new config structure for target output.
# - Added Variable Importance helper.
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr, blockCV, ggplot2, ecospat, spThin, geosphere) # Added readr

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





#' Thin Occurrence Records Based on Spatial Autocorrelation (Mantel Correlogram)
#'
#' Calculates spatial autocorrelation using Mantel correlograms and optionally
#' thins occurrence points to the minimum non-significant distance found.
#'
#' @param occs_coords Data frame or matrix with occurrence coordinates (colnames 'longitude', 'latitude').
#' @param predictor_stack SpatRaster stack covering the occurrence area. Used for extracting env data.
#' @param config Configuration list with SAC parameters.
#' @param logger Logger object.
#' @param species_log_file Path to species log.
#' @return A list containing:
#'         `coords_thinned`: Data frame of thinned coordinates or original coordinates if thinning not applied/failed.
#'         `n_thinned`: Number of points after thinning.
#'         `thinning_distance_km`: The distance (km) used for thinning (or NA if not thinned).
#'         or NULL if a fatal error occurs.
#' @export
thin_occurrences_by_sac <- function(occs_coords, predictor_stack, config, logger, species_log_file) {
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SacThinHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  if (!config$apply_sac_thinning) {
    hlog("INFO", "Spatial autocorrelation thinning disabled in config.")
    return(list(coords_thinned = as.data.frame(occs_coords),
                n_thinned = nrow(occs_coords),
                thinning_distance_km = NA))
  }
  
  if (!requireNamespace("ecospat", quietly = TRUE) || !requireNamespace("spThin", quietly = TRUE)) {
    hlog("ERROR", "Packages 'ecospat' and 'spThin' required for SAC thinning.")
    return(NULL)
  }
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    hlog("ERROR", "Package 'geosphere' required for distance calculations.")
    return(NULL)
  }
  
  
  hlog("INFO", "Applying spatial autocorrelation thinning (Mantel method)...")
  
  # --- Prepare Data ---
  pts_df <- as.data.frame(occs_coords)
  colnames(pts_df) <- c("longitude", "latitude") # Ensure correct names
  
  # Limit points for correlogram calculation if dataset is very large
  if (nrow(pts_df) > 1000) {
    hlog("DEBUG", "Subsampling to 1000 points for Mantel correlogram calculation.")
    set.seed(123) # for reproducibility of subsampling
    pts_for_mantel <- pts_df[sample(1:nrow(pts_df), 1000, replace = FALSE), ]
  } else {
    pts_for_mantel <- pts_df
  }
  pts_for_mantel <- pts_for_mantel[!duplicated(pts_for_mantel), ] # Remove duplicates just in case
  
  # Extract environmental data at points
  env_extract <- tryCatch({
    terra::extract(predictor_stack, pts_for_mantel, ID = FALSE, na.rm = FALSE) # Keep NAs for now
  }, error = function(e){
    hlog("ERROR", paste("Failed extracting env data for SAC:", e$message)); NULL
  })
  if(is.null(env_extract)) return(NULL)
  
  # Combine coords and env data, remove rows with any NAs in env data
  data_df_mantel <- cbind(pts_for_mantel, env_extract)
  valid_rows_mantel <- complete.cases(data_df_mantel)
  data_df_mantel <- data_df_mantel[valid_rows_mantel, ]
  
  if (nrow(data_df_mantel) < 10) { # Need sufficient points for calculation
    hlog("WARN", "Less than 10 valid points with env data for Mantel test. Skipping thinning.")
    return(list(coords_thinned = pts_df, n_thinned = nrow(pts_df), thinning_distance_km = NA))
  }
  hlog("DEBUG", paste("Using", nrow(data_df_mantel), "points for Mantel correlogram."))
  
  # --- Calculate Mantel Correlogram ---
  # ecospat requires matrix/dataframe input
  coords_for_mantel <- data_df_mantel[, c("longitude", "latitude")]
  env_for_mantel <- data_df_mantel[, -(1:2)] # Exclude lon/lat columns
  
  # Calculate distances in meters using geosphere for accuracy with lon/lat
  # Note: ecospat often works with km, but let's calculate in m and convert class dist
  geo_dist_m <- tryCatch(geosphere::distm(coords_for_mantel), error=function(e){
    hlog("ERROR", paste("Distance calculation failed:", e$message)); NULL
  })
  if(is.null(geo_dist_m)) return(NULL)
  
  env_dist <- tryCatch(dist(scale(env_for_mantel)), error=function(e){ # Scale env vars
    hlog("ERROR", paste("Env distance calculation failed:", e$message)); NULL
  })
  if(is.null(env_dist)) return(NULL)
  
  # Convert classdist and maxdist from config (assumed meters) to km if needed by ecospat, or keep as is
  # ecospat.mantel.correlogram expects max distance, number of classes
  # Let's use meters internally and convert at the end.
  class_dist_m <- config$autocor_classdist
  max_dist_m <- config$autocor_maxdist
  n_classes <- floor(max_dist_m / class_dist_m)
  
  if(n_classes < 2) {
    hlog("WARN", "Max distance or class distance results in < 2 classes. Skipping thinning.")
    return(list(coords_thinned = pts_df, n_thinned = nrow(pts_df), thinning_distance_km = NA))
  }
  
  mantel_res <- NULL
  tryCatch({
    mantel_res <- ecospat::ecospat.mantel.correlogram(
      dfvar = as.data.frame(scale(env_for_mantel)), # Needs dataframe
      colxy = coords_for_mantel,
      nclass = n_classes,
      max = max_dist_m, # Provide max distance in meters
      nperm = 100 # Number of permutations for significance
    )
  }, error = function(e) {
    hlog("ERROR", paste("ecospat.mantel.correlogram failed:", e$message))
    # Provide more detail if possible
    if(grepl("smaller than the radius of the neighbourhood", e$message)) {
      hlog("ERROR", "  -> Hint: Check if autocor_classdist is large relative to data spread, or if points are too clustered.")
    }
    NULL
  })
  
  if (is.null(mantel_res) || is.null(mantel_res$mantel.res)) {
    hlog("WARN", "Mantel correlogram calculation failed or returned no results. Skipping thinning.")
    return(list(coords_thinned = pts_df, n_thinned = nrow(pts_df), thinning_distance_km = NA))
  }
  
  # Find first non-significant distance class (upper bound)
  # ecospat output: mantel.res has columns: n.tests, Mantel.cor, p.value, p.val.adj
  results_df_mantel <- as.data.frame(mantel_res$mantel.res)
  results_df_mantel$dist.class.m <- mantel_res$breaks[-1] # Upper bound of distance class in meters
  
  first_nonsig_idx <- which(results_df_mantel$p.value > config$autocor_signif)[1]
  
  distance_m <- NA
  if (!is.na(first_nonsig_idx)) {
    distance_m <- results_df_mantel$dist.class.m[first_nonsig_idx]
    hlog("INFO", paste("First non-significant distance class ends at:", round(distance_m / 1000, 1), "km"))
  } else {
    distance_m <- max_dist_m # If all are significant, use max distance tested
    hlog("WARN", paste("All distance classes up to", round(max_dist_m / 1000, 1), "km showed significant autocorrelation. Using max distance for potential thinning."))
  }
  
  # --- Apply Thinning using spThin ---
  thinned_coords <- pts_df # Default to original if thinning not applied/needed
  n_final <- nrow(thinned_coords)
  thin_dist_km_final <- NA
  
  prune_thresh_m <- config$sac_prune_threshold
  if (distance_m >= prune_thresh_m) {
    hlog("INFO", paste("Non-significant distance (", round(distance_m / 1000, 1), "km) meets threshold (", round(prune_thresh_m / 1000, 1), "km). Applying thinning..."))
    thin_dist_km_final <- round(distance_m / 1000, 1) # Distance used for thinning
    
    # spThin expects specific column names and a species column
    pts_for_thinning <- pts_df
    colnames(pts_for_thinning) <- c("LONG", "LAT") # Match spThin defaults if possible
    pts_for_thinning$SPECIES <- "MySpecies"
    
    # Check number of points vs spThin limit
    if(nrow(pts_for_thinning) > 15000) {
      hlog("WARN", "More than 15000 points, spThin might be slow or error-prone. Consider grid thinning as alternative if needed.")
      # Add grid thinning alternative here if desired
    }
    
    thinned_list <- NULL
    tryCatch({
      # Run thinning
      thinned_list <- spThin::thin(
        loc.data = pts_for_thinning,
        lat.col = "LAT",
        long.col = "LONG",
        spec.col = "SPECIES",
        thin.par = thin_dist_km_final, # Use distance in KM
        reps = 1, # Only need one rep
        locs.thinned.list.return = TRUE,
        write.files = FALSE,
        write.log.file = FALSE,
        verbose = FALSE # Control verbosity internally
      )
    }, error = function(e) {
      hlog("ERROR", paste("spThin::thin failed:", e$message))
      NULL # Return NULL on error
    })
    
    if (!is.null(thinned_list) && length(thinned_list) > 0 && inherits(thinned_list[[1]], "data.frame")) {
      thinned_coords <- thinned_list[[1]]
      colnames(thinned_coords) <- c("longitude", "latitude") # Revert names
      n_final <- nrow(thinned_coords)
      hlog("INFO", paste("Thinning resulted in", n_final, "records (using", thin_dist_km_final, "km distance)."))
    } else {
      hlog("WARN", "spThin failed or returned empty list. Proceeding with unthinned data.")
      thin_dist_km_final <- NA # Reset distance if thinning failed
    }
    # Remove spThin log if created
    if(file.exists("spatial_thin_log.txt")) file.remove("spatial_thin_log.txt")
    
  } else {
    hlog("INFO", paste("First non-significant distance (", round(distance_m / 1000, 1), "km) is below threshold (", round(prune_thresh_m / 1000, 1), "km). No thinning applied."))
    thinned_coords <- pts_df # Keep original points
    n_final <- nrow(thinned_coords)
    thin_dist_km_final <- NA
  }
  
  return(list(coords_thinned = thinned_coords,
              n_thinned = n_final,
              thinning_distance_km = thin_dist_km_final))
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
  
  # --- 2. Depth Filtering (Using Fixed Limits from Config) ---
  depth_mask <- NULL
  limit_by_depth_flag <- config$limit_by_depth_obis # Check if depth filtering is enabled at all
  
  if (limit_by_depth_flag) {
    hlog("INFO", "Applying FIXED depth filter based on config$depth_min and config$depth_max...")
    
    # Check if fixed limits are valid
    min_depth_fixed <- config$depth_min
    max_depth_fixed <- config$depth_max
    if (is.null(min_depth_fixed) || is.null(max_depth_fixed) || !is.numeric(min_depth_fixed) || !is.numeric(max_depth_fixed)) {
      hlog("ERROR", "config$depth_min or config$depth_max not found or not numeric. Cannot apply fixed depth filter.")
      limit_by_depth_flag <- FALSE # Disable filtering if config values are bad
    } else {
      hlog("INFO", paste("  Using fixed depth range:", min_depth_fixed, "to", max_depth_fixed, "m"))
      
      # --- Load and prepare bathymetry raster ---
      bath_global <- tryCatch(terra::rast(config$bathymetry_file), error = function(e) {
        hlog("ERROR", paste("Failed load bathymetry:", e$message))
        NULL
      })
      
      if (is.null(bath_global)) {
        limit_by_depth_flag <- FALSE
        hlog("WARN", "Disabling depth filter because bathymetry loading failed.")
      } else {
        # Optional: Crop global bathy to Indo-Pacific (consistent with PCA prep)
        if (config$apply_indo_pacific_crop) {
          hlog("DEBUG", "Applying Indo-Pacific crop to global bathymetry layer...")
          ip_extent <- terra::ext(config$indo_pacific_bbox)
          bath_global <- tryCatch(terra::crop(bath_global, ip_extent), error=function(e){
            hlog("WARN", paste("Failed cropping global bathy:", e$message, "- proceeding without crop."))
            bath_global # Use uncropped if error
          })
        }
        
        # Ensure CRS matches predictor stack
        if (terra::crs(bath_global) != raster_crs_terra) {
          hlog("DEBUG", "Projecting global bathymetry...")
          bath_global <- tryCatch(terra::project(bath_global, raster_crs_terra), error = function(e) {
            hlog("ERROR", paste("Failed project bathymetry:", e$message))
            NULL
          })
          if (is.null(bath_global)) {
            limit_by_depth_flag <- FALSE
            hlog("WARN", "Disabling depth filter because projection failed.")
          }
        }
        
        # --- Create mask based on fixed depth within OBIS extent ---
        if (limit_by_depth_flag) {
          # Crop bathymetry to the species' OBIS extent first
          bath_species_extent <- tryCatch(terra::crop(bath_global, obis_calibration_vect, snap="near"), error = function(e) {
            hlog("WARN", paste("Failed cropping bathymetry to OBIS extent for depth filter:", e$message))
            NULL
          })
          
          if (is.null(bath_species_extent)) {
            hlog("WARN", "Proceeding without depth filter due to cropping error.")
            limit_by_depth_flag <- FALSE
          } else {
            # Create the mask using fixed limits
            depth_mask <- bath_species_extent
            depth_mask[depth_mask < min_depth_fixed | depth_mask > max_depth_fixed] <- NA
            depth_mask[!is.na(depth_mask)] <- 1 # Set valid cells to 1
            
            # Check if the mask removed all valid cells
            max_val_check <- tryCatch(terra::global(depth_mask, "max", na.rm = TRUE)$max, error = function(e) NA) # Handle potential errors in global
            if (is.na(max_val_check) || max_val_check == 0) {
              hlog("WARN", "FIXED depth filter masked all cells within the OBIS ecoregion extent. Skipping depth filter.")
              depth_mask <- NULL
              limit_by_depth_flag <- FALSE # Disable filter if it masks everything
            } else {
              hlog("INFO", "FIXED depth filter mask created.")
            }
            rm(bath_species_extent); gc() # Clean up intermediate raster
          }
        }
        if(exists("bath_global")) rm(bath_global); gc() # Clean up global bathy
      }
    } # End else block (config values were valid)
  } else {
    hlog("INFO", "Depth filter disabled by config$limit_by_depth_obis.")
  }
  # --- End Depth Filtering Block ---
  
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


generate_sdm_background_obis_original <- function(occs_sf, global_predictor_stack, config, logger, species_log_file = NULL, seed_offset = 4) {
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
  
  
  # --- 6. Return Results ---
  hlog("INFO", "Returning background points, species stack, and polygon.\n")
  return(list(
    background_points = bg_coords_df,
    species_specific_stack = predictor_stack_obis_part2,
    calibration_polygon = obis_calibration_vect, # Polygon *before* depth masking
    quad_n_calculated = quad_n_final # Return the actual number sampled/calculated
  ))
  
  
  # return(bg_coords_df)
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
  
  
  hlog("INFO", "OKAY WE ARE DOING BLOCKCV STUFF")
  
  
  # --- Extract Config Parameters ---
  cv_method      <- config$sdm_spatial_cv_type_to_use
  
  nfolds         <- config$sdm_n_folds
  auto_range     <- config$blockcv_auto_range
  range          <- config$blockcv_range_default
  range_max      <- config$blockcv_range_max
  use_hexagon    <- config$blockcv_hexagon
  selection_type <- config$blockcv_selection
  n_iterate      <- config$blockcv_n_iterate
  lat_blocks     <- config$blockcv_lat_blocks
  
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
    hlog("INFO", paste("Auto range result: ", auto_cor$range))
    
    if (auto_cor$range < range_max) {
      range <- auto_cor$range
    } else {
      hlog("INFO", "cv_spatial_autocor result exceeded max range. Setting to max allowed range")
      range <- range_max
    }
  }
  
  
  if (cv_method == "spatial_grid") {
    hlog("INFO", "Generating grid blocks.")
    hlog("INFO", paste("Final range to use: ", range))
    spatial_folds <- blockCV::cv_spatial(r = predictor_stack, x = swd_sf, column = "pa", iteration = n_iterate, size = range,
                                         hexagon = use_hexagon, k = nfolds, progress = T, report = T, plot = T, selection = selection_type, extend = 0.5)
    blockCV::cv_plot(cv = spatial_folds, x = swd_sf, r = predictor_stack) + geom_sf(data = swd_sf, 
                                                                                    alpha = .5)
  }
  
  if (cv_method == "spatial_lat") {
    hlog("INFO", "Generating latitudinal blocks.")
    spatial_folds <- blockCV::cv_spatial(r = predictor_stack, x = swd_sf, column = "pa", iteration = n_iterate, rows_cols = c(lat_blocks, 1),
                                         hexagon = use_hexagon, k = nfolds, progress = T, report = T, plot = T, extend = 0.5)
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







# Corrected function for scripts/helpers/sdm_modeling_helpers.R

#' Log Final Model Evaluation Metrics (v4 - Corrected for Available CV Metrics)
#'
#' Calculates training AUC/TSS, AICc for the final model, and extracts
#' the mean test AUC from the cross-validation tuning results table.
#' Appends all metrics to a CSV log file.
#'
#' @param final_model The trained `SDMmodel` object.
#' @param full_swd_data The `SWD` object used to train the final model.
#' @param tuning_predictor_stack The `SpatRaster` stack used for training/AICc.
#' @param tuning_output The `SDMtune` object returned by `gridSearch`.
#' @param species_name_sanitized Sanitized species name.
#' @param group_name Group name ("anemone" or "anemonefish").
#' @param predictor_type_suffix Suffix indicating model type.
#' @param config The configuration list.
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
log_final_model_metrics <- function(final_model, full_swd_data, tuning_predictor_stack, tuning_output,
                                    species_name_sanitized, group_name, predictor_type_suffix,
                                    config, logger = NULL, species_log_file = NULL) {
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[LogMetricsHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) log4r::log(logger, level, msg) else {}}
  
  # --- Input validation ---
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Invalid SDMmodel provided."); return(FALSE) }
  if (is.null(full_swd_data) || !inherits(full_swd_data, "SWD")) { hlog("WARN", "SWD object missing, cannot calculate training metrics.") } # Warning instead of error
  if (is.null(tuning_predictor_stack) || !inherits(tuning_predictor_stack, "SpatRaster")) { hlog("ERROR", "Predictor stack required for AICc."); return(FALSE) }
  if (is.null(tuning_output) || !inherits(tuning_output, "SDMtune")) { hlog("ERROR", "Valid tuning_output object required for test metrics."); return(FALSE) }
  
  # --- Define Output File ---
  output_log_file <- file.path(config$log_dir_base, paste0("sdm_final_model_metrics_", group_name, ".csv"))
  hlog("DEBUG", paste("Logging final model metrics to:", output_log_file))
  
  # --- Initialize Metrics List ---
  # Only include metrics we can actually get
  metrics <- list(
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    GroupName = group_name,
    SpeciesName = species_name_sanitized,
    PredictorSuffix = predictor_type_suffix,
    AUC_train = NA_real_,
    TSS_train = NA_real_,
    AUC_test_CV = NA_real_, # Mean test AUC from CV tuning results
    AICc = NA_real_
    # Removed TSS_test_CV as it's not in tuning_output@results by default
  )
  
  # --- Calculate Training Metrics (from final_model using full_swd_data) ---
  if (!is.null(full_swd_data)) {
    tryCatch({ metrics$AUC_train <- SDMtune::auc(final_model, test = NULL) }, error = function(e) { hlog("WARN", paste("Could not calculate Training AUC:", e$message)) })
    tryCatch({ metrics$TSS_train <- SDMtune::tss(final_model, test = NULL) }, error = function(e) { hlog("WARN", paste("Could not calculate Training TSS:", e$message)) })
  } else {
    hlog("WARN", "Skipping training metric calculation due to missing SWD data.")
  }
  
  # --- Calculate AICc (from final_model using tuning_predictor_stack) ---
  tryCatch({
    metrics$AICc <- SDMtune::aicc(final_model, env = tuning_predictor_stack)
    hlog("INFO", metrics$AICc)
  }, error = function(e) { hlog("WARN", paste("Could not calculate AICc:", e$message)) })
  hlog("INFO", metrics$AICc)
  
  
  # --- Extract Test AUC Metric (from tuning_output results for the best model) ---
  tryCatch({
    res_df <- tuning_output@results
    if(!is.null(res_df) && nrow(res_df) > 0) {
      # Find the best row based on the metric used for tuning
      tuning_metric <- config$sdm_evaluation_metric %||% "AUC"
      metric_base_upper <- toupper(tuning_metric)
      target_metric_col <- if (tolower(tuning_metric) == "aicc") "AICc" else paste0("test_", metric_base_upper) # Should be test_AUC here
      
      best_row_index <- NULL
      if (target_metric_col %in% colnames(res_df)) {
        valid_metric_indices <- which(!is.na(res_df[[target_metric_col]]))
        if(length(valid_metric_indices) > 0) {
          select_fun <- if (target_metric_col == "AICc") which.min else which.max
          best_row_relative_index <- select_fun(res_df[[target_metric_col]][valid_metric_indices])
          best_row_index <- valid_metric_indices[best_row_relative_index]
        }
      }
      
      if(is.null(best_row_index)) {
        hlog("WARN", "Could not reliably determine best row index from tuning results to extract test AUC. Attempting fallback to first row.")
        fallback_idx <- which(!is.na(res_df$test_AUC))[1] # Fallback relies on test_AUC existing
        if(!is.na(fallback_idx)) {
          best_row_index <- fallback_idx
          hlog("WARN", "Using first valid row (index ", best_row_index, ") for test AUC.")
        } else {
          hlog("ERROR", "No valid rows/test_AUC found in tuning results table.")
          best_row_index <- NA
        }
      }
      
      # Extract test AUC from the best row if index is valid
      if(!is.na(best_row_index)){
        if ("test_AUC" %in% colnames(res_df)) {
          metrics$AUC_test_CV <- res_df$test_AUC[best_row_index]
          hlog("DEBUG", paste("  Extracted test_AUC (CV):", round(metrics$AUC_test_CV, 4)))
        } else {
          hlog("WARN", "test_AUC column not found in tuning results table.") # Should not happen if tuning metric is AUC
          metrics$AUC_test_CV <- NA # Ensure it's NA if column missing
        }
        # Cannot extract test_TSS as it's not present
        # metrics$TSS_test_CV <- NA
      } else {
        metrics$AUC_test_CV <- NA # Ensure NA if no valid row found
      }
      
    } else {
      hlog("WARN", "Tuning results table is empty or missing. Cannot extract test AUC.")
    }
  }, error = function(e) {
    hlog("ERROR", paste("Failed to extract test AUC from tuning_output:", e$message))
  })
  
  # --- Append to CSV File ---
  metrics_df <- as.data.frame(metrics)
  # Round numeric columns
  numeric_cols <- sapply(metrics_df, is.numeric)
  metrics_df[numeric_cols] <- lapply(metrics_df[numeric_cols], round, digits = 4)
  if("AICc" %in% names(metrics_df)) metrics_df$AICc <- round(metrics_df$AICc, 2)
  
  tryCatch({
    write_headers <- !file.exists(output_log_file)
    readr::write_csv(metrics_df, output_log_file, append = !write_headers, col_names = write_headers)
    hlog("INFO", paste("Successfully logged final model metrics for", species_name_sanitized))
    return(TRUE)
  }, error = function(e) {
    hlog("ERROR", paste("Failed to write metrics to log file:", output_log_file, "Error:", e$message))
    return(FALSE)
  })
}