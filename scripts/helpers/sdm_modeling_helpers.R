# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling using SDMtune (Revised Workflow v2)
# - Saving logic moved to dedicated helper functions.
# - Path construction uses new config structure for target output.
# - Added Variable Importance helper.
# - Incorporated SAC thinning and OBIS-style background generation returning list.
# - Added spatial CV fold creation and metric logging.
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr, blockCV, ggplot2, ecospat, spThin, geosphere,
               biomod2, presenceabsence, randomForest, gbm, mda, gam) # Added biomod2 and dependencies

#' Load, Clean, and Get Coordinates for Individual Species Occurrences
#' (Returns list with data and count) - Modified to NOT do cell thinning internally.
#'
#' @param species_aphia_id AphiaID of the target species.
#' @param occurrence_dir Directory containing individual species CSV files.
#' @param config Project configuration list (needs `occurrence_crs`, `min_occurrences_sdm`). Cell thinning config is ignored here.
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
      # Use standard logging levels from log4r if logger is provided
      if(level == "DEBUG") log4r::debug(logger, msg)
      else if(level == "INFO") log4r::info(logger, msg)
      else if(level == "WARN") log4r::warn(logger, msg)
      else log4r::error(logger, msg)
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
    hlog("DEBUG", paste("  Retained", count_after_clean, "records after basic coordinate cleaning (NA removal)."))
    
    if (count_after_clean == 0) {
      hlog("WARN", "No valid coordinates after cleaning.")
      return(list(coords = NULL, count = 0))
    }
    
    # --- NO SPATIAL THINNING HERE ---
    # Thinning is now handled by a separate function `thin_occurrences_by_sac`
    hlog("DEBUG", "  Skipping internal spatial thinning. Raw cleaned coordinates returned.")
    occ_final_coords <- as.matrix(occ_clean[, c("decimalLongitude", "decimalLatitude")])
    colnames(occ_final_coords) <- c("longitude", "latitude") # Ensure consistent names
    
    return(list(coords = occ_final_coords, count = count_after_clean))
    
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
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SacThinHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)} }
  
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
    mantel_res <- NULL # Ensure it's NULL on error
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
      thinned_list <- NULL # Return NULL on error
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


#' Generate Background Points using OBIS/MPAEU Ecoregion/Depth Method (Returns List v4)
#'
#' Replicates the OBIS/MPAEU Part 2 spatial extent definition AND background point sampling.
#' Returns a list containing background points and the species-specific masked stack.
#' Uses config values for depth limits directly if limit_by_depth_obis is TRUE.
#'
#' @param occs_sf sf object with cleaned, thinned occurrence points (needs CRS).
#' @param global_predictor_stack SpatRaster stack covering the global study area (needs CRS).
#' @param config Project configuration list. Requires relevant paths, flags, and
#'        `background_points_n` (used as the base target number).
#' @param logger A log4r logger object (can be NULL).
#' @param species_log_file Optional path to species-specific log file.
#' @param seed_offset Integer offset for background sampling seed.
#' @return A list containing `$background_points` (data frame of coordinates),
#'         `$species_specific_stack` (SpatRaster, masked to final extent),
#'         `$calibration_polygon` (SpatVector, pre-depth mask extent),
#'         `$quad_n_calculated` (Numeric, the number of background points sampled),
#'         or NULL on error.
generate_sdm_background_obis <- function(occs_sf, global_predictor_stack, config, logger, species_log_file = NULL, seed_offset = 4) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[GenerateBG_OBIS]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)} }
  
  # --- Input Checks ---
  if (!inherits(occs_sf, "sf") || nrow(occs_sf) < 3) { hlog("ERROR", "Valid sf occurrence object with >= 3 points required.\n"); return(NULL) }
  if (is.na(sf::st_crs(occs_sf))) { hlog("ERROR", "Occurrence sf object needs a CRS.\n"); return(NULL) }
  if (!inherits(global_predictor_stack, "SpatRaster") || terra::nlyr(global_predictor_stack) < 1) { hlog("ERROR", "Valid global SpatRaster required.\n"); return(NULL) }
  if (terra::crs(global_predictor_stack) == "") { hlog("ERROR", "Global predictor stack needs a CRS.\n"); return(NULL) }
  if (!file.exists(config$ecoregion_shapefile)) { hlog("ERROR", "Ecoregion shapefile not found:", config$ecoregion_shapefile, "\n"); return(NULL) }
  if (config$limit_by_depth_obis && !file.exists(config$bathymetry_file)) { hlog("ERROR", "Bathymetry file required for depth filtering not found:", config$bathymetry_file, "\n"); return(NULL) }
  
  # Get AphiaID if present, otherwise use a placeholder for seed offset
  species_aphia_id_internal <- tryCatch(occs_sf$AphiaID[1], error = function(e) { NA })
  if(is.na(species_aphia_id_internal)) {
    hlog("WARN", "AphiaID missing from occs_sf, using random seed component.")
    species_aphia_id_internal <- sample.int(1e6, 1)
  }
  
  # --- 1. Ecoregion Extent Definition ---
  # ... (Keep the ecoregion loading, intersection, adjacency, buffering logic) ...
  hlog("DEBUG", "Loading ecoregions...")
  ecoregions <- tryCatch(terra::vect(config$ecoregion_shapefile), error = function(e) { hlog("ERROR", "Failed load ecoregions:", e$message, "\n"); NULL })
  if (is.null(ecoregions)) return(NULL)
  occs_vect <- tryCatch(terra::vect(occs_sf), error = function(e) { hlog("ERROR", "Failed convert occs sf to vect:", e$message, "\n"); NULL })
  if(is.null(occs_vect)) return(NULL)
  if (terra::crs(occs_vect) != terra::crs(ecoregions)) {
    hlog("DEBUG", "Projecting occurrences to match ecoregions CRS...")
    occs_vect <- tryCatch(terra::project(occs_vect, terra::crs(ecoregions)), error = function(e) { hlog("ERROR", "Failed project occurrences:", e$message, "\n"); NULL })
    if(is.null(occs_vect)) return(NULL)
  }
  hlog("DEBUG", "Identifying intersecting ecoregions...")
  intersect_idx <- tryCatch(terra::is.related(ecoregions, occs_vect, "intersects"), error=function(e){hlog("WARN", "is.related failed:",e$message,"\n");NULL})
  if (is.null(intersect_idx) || !any(intersect_idx)) { hlog("WARN", "No occurrences intersect with ecoregions. Cannot define extent this way.\n"); return(NULL) }
  ecoreg_occ <- ecoregions[intersect_idx, ]
  hlog("DEBUG", "Finding adjacent ecoregions using buffered points...")
  sf::sf_use_s2(FALSE)
  occs_buffered_sf <- tryCatch(sf::st_buffer(sf::st_as_sf(occs_vect), dist = config$poly_buffer_obis), error = function(e) { hlog("ERROR", "Failed buffer points:", e$message, "\n"); NULL })
  if(is.null(occs_buffered_sf)) { sf::sf_use_s2(TRUE); return(NULL) }
  adjacent_idx <- tryCatch(terra::is.related(ecoregions, terra::vect(occs_buffered_sf), "intersects"), error=function(e){hlog("WARN", "is.related failed for adjacent:",e$message,"\n");NULL})
  if(is.null(adjacent_idx)) adjacent_idx <- intersect_idx # Fallback
  sf::sf_use_s2(TRUE)
  final_indices <- unique(c(which(intersect_idx), which(adjacent_idx)))
  ecoreg_sel <- ecoregions[final_indices, ]
  hlog("DEBUG", "Buffering final selected ecoregion polygon...")
  sf::sf_use_s2(FALSE)
  ecoreg_sel_buffered_sf <- tryCatch(sf::st_buffer(sf::st_as_sf(ecoreg_sel), dist = config$poly_buffer_final), error = function(e) { hlog("ERROR", "Failed buffer final polygon:", e$message, "\n"); NULL })
  if(is.null(ecoreg_sel_buffered_sf)) { sf::sf_use_s2(TRUE); return(NULL) }
  sf::sf_use_s2(TRUE)
  raster_crs_terra <- terra::crs(global_predictor_stack)
  obis_calibration_poly_sf <- sf::st_transform(ecoreg_sel_buffered_sf, crs = sf::st_crs(raster_crs_terra))
  obis_calibration_vect <- terra::vect(obis_calibration_poly_sf)
  hlog("INFO", "Ecoregion-based calibration polygon created (before depth filter).\n")
  
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
  
  quad_samp_target <- config$background_points_n # Base target number
  n_presences <- nrow(occs_sf) # Number of presence points (after cleaning/thinning)
  
  env_size_t <- tryCatch({valid_cells <- terra::global(predictor_stack_obis_part2[[1]], fun="notNA"); valid_cells$notNA}, error=function(e){ hlog("WARN", "Could not count valid cells for quad_n calc.\n"); NA})
  
  if(is.na(env_size_t) || env_size_t == 0){ hlog("ERROR", "No valid cells in final masked stack to sample background from.\n"); return(NULL) }
  hlog("DEBUG", "  Number of valid cells in final masked stack:", env_size_t, "\n")
  
  quad_n_final <- quad_samp_target
  if (n_presences >= quad_samp_target) {
    est_tsize <- n_presences * 2 # Potential doubling
    if (est_tsize > env_size_t) {
      quad_n_final <- env_size_t - n_presences
      if (quad_n_final < 10) quad_n_final <- min(10, env_size_t) # Ensure at least some points, capped by available cells
      hlog("WARN", "  Doubled presences exceed available cells. Capping background points to:", quad_n_final, "(Total cells:", env_size_t, ", Presences:", n_presences, ")\n")
    } else {
      quad_n_final <- n_presences * 2
      hlog("DEBUG", "  Using doubled presences as background target:", quad_n_final, "\n")
    }
  }
  quad_n_final <- min(quad_n_final, env_size_t)
  if(quad_n_final <= 0) { quad_n_final <- min(10, env_size_t); hlog("WARN", "Calculated quad_n <= 0, setting to 10 or available cells.")}
  hlog("INFO", "  Final calculated number of background points (quad_n):", quad_n_final, "\n")
  
  # --- 5. Sample Background Points ---
  hlog("INFO", "Generating background points using OBIS Part2 logic & stack...\n")
  # Ensure seed includes species ID and offset for reproducibility
  seed_value <- config$background_seed %||% 111
  set.seed(seed_value + species_aphia_id_internal + seed_offset)
  
  bg_obis_part2_df <- tryCatch({
    terra::spatSample(predictor_stack_obis_part2, size = quad_n_final, method = "random", na.rm = TRUE, xy = TRUE, warn = FALSE)
  }, error = function(e) { hlog("ERROR", "spatSample failed for background points:", e$message, "\n"); NULL })
  
  if (is.null(bg_obis_part2_df)) { hlog("ERROR", "  Background point generation failed.\n"); return(NULL) }
  if (nrow(bg_obis_part2_df) < quad_n_final) { hlog("WARN", "  Could only sample", nrow(bg_obis_part2_df), "background points (requested", quad_n_final, ").\n") }
  if (nrow(bg_obis_part2_df) == 0) { hlog("ERROR", "  Failed to generate ANY background points.\n"); return(NULL) }
  
  bg_coords_df <- as.data.frame(bg_obis_part2_df[, c("x", "y")])
  colnames(bg_coords_df) <- c("longitude", "latitude") # Standardize names
  hlog("INFO", "  Generated", nrow(bg_coords_df), "OBIS Part2 background points.\n")
  
  # --- 6. Return Results ---
  hlog("INFO", "Returning background points, species stack, and polygon.\n")
  return(list(
    background_points = bg_coords_df,
    species_specific_stack = predictor_stack_obis_part2,
    calibration_polygon = obis_calibration_vect, # Polygon *before* depth masking
    quad_n_calculated = nrow(bg_coords_df) # Return the actual number sampled
  ))
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
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[CreateCVFoldsSimp]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
  # --- Extract Config Parameters ---
  cv_method      <- config$sdm_spatial_cv_type_to_use
  nfolds         <- config$sdm_n_folds
  auto_range     <- config$blockcv_auto_range
  range_default  <- config$blockcv_range_default # Default fixed range
  range_max      <- config$blockcv_range_max     # Max allowed auto-range
  use_hexagon    <- config$blockcv_hexagon
  selection_type <- config$blockcv_selection
  n_iterate      <- config$blockcv_n_iterate
  # lat_blocks     <- config$blockcv_lat_blocks # Only needed if cv_method=="spatial_lat"
  
  hlog("INFO", paste("Attempting to create", nfolds, "spatial folds (Method:", cv_method, ", blockCV simplified logic)..."))
  
  # --- Prepare sf object from SWD data ---
  swd_coords_df <- as.data.frame(full_swd_data@coords); colnames(swd_coords_df) <- c("X", "Y")
  swd_coords_df$pa <- full_swd_data@pa
  target_crs_str <- terra::crs(predictor_stack, proj=TRUE); if (is.null(target_crs_str)|target_crs_str=="") target_crs_str <- "EPSG:4326"
  swd_sf <- tryCatch({ sf::st_as_sf(swd_coords_df, coords = c("X", "Y"), crs = target_crs_str) },
                     error = function(e){ hlog("ERROR", paste("Failed create sf from SWD:", e$message)); NULL })
  if(is.null(swd_sf)) return(NULL)
  hlog("DEBUG", paste("Created sf object with", nrow(swd_sf), "rows."))
  
  # --- Determine Range/Block Size ---
  final_range_m <- range_default # Start with default fixed range
  
  if (auto_range) {
    hlog("INFO", "Calculating spatial autocorrelation range...")
    sf::sf_use_s2(FALSE) # Disable S2 for blockCV functions if issues arise
    auto_cor_result <- tryCatch({
      blockCV::cv_spatial_autocor(x = swd_sf, column = "pa", plot = FALSE) # Avoid plotting here
    }, error = function(e) {
      hlog("WARN", paste("cv_spatial_autocor failed:", e$message, "- using default range."))
      NULL
    })
    sf::sf_use_s2(TRUE) # Re-enable S2 if disabled
    
    if (!is.null(auto_cor_result) && !is.na(auto_cor_result$range)) {
      auto_range_val <- auto_cor_result$range
      hlog("INFO", paste("  Autocorrelation range result:", round(auto_range_val, 1), "m"))
      if (auto_range_val < range_max) {
        final_range_m <- auto_range_val
        hlog("INFO", paste("  Using auto-calculated range:", round(final_range_m, 1), "m"))
      } else {
        final_range_m <- range_max
        hlog("WARN", paste("  Autocorrelation range exceeded max allowed (", range_max, "m). Using max range instead."))
      }
    } else {
      hlog("WARN", "  Autocorrelation calculation failed or returned NA. Using default range:", range_default, "m")
      final_range_m <- range_default
    }
  } else {
    hlog("INFO", paste("Using fixed default range:", range_default, "m"))
    final_range_m <- range_default
  }
  
  # Ensure range is positive
  if(final_range_m <= 0) {
    hlog("WARN", paste("Calculated/Default range is not positive (", final_range_m, "). Setting to a small positive value (1000m)."))
    final_range_m <- 1000
  }
  
  # --- Generate Folds based on Method ---
  spatial_folds <- NULL
  tryCatch({
    if (cv_method == "spatial_grid") {
      hlog("INFO", paste("Generating spatial grid blocks with size:", round(final_range_m, 1), "m"))
      spatial_folds <- blockCV::cv_spatial(
        x = swd_sf,
        column = "pa",
        r = predictor_stack[[1]], # Provide first layer for extent/res info if needed
        k = nfolds,
        size = final_range_m,
        hexagon = use_hexagon,
        selection = selection_type,
        iteration = n_iterate,
        progress = FALSE, # Disable internal progress bar
        report = FALSE,   # Disable internal report
        plot = FALSE      # Disable internal plot
      )
      # Optionally plot separately
      # blockCV::cv_plot(cv = spatial_folds, x = swd_sf, r = predictor_stack)
      
    } else if (cv_method == "spatial_lat") {
      hlog("INFO", "Generating latitudinal blocks.")
      lat_blocks_n <- config$blockcv_lat_blocks # Get from config
      if(is.null(lat_blocks_n) || !is.numeric(lat_blocks_n) || lat_blocks_n < 2) {
        hlog("WARN", "blockcv_lat_blocks invalid in config, defaulting to 5.")
        lat_blocks_n <- 5
      }
      spatial_folds <- blockCV::cv_spatial(
        x = swd_sf,
        column = "pa",
        r = predictor_stack[[1]],
        k = nfolds,
        rows_cols = c(lat_blocks_n, 1), # Use lat blocks from config
        selection = selection_type,
        iteration = n_iterate,
        progress = FALSE, report = FALSE, plot = FALSE
      )
    } else if (cv_method == "random") {
      hlog("WARN", "'random' CV method selected. This is standard k-fold, not spatial. Consider 'spatial_grid' or 'spatial_lat'.")
      # SDMtune handles random folds internally if folds=NULL or folds=integer 'k'
      # To explicitly return folds compatible with SDMtune's expectation for spatial CV:
      # We can create random folds using blockCV as well, but it's less standard.
      # For consistency, let SDMtune handle random k-fold by returning NULL here.
      # SDMtune will use its default random k-fold when `folds` is NULL in `train`.
      # However, the tuning function expects folds. Let's create standard random folds.
      set.seed(123) # for reproducibility
      n_pts <- nrow(swd_sf)
      fold_indices <- sample(rep(1:nfolds, length.out = n_pts))
      # Create the list structure SDMtune expects for folds
      spatial_folds <- list(
        train = lapply(1:nfolds, function(k) which(fold_indices != k)),
        test = lapply(1:nfolds, function(k) which(fold_indices == k))
      )
      attr(spatial_folds, "method") <- "Random k-fold" # Add attribute for clarity
      
    } else {
      hlog("ERROR", paste("Unsupported sdm_spatial_cv_type_to_use:", cv_method))
      return(NULL)
    }
    
    if(is.null(spatial_folds)) { hlog("ERROR", "blockCV fold generation returned NULL."); return(NULL) }
    
    hlog("DEBUG", paste("Spatial folds object created using blockCV (k=", nfolds, ")."))
    return(spatial_folds) # Return the blockCV object or random fold list
    
  }, error = function(e) {
    hlog("ERROR", paste("Error during spatial fold creation:", e$message))
    return(NULL)
  })
}


#' Tune SDM Hyperparameters using SDMtune gridSearch with SPATIAL k-fold CV
#' (Wrapper calling simplified blockCV helper)
run_sdm_tuning_scv <- function(occs_coords, predictor_stack, background_df, config, logger, species_name = "species", species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SCV_TuneHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
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
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[TrainHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
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
  
  base_filename <- paste0("mean_pred_", species_name_sanitized) # Using SDMtune mean convention implicitly
  target_dir <- NULL
  target_filename_stem <- NULL # Filename without extension
  
  if (scenario_name == "current") {
    target_dir <- config$target_predictions_current_dir
    target_filename_stem <- paste0(base_filename, predictor_type_suffix) # Add suffix
  } else {
    ssp_dir_name <- config$ssp_scenario_map[[scenario_name]]
    if (is.null(ssp_dir_name)) { warning("No SSP directory mapping for scenario: ", scenario_name, call. = FALSE); return(NULL) }
    target_dir <- file.path(config$target_predictions_future_base, ssp_dir_name)
    time_tag_clean <- ifelse(grepl("2050$", scenario_name), "dec50", ifelse(grepl("2100$", scenario_name), "dec100", "UNKNOWN"))
    if(time_tag_clean == "UNKNOWN") { warning("Cannot extract time tag from scenario: ", scenario_name, call. = FALSE); return(NULL) }
    target_filename_stem <- paste0(base_filename, "_", ssp_dir_name, "_", time_tag_clean, predictor_type_suffix) # Add suffix
  }
  
  if(is.null(target_dir) || is.null(target_filename_stem)) { warning("Could not construct target path components.", call. = FALSE); return(NULL) }
  target_filename <- paste0(target_filename_stem, ".tif")
  return(file.path(target_dir, target_filename))
}


#' Predict SDM Suitability (Returns Raster or Error Message)
#' @return A SpatRaster object with the prediction, or a character string with the error message on failure.
predict_sdm_suitability <- function(final_sdm_model, predictor_stack, config, logger, species_log_file = NULL, output_type = "cloglog") {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[PredictHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
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
#' @param tuning_output The object returned by `run_sdm_tuning_scv` (contains `@results`).
#' @param species_name_sanitized Sanitized species name (e.g., "Genus_species").
#' @param predictor_type_suffix Suffix indicating model type (e.g., "_pca", "_combined_pca").
#' @param config The configuration list (needs target_results_base, model_output_subdir_map).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_tuning_results <- function(tuning_output, species_name_sanitized, predictor_type_suffix, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveTuneHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
  if (is.null(tuning_output) || !inherits(tuning_output, "SDMtune")) { hlog("ERROR", "Invalid tuning object provided."); return(FALSE) }
  
  # Determine target subdirectory using the map from config
  subdir_name <- config$model_output_subdir_map[[predictor_type_suffix]]
  if (is.null(subdir_name)) { hlog("ERROR", paste("No output subdirectory mapping found for suffix:", predictor_type_suffix)); return(FALSE) }
  target_subdir <- file.path(config$target_results_base, subdir_name)
  dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
  
  # Construct filenames matching the target repo style
  target_base_name <- paste0("CV_Results_", species_name_sanitized) # Using sanitized name
  
  # 1. Save full tuning object (RDS - intermediate/diagnostic)
  # Let's save this to the intermediate results dir instead of target results
  intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(basename(config$anemone_occurrence_dir), predictor_type_suffix)) # Example based on group name
  dir.create(intermediate_results_dir, recursive = TRUE, showWarnings = FALSE)
  rds_file <- file.path(intermediate_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, "_object.rds"))
  tryCatch({ saveRDS(tuning_output, file = rds_file); hlog("DEBUG", "Full tuning object saved (RDS to intermediate):", basename(rds_file)) },
           error = function(e) { hlog("ERROR", paste("Failed save tuning RDS:", e$message)) }) # Log error but continue
  
  # 2. Save results table (CSV - for target analysis mirroring target)
  csv_file <- file.path(target_subdir, paste0(target_base_name, ".csv"))
  results_df <- tuning_output@results
  if (is.null(results_df) || nrow(results_df) == 0) { hlog("WARN", "No results table found in tuning object."); return(TRUE) } # Return TRUE as RDS might have saved
  
  tryCatch({ readr::write_csv(results_df, file = csv_file); hlog("DEBUG", "Tuning results table saved (CSV to target):", basename(csv_file)); TRUE },
           error = function(e) { hlog("ERROR", paste("Failed save tuning CSV:", e$message)); FALSE })
}


#' Save Final SDM Model Object (Intermediate Location)
#' Saves the trained SDMmodel object as an RDS file to the intermediate models directory.
#'
#' @param final_model The trained `SDMmodel` object.
#' @param species_name_sanitized Sanitized species name.
#' @param predictor_type_suffix Suffix indicating model type.
#' @param group_name Group name ("anemone" or "anemonefish") - used to determine subdir.
#' @param config The configuration list (needs `models_dir_intermediate`).
#' @param logger A log4r logger object.
#' @param species_log_file Optional path for detailed species logs.
#' @return TRUE on success, FALSE on failure.
save_final_model <- function(final_model, species_name_sanitized, predictor_type_suffix, group_name, config, logger, species_log_file = NULL) {
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SaveModelHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Invalid model object provided."); return(FALSE) }
  
  # Create model type subdirectory within intermediate models dir
  model_subdir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
  dir.create(model_subdir, recursive = TRUE, showWarnings = FALSE)
  
  model_file <- file.path(model_subdir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  
  tryCatch({ saveRDS(final_model, file = model_file); hlog("DEBUG", "Final model object saved (RDS to intermediate):", basename(model_file)); TRUE },
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
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[SavePredHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
  if (is.null(prediction_raster) || !inherits(prediction_raster, "SpatRaster")) { hlog("ERROR", "Invalid prediction raster provided."); return(FALSE) }
  
  # Construct the target filename using the dedicated helper
  pred_file_path <- construct_prediction_filename(species_name_sanitized, scenario_name, predictor_type_suffix, config)
  if (is.null(pred_file_path)) { hlog("ERROR", "Failed to construct target prediction path."); return(FALSE) }
  
  # Ensure target directory exists
  target_dir <- dirname(pred_file_path)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    terra::writeRaster(prediction_raster, filename = pred_file_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    hlog("DEBUG", paste("Prediction raster saved:", basename(pred_file_path), "to", target_dir))
    TRUE
  }, error = function(e) {
    hlog("ERROR", paste("Failed save prediction raster:", e$message))
    FALSE
  })
}


#' Calculate and Save Variable Importance (Handles MaxNet Directly)
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
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[VarImpHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
  # --- Input validation ---
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) {
    hlog("ERROR", "Invalid SDMmodel object provided for VI.")
    return(FALSE)
  }
  
  hlog("INFO", "Calculating/Extracting variable importance...")
  vi_results_df <- NULL
  
  tryCatch({
    model_method <- final_model@method # Get the method used (e.g., "Maxnet")
    
    if (model_method == "Maxnet") {
      hlog("DEBUG", "Maxnet model detected. Extracting pre-calculated permutation importance.")
      raw_maxnet_model <- final_model@model
      if (!is.null(raw_maxnet_model) && !is.null(raw_maxnet_model$variable.importance)) {
        importance_vector <- raw_maxnet_model$variable.importance
        vi_results_df <- data.frame(Variable = names(importance_vector), Importance = as.numeric(importance_vector), stringsAsFactors = FALSE)
        vi_results_df <- vi_results_df[order(-vi_results_df$Importance), ] # Order descending
        hlog("DEBUG", "Successfully extracted MaxNet variable importance.")
      } else { hlog("WARN", "Could not find pre-calculated variable importance in the MaxNet model object.") }
    } else {
      hlog("DEBUG", paste("Model method is", model_method, ". Attempting SDMtune::varImp..."))
      if (is.null(training_swd) || !inherits(training_swd, "SWD")) { hlog("ERROR", "Training SWD object required for varImp calculation with non-MaxNet models."); return(FALSE) }
      show_progress <- FALSE
      if (!is.null(logger)) { tryCatch({current_log_level_num <- log4r::log_level(config$log_level %||% "INFO"); debug_level_num <- log4r::log_level("DEBUG"); show_progress <- current_log_level_num <= debug_level_num}, error=function(e){show_progress<-FALSE}) }
      vi_results <- SDMtune::varImp(model = final_model, permut = 10, progress = show_progress, test = training_swd)
      if (!is.null(vi_results) && nrow(vi_results) > 0) { vi_results_df <- vi_results; hlog("DEBUG", "SDMtune::varImp successful for non-MaxNet model.") }
      else { hlog("WARN", "SDMtune::varImp returned empty results for non-MaxNet model.") }
    }
    
    # --- Saving Logic ---
    if (is.null(vi_results_df)) { hlog("ERROR", "Variable importance could not be obtained."); return(FALSE) }
    
    # Determine target subdirectory using the map from config
    subdir_name <- config$model_output_subdir_map[[predictor_type_suffix]] # Uses same map as CV results
    if (is.null(subdir_name)) { hlog("ERROR", paste("No output subdirectory mapping found for suffix:", predictor_type_suffix)); return(FALSE) }
    target_subdir <- file.path(config$target_vi_base, subdir_name) # Use config$target_vi_base
    dir.create(target_subdir, recursive = TRUE, showWarnings = FALSE)
    
    vi_filename <- paste0("VI_", species_name_sanitized, ".csv") # Filename consistent with target repo
    vi_file_path <- file.path(target_subdir, vi_filename)
    
    readr::write_csv(vi_results_df, vi_file_path)
    hlog("INFO", paste("Variable importance saved to:", vi_file_path))
    return(TRUE)
    
  }, error = function(e) {
    hlog("ERROR", paste("Variable importance processing/saving failed:", e$message))
    return(FALSE)
  })
}


#' Log Final Model Evaluation Metrics
#'
#' Calculates training AUC/TSS, AICc for the final model, and extracts
#' the mean test AUC from the cross-validation tuning results table.
#' Appends all metrics to a CSV log file.
#'
#' @param final_model The trained `SDMmodel` object.
#' @param full_swd_data The `SWD` object used to train the final model (can be NULL if model loaded).
#' @param tuning_predictor_stack The `SpatRaster` stack used for training/AICc.
#' @param tuning_output The `SDMtune` object returned by `gridSearch` (can be NULL if model loaded).
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
  
  hlog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[LogMetricsHelper]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else if (!is.null(logger)) {if(level=="INFO") log4r::info(logger, msg) else if(level=="WARN") log4r::warn(logger, msg) else if(level=="ERROR") log4r::error(logger, msg) else log4r::debug(logger, msg)}}
  
  # --- Input validation ---
  if (is.null(final_model) || !inherits(final_model, "SDMmodel")) { hlog("ERROR", "Invalid SDMmodel provided."); return(FALSE) }
  # SWD data is needed for training metrics, warn if missing
  if (is.null(full_swd_data) || !inherits(full_swd_data, "SWD")) { hlog("WARN", "SWD object missing, cannot calculate training metrics.") }
  # Stack is needed for AICc
  if (is.null(tuning_predictor_stack) || !inherits(tuning_predictor_stack, "SpatRaster")) { hlog("ERROR", "Predictor stack required for AICc."); return(FALSE) }
  # Tuning output is needed for test metrics
  if (is.null(tuning_output) || !inherits(tuning_output, "SDMtune")) { hlog("WARN", "Tuning_output object missing or invalid, cannot extract test metrics.") }
  
  # --- Define Output File ---
  # Log metrics centrally per group
  output_log_file <- file.path(config$log_dir_base, paste0("sdm_final_model_metrics_", group_name, ".csv"))
  hlog("DEBUG", paste("Logging final model metrics to:", output_log_file))
  
  # --- Initialize Metrics List ---
  metrics <- list(
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    GroupName = group_name,
    SpeciesName = species_name_sanitized,
    PredictorSuffix = predictor_type_suffix,
    AUC_train = NA_real_,
    TSS_train = NA_real_,
    AUC_test_CV = NA_real_, # Mean test AUC from CV tuning results
    TSS_test_CV = NA_real_, # Mean test TSS from CV tuning results
    AICc = NA_real_
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
  }, error = function(e) { hlog("WARN", paste("Could not calculate AICc:", e$message)) })
  
  # --- Extract Test Metrics (from tuning_output results for the best model) ---
  if (!is.null(tuning_output) && inherits(tuning_output, "SDMtune")) {
    tryCatch({
      res_df <- tuning_output@results
      if(!is.null(res_df) && nrow(res_df) > 0) {
        # Find the best row based on the metric used for tuning
        tuning_metric <- config$sdm_evaluation_metric %||% "AUC" # Default to AUC if not set
        metric_base_upper <- toupper(tuning_metric)
        target_metric_col <- if (tolower(tuning_metric) == "aicc") "AICc" else paste0("test_", metric_base_upper)
        
        best_row_index <- NA # Initialize
        if (target_metric_col %in% colnames(res_df)) {
          valid_metric_indices <- which(!is.na(res_df[[target_metric_col]]))
          if(length(valid_metric_indices) > 0) {
            select_fun <- if (target_metric_col == "AICc") which.min else which.max
            best_row_relative_index <- select_fun(res_df[[target_metric_col]][valid_metric_indices])
            best_row_index <- valid_metric_indices[best_row_relative_index]
          } else {
            hlog("WARN", paste("No valid values found for tuning metric", target_metric_col, "in tuning results."))
          }
        } else {
          hlog("WARN", paste("Tuning metric column", target_metric_col, "not found in tuning results table."))
        }
        
        # If best row found, extract test AUC and TSS if available
        if(!is.na(best_row_index)) {
          if ("test_AUC" %in% colnames(res_df)) {
            metrics$AUC_test_CV <- res_df$test_AUC[best_row_index]
            hlog("DEBUG", paste("  Extracted test_AUC (CV):", round(metrics$AUC_test_CV, 4)))
          } else { hlog("WARN", "test_AUC column not found in tuning results table.") }
          
          if ("test_TSS" %in% colnames(res_df)) {
            metrics$TSS_test_CV <- res_df$test_TSS[best_row_index]
            hlog("DEBUG", paste("  Extracted test_TSS (CV):", round(metrics$TSS_test_CV, 4)))
          } else { hlog("WARN", "test_TSS column not found in tuning results table.") }
          
        } else {
          hlog("WARN", "Could not determine best row index from tuning results to extract test metrics.")
        }
      } else { hlog("WARN", "Tuning results table is empty or missing. Cannot extract test metrics.") }
    }, error = function(e) { hlog("ERROR", paste("Failed to extract test metrics from tuning_output:", e$message)) })
  } else { hlog("WARN", "Skipping test metric extraction due to missing/invalid tuning_output.") }
  
  # --- Append to CSV File ---
  metrics_df <- as.data.frame(metrics)
  numeric_cols <- sapply(metrics_df, is.numeric)
  metrics_df[numeric_cols] <- lapply(metrics_df[numeric_cols], round, digits = 4) # Round numerics
  if("AICc" %in% names(metrics_df)) metrics_df$AICc <- round(metrics_df$AICc, 2) # Different rounding for AICc
  
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









# ===========================================================================
# --- BIOMOD2 Helper Functions (New Additions) ---
# ===========================================================================
#' Format Data for BIOMOD2
#' (Keep this function as previously defined - it's good)



# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling (SDMtune & BIOMOD2)
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr, 
               blockCV, ggplot2, ecospat, spThin, geosphere,
               biomod2, maxnet, # For BIOMOD2 MAXNET
               # Common dependencies for other biomod2 models if used later:
               presenceabsence, randomForest, gbm, mda, gam, earth, xgboost 
)

# --- SDMtune Specific Helpers (Keep all your existing ones) ---
# load_clean_individual_occ_coords
# thin_occurrences_by_sac
# generate_sdm_background_obis
# create_spatial_cv_folds_simplified (This is for SDMtune)
# run_sdm_tuning_scv
# train_final_sdm
# construct_prediction_filename # IMPORTANT: We will reuse this for BIOMOD2 naming
# predict_sdm_suitability
# save_tuning_results
# save_final_model
# save_sdm_prediction # IMPORTANT: We will adapt a BIOMOD2 version based on this
# calculate_and_save_vi
# log_final_model_metrics
# --- End of SDMtune Specific Helpers ---


# ===========================================================================
# --- BIOMOD2 Helper Functions ---
# ===========================================================================
# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for SDM Modeling (SDMtune & BIOMOD2)
#-------------------------------------------------------------------------------
pacman::p_load(SDMtune, terra, sf, dplyr, tools, stringr, log4r, readr, 
               blockCV, ggplot2, ecospat, spThin, geosphere,
               biomod2, maxnet, # For BIOMOD2 MAXNET
               presenceabsence, randomForest, gbm, mda, gam, earth, xgboost 
)

# --- SDMtune Specific Helpers (Keep all your existing ones) ---
# load_clean_individual_occ_coords
# thin_occurrences_by_sac
# generate_sdm_background_obis
# create_spatial_cv_folds_simplified 
# run_sdm_tuning_scv
# train_final_sdm
# construct_prediction_filename 
# predict_sdm_suitability
# save_tuning_results
# save_final_model
# save_sdm_prediction 
# calculate_and_save_vi
# log_final_model_metrics
# --- End of SDMtune Specific Helpers ---


# ===========================================================================
# --- BIOMOD2 Helper Functions ---
# ===========================================================================

#' Format Data for BIOMOD2
format_data_for_biomod2 <- function(species_name_for_biomod, pres_coords, bg_coords, env_stack, species_log_file = NULL) {
  log_prefix_b2_fmt <- paste(Sys.time(), paste0("[", species_name_for_biomod, "]"))
  b2_hlog_fmt <- function(level, ...) {
    msg_content <- paste0(..., collapse = " ")
    full_msg <- paste(log_prefix_b2_fmt, paste0("[BIOMOD2_FormatData]"), level, "-", msg_content)
    if (!is.null(species_log_file)) cat(full_msg, "\n", file = species_log_file, append = TRUE) else cat(full_msg, "\n")
  }
  b2_hlog_fmt("INFO", "Formatting data for BIOMOD2...")
  if (is.null(pres_coords) || nrow(pres_coords) == 0) { b2_hlog_fmt("ERROR", "Presence coords missing."); return(NULL) }
  if (is.null(bg_coords) || nrow(bg_coords) == 0) { b2_hlog_fmt("ERROR", "Background coords missing."); return(NULL) }
  if (is.null(env_stack) || !inherits(env_stack, "SpatRaster")) { b2_hlog_fmt("ERROR", "Env stack missing/invalid."); return(NULL) }
  tryCatch({
    myResp <- c(rep(1, nrow(pres_coords)), rep(0, nrow(bg_coords)))
    myXY <- rbind(as.data.frame(pres_coords), as.data.frame(bg_coords)) 
    colnames(myXY) <- c('x', 'y') # Ensure standard names biomod2 might expect
    myBiomodData <- biomod2::BIOMOD_FormatingData(
      resp.var = myResp, expl.var = env_stack, resp.xy = myXY, 
      resp.name = species_name_for_biomod, PA.nb.rep = 0 # No pseudo-absences generated here
    )
    b2_hlog_fmt("INFO", "BIOMOD_FormatingData successful.")
    return(myBiomodData)
  }, error = function(e) {
    b2_hlog_fmt("ERROR", paste("Error in BIOMOD_FormatingData:", e$message))
    return(NULL)
  })
}

#' Create Spatial CV Folds Table for BIOMOD2 using blockCV
create_biomod2_block_cv_table <- function(biomod_formatted_data, predictor_stack, config, species_log_file = NULL) {
  log_prefix_b2_cv <- paste(Sys.time(), paste0("[", biomod_formatted_data@sp.name, "]"))
  b2_hlog_cv <- function(level, ...) {
    msg_content <- paste0(..., collapse = " ")
    full_msg <- paste(log_prefix_b2_cv, "[BIOMOD2_BlockCV]", level, "-", msg_content)
    if (!is.null(species_log_file)) cat(full_msg, "\n", file = species_log_file, append = TRUE) else cat(full_msg, "\n")
  }
  b2_hlog_cv("INFO", "Creating blockCV spatial folds for BIOMOD2...")
  if (!requireNamespace("blockCV", quietly = TRUE)) { b2_hlog_cv("ERROR", "Package 'blockCV' is required."); return(NULL) }
  if (is.null(biomod_formatted_data)) { b2_hlog_cv("ERROR", "biomod_formatted_data is required."); return(NULL) }
  if (is.null(predictor_stack)) { b2_hlog_cv("ERROR", "predictor_stack is required for blockCV spatial context."); return(NULL) }
  
  tryCatch({
    pa_data_df <- data.frame(x = biomod_formatted_data@coord[,1], y = biomod_formatted_data@coord[,2], occ = biomod_formatted_data@data.species)
    target_crs_val <- terra::crs(predictor_stack, proj=TRUE); if(is.null(target_crs_val) || target_crs_val == "") target_crs_val <- "EPSG:4326"
    pa_data_sf <- sf::st_as_sf(pa_data_df, coords = c("x", "y"), crs = target_crs_val)
    
    block_range_val <- config$blockcv_range_default 
    if (config$blockcv_auto_range) {
      sf::sf_use_s2(FALSE); # Disable s2 for blockCV's internal distance calcs if issues arise
      auto_cor_info <- tryCatch(blockCV::cv_spatial_autocor(x = pa_data_sf, column = "occ", plot = FALSE, progress = FALSE), 
                                error = function(e){ b2_hlog_cv("WARN", paste("cv_spatial_autocor failed:", e$message)); NULL})
      sf::sf_use_s2(TRUE); # Re-enable s2
      if (!is.null(auto_cor_info) && !is.na(auto_cor_info$range) && auto_cor_info$range > 0) {
        block_range_val <- min(auto_cor_info$range, config$blockcv_range_max %||% Inf)
        b2_hlog_cv("INFO", paste("  Using auto-calculated block range for blockCV:", round(block_range_val), "m"))
      } else {
        b2_hlog_cv("WARN", paste("  Auto-range for blockCV failed or returned invalid. Using default:", config$blockcv_range_default, "m"))
        block_range_val <- config$blockcv_range_default
      }
    } else {
      b2_hlog_cv("INFO", paste("  Using fixed block range for blockCV:", block_range_val, "m"))
    }
    if(block_range_val <= 0) block_range_val <- 1000 # Ensure positive range
    
    # Generate folds for BIOMOD2
    scv_results <- blockCV::cv_spatial(x = pa_data_sf, column = "occ", r = predictor_stack[[1]], 
                                       k = config$sdm_n_folds, size = block_range_val, 
                                       selection = config$blockcv_selection %||% "random", 
                                       iteration = config$blockcv_n_iterate %||% 50, 
                                       biomod2 = TRUE, # Critical for BIOMOD2 table format
                                       progress = FALSE, plot = FALSE)
    
    biomod_cv_table_from_blockCV <- NULL
    if (!is.null(scv_results)) {
      if (!is.null(scv_results$biomodTable)) biomod_cv_table_from_blockCV <- scv_results$biomodTable
      else if (!is.null(scv_results$biomod_table)) biomod_cv_table_from_blockCV <- scv_results$biomod_table # Older blockCV version
    }
    
    if (is.null(biomod_cv_table_from_blockCV)) {
      b2_hlog_cv("ERROR", "blockCV::cv_spatial did not return a biomodTable component.");
      if(!is.null(scv_results)) b2_hlog_cv("DEBUG", paste("Names in scv_results:", paste(names(scv_results), collapse=", ")))
      return(NULL)
    }
    
    # BIOMOD2 expects column names like _allData_RUN1, _allData_RUN2, etc.
    colnames(biomod_cv_table_from_blockCV) <- paste0("_allData_RUN", seq_len(ncol(biomod_cv_table_from_blockCV)))
    
    b2_hlog_cv("INFO", "blockCV spatial folds table for BIOMOD2 created.");
    return(as.data.frame(biomod_cv_table_from_blockCV)) # Return as data.frame
    
  }, error = function(e) {
    b2_hlog_cv("ERROR", paste("Error creating BIOMOD2 block CV table:", e$message)); return(NULL)
  })
}



#' Prepare User-Defined Options for BIOMOD2's MAXNET (R version)
#' Translates SDMtune Maxnet hyperparameters to BIOMOD2 MAXNET format.
#' @param sdmtune_maxnet_hypers Data frame of best hyperparameters from SDMtune for Maxnet.
#' @param species_log_file Optional path to species-specific log file.
#' @return A list formatted for `bm_ModelingOptions(user.val = ...)` for MAXNET.
prepare_biomod2_maxnet_user_val <- function(sdmtune_maxnet_hypers, species_log_file = NULL) {
  hlog_b2 <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[BIOMOD2_PrepOpt_MAXNET]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else cat(msg, "\n")}
  
  if (is.null(sdmtune_maxnet_hypers) || !inherits(sdmtune_maxnet_hypers, "data.frame") || nrow(sdmtune_maxnet_hypers) != 1) {
    hlog_b2("WARN", "SDMtune Maxnet hypers NULL or invalid. Cannot prepare MAXNET user options."); return(NULL)
  }
  hlog_b2("INFO", "Preparing user-defined options for BIOMOD2 MAXNET (R)...")
  
  maxnet_r_args <- list()
  if ("reg" %in% names(sdmtune_maxnet_hypers)) {
    maxnet_r_args$regmult <- as.numeric(sdmtune_maxnet_hypers$reg)
  }
  if ("fc" %in% names(sdmtune_maxnet_hypers)) {
    maxnet_r_args$classes <- as.character(sdmtune_maxnet_hypers$fc)
  }
  
  if (length(maxnet_r_args) > 0) {
    hlog_b2("DEBUG", paste("  MAXNET params:", paste(names(maxnet_r_args), unlist(maxnet_r_args), collapse=", ")))
    return(list('_allData_allRun' = maxnet_r_args))
  } else {
    hlog_b2("WARN", "No 'reg' or 'fc' found in SDMtune Maxnet hypers. MAXNET will use defaults."); return(NULL)
  }
}

#' Prepare User-Defined Options for BIOMOD2's RF (Random Forest)
#' Translates SDMtune RF hyperparameters to BIOMOD2 RF format.
#' @param sdmtune_rf_hypers Data frame of best hyperparameters from SDMtune for RF.
#' @param species_log_file Optional path to species-specific log file.
#' @return A list formatted for `bm_ModelingOptions(user.val = ...)` for RF.
prepare_biomod2_rf_user_val <- function(sdmtune_rf_hypers, species_log_file = NULL) {
  hlog_b2 <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[BIOMOD2_PrepOpt_RF]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else cat(msg, "\n")}
  
  if (is.null(sdmtune_rf_hypers) || !inherits(sdmtune_rf_hypers, "data.frame") || nrow(sdmtune_rf_hypers) != 1) {
    hlog_b2("WARN", "SDMtune RF hypers NULL or invalid. Cannot prepare RF user options."); return(NULL)
  }
  hlog_b2("INFO", "Preparing user-defined options for BIOMOD2 RF...")
  
  rf_args <- list()
  if ("mtry" %in% names(sdmtune_rf_hypers)) {
    rf_args$mtry <- as.integer(sdmtune_rf_hypers$mtry)
  }
  # Add other RF parameters if tuned, e.g., ntree
  # rf_args$ntree <- 500 # Default for biomod2 RF if not specified
  
  if (length(rf_args) > 0) {
    hlog_b2("DEBUG", paste("  RF params:", paste(names(rf_args), unlist(rf_args), collapse=", ")))
    return(list('_allData_allRun' = rf_args))
  } else {
    hlog_b2("WARN", "No 'mtry' found in SDMtune RF hypers. RF will use defaults."); return(NULL)
  }
}

#' Prepare User-Defined Options for BIOMOD2's ANN (nnet)
#' Translates SDMtune ANN hyperparameters to BIOMOD2 ANN format.
#' @param sdmtune_ann_hypers Data frame of best hyperparameters from SDMtune for ANN.
#' @param species_log_file Optional path to species-specific log file.
#' @return A list formatted for `bm_ModelingOptions(user.val = ...)` for ANN.
prepare_biomod2_ann_user_val <- function(sdmtune_ann_hypers, species_log_file = NULL) {
  hlog_b2 <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[BIOMOD2_PrepOpt_ANN]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else cat(msg, "\n")}
  
  if (is.null(sdmtune_ann_hypers) || !inherits(sdmtune_ann_hypers, "data.frame") || nrow(sdmtune_ann_hypers) != 1) {
    hlog_b2("WARN", "SDMtune ANN hypers NULL or invalid. Cannot prepare ANN user options."); return(NULL)
  }
  hlog_b2("INFO", "Preparing user-defined options for BIOMOD2 ANN (nnet)...")
  
  ann_args <- list()
  if ("size" %in% names(sdmtune_ann_hypers)) {
    ann_args$size <- as.integer(sdmtune_ann_hypers$size)
  }
  if ("decay" %in% names(sdmtune_ann_hypers)) {
    ann_args$decay <- as.numeric(sdmtune_ann_hypers$decay)
  }
  # ann_args$maxit <- 100 # Default for biomod2 ANN if not specified
  # ann_args$MaxNWts <- 10000 # Default in biomod2
  
  if (length(ann_args) > 0) {
    hlog_b2("DEBUG", paste("  ANN params:", paste(names(ann_args), unlist(ann_args), collapse=", ")))
    return(list('_allData_allRun' = ann_args))
  } else {
    hlog_b2("WARN", "No 'size' or 'decay' found in SDMtune ANN hypers. ANN will use defaults."); return(NULL)
  }
}

#' Prepare User-Defined Options for BIOMOD2's GBM (Boosted Regression Trees)
#' Translates SDMtune BRT hyperparameters to BIOMOD2 GBM format.
#' @param sdmtune_brt_hypers Data frame of best hyperparameters from SDMtune for BRT.
#' @param species_log_file Optional path to species-specific log file.
#' @return A list formatted for `bm_ModelingOptions(user.val = ...)` for GBM.
prepare_biomod2_gbm_user_val <- function(sdmtune_brt_hypers, species_log_file = NULL) {
  hlog_b2 <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), "[BIOMOD2_PrepOpt_GBM]", paste0(..., collapse = " ")); if (!is.null(species_log_file)) cat(msg, "\n", file = species_log_file, append = TRUE) else cat(msg, "\n")}
  
  if (is.null(sdmtune_brt_hypers) || !inherits(sdmtune_brt_hypers, "data.frame") || nrow(sdmtune_brt_hypers) != 1) {
    hlog_b2("WARN", "SDMtune BRT hypers NULL or invalid. Cannot prepare GBM user options."); return(NULL)
  }
  hlog_b2("INFO", "Preparing user-defined options for BIOMOD2 GBM...")
  
  gbm_args <- list()
  if ("interaction.depth" %in% names(sdmtune_brt_hypers)) {
    gbm_args$interaction.depth <- as.integer(sdmtune_brt_hypers$interaction.depth)
  }
  if ("n.trees" %in% names(sdmtune_brt_hypers)) {
    gbm_args$n.trees <- as.integer(sdmtune_brt_hypers$n.trees)
  }
  if ("shrinkage" %in% names(sdmtune_brt_hypers)) {
    gbm_args$shrinkage <- as.numeric(sdmtune_brt_hypers$shrinkage)
  }
  if ("bag.fraction" %in% names(sdmtune_brt_hypers)) {
    gbm_args$bag.fraction <- as.numeric(sdmtune_brt_hypers$bag.fraction)
  }
  # gbm_args$distribution <- "bernoulli" # Default for binary data in biomod2
  
  if (length(gbm_args) > 0) {
    hlog_b2("DEBUG", paste("  GBM params:", paste(names(gbm_args), unlist(gbm_args), collapse=", ")))
    return(list('_allData_allRun' = gbm_args))
  } else {
    hlog_b2("WARN", "No relevant hypers found in SDMtune BRT results. GBM will use defaults."); return(NULL)
  }
}




#' Run BIOMOD2 Modeling with Specified Algorithms and Options
run_biomod2_models_with_blockcv <- function(biomod_formatted_data, biomod_model_options, models_to_run, biomod_cv_table,
                                            species_name_for_files, predictor_type_suffix, group_name_for_paths,
                                            config, species_log_file = NULL) {
  log_prefix_b2_run <- paste(Sys.time(), paste0("[", species_name_for_files, "]"))
  b2_hlog_run <- function(level, ...) { 
    msg_content <- paste0(..., collapse = " ")
    full_msg <- paste(log_prefix_b2_run, paste0("[BIOMOD2_RunModels:", paste(models_to_run, collapse="+"), "]"), level, "-", msg_content)
    if (!is.null(species_log_file)) cat(full_msg, "\n", file = species_log_file, append = TRUE) else cat(full_msg, "\n")
  }
  
  b2_hlog_run("INFO", paste("Running BIOMOD2 modeling for algorithms:", paste(models_to_run, collapse=", ")))
  
  cv_strat_to_use <- if (is.null(biomod_cv_table)) 'random' else 'user.defined'
  num_reps_to_use <- if (is.null(biomod_cv_table)) (config$sdm_n_folds %||% 2) else ncol(biomod_cv_table)
  cv_perc_val_to_use <- if (is.null(biomod_cv_table)) 0.8 else NULL # Only for 'random'
  b2_hlog_run("INFO", paste("  Using CV Strategy:", cv_strat_to_use, "with", num_reps_to_use, "reps/folds."))
  
  species_biomod_run_dir <- file.path(config$sdm_output_dir_intermediate, "biomod2_outputs", 
                                      paste0(group_name_for_paths, predictor_type_suffix, "_biomod2"), 
                                      species_name_for_files) 
  
  #TODO: Experiment with this!!
  config$species_biomod_run_dir <- species_biomod_run_dir
  dir.create(species_biomod_run_dir, recursive = TRUE, showWarnings = FALSE)
  
  # current_wd_b2 <- getwd()
  # setwd(species_biomod_run_dir) 
  b2_hlog_run("DEBUG", paste("Temp WD for BIOMOD_Modeling:", getwd()))
  
  myBiomodModelOut <- NULL
  tryCatch({
    myBiomodModelOut <- biomod2::BIOMOD_Modeling(
      bm.format = biomod_formatted_data,
      modeling.id = paste0('CombinedRun_', format(Sys.time(), "%Y%m%d%H%M%S")),
      models = models_to_run, 
      bm.options = biomod_model_options,
      CV.strategy = cv_strat_to_use,
      CV.nb.rep = num_reps_to_use, 
      CV.perc = cv_perc_val_to_use,
      CV.user.table = biomod_cv_table, 
      var.import = 0, # Set to >0 if you want BIOMOD2's VI
      metric.eval = c('ROC', 'TSS'), 
      do.full.models = TRUE, 
      seed.val = config$AphiaID_for_seed %||% 42 
    )
  }, error = function(e) {
    b2_hlog_run("ERROR", paste("Error in BIOMOD_Modeling:", e$message))
    # If error is related to MAXENT.Phillips / maxent.jar, it will show here
    if (grepl("maxent.jar", e$message, ignore.case = TRUE)) {
      b2_hlog_run("ERROR", "This error might be related to MAXENT.Phillips. Ensure it's fully removed if only MAXNET is intended.")
    }
  }, finally = {
    # setwd(current_wd_b2) 
    # b2_hlog_run("DEBUG", paste("Restored WD to:", current_wd_b2))
  })
  
  if (!is.null(myBiomodModelOut)) b2_hlog_run("INFO", "BIOMOD_Modeling completed.")
  else b2_hlog_run("ERROR", "BIOMOD_Modeling failed.")
  return(myBiomodModelOut)
}

#' Project BIOMOD2 Models to Current Scenario
project_biomod2_models_current <- function(biomod_model_out, env_stack_current,
                                           species_name_for_saving, predictor_type_suffix, 
                                           selected_algo, 
                                           config, species_log_file = NULL) {
  log_prefix_b2_proj <- paste(Sys.time(), paste0("[", species_name_for_saving, "]")) 
  b2_hlog_proj <- function(level, ...) { 
    msg_content <- paste0(..., collapse = " ")
    full_msg <- paste(log_prefix_b2_proj, paste0("[BIOMOD2_Project:", selected_algo, "]"), level, "-", msg_content)
    if (!is.null(species_log_file)) cat(full_msg, "\n", file = species_log_file, append = TRUE) else cat(full_msg, "\n")
  }
  b2_hlog_proj("INFO", "Attempting to project model for current scenario...")
  if (is.null(biomod_model_out) || !inherits(biomod_model_out, "BIOMOD.models.out")) { b2_hlog_proj("ERROR", "Invalid biomod_model_out."); return(NULL) }
  if (is.null(env_stack_current) || !inherits(env_stack_current, "SpatRaster")) { b2_hlog_proj("ERROR", "Invalid env_stack_current."); return(NULL) }
  
  biomod_internal_resp_name <- biomod_model_out@sp.name
  target_full_model_name <- paste0(biomod_internal_resp_name, "_allData_allRun_", selected_algo)
  model_to_project_path <- NULL
  
  if (target_full_model_name %in% biomod_model_out@models.computed) {
    model_to_project_path <- target_full_model_name
    b2_hlog_proj("INFO", paste("Identified full model for projection:", model_to_project_path))
  } else {
    # Fallback if _allData_allRun_ALGONAME isn't present (e.g. if only CV runs were kept by some setting)
    available_models_for_algo <- grep(paste0("_", selected_algo, "$"), biomod_model_out@models.computed, value = TRUE)
    if (length(available_models_for_algo) > 0) {
      model_to_project_path <- available_models_for_algo[1] # Take the first one found
      b2_hlog_proj("WARN", paste0("'_allData_allRun_' model for ", selected_algo, " not found. Using first available: ", model_to_project_path))
    } else {
      b2_hlog_proj("ERROR", paste("No models found for algorithm:", selected_algo, ". Available:", paste(biomod_model_out@models.computed, collapse=", ")))
      return(NULL)
    }
  }
  
  projection_run_name <- paste0("CurrentProjection_", selected_algo, "_", format(Sys.time(), "%Y%m%d%H%M%S"))
  
  # Get group name from config (set by calling script)
  current_group_name_from_config <- config$group_name_for_biomod_output %||% "unknown_group"
  # Construct expected path where BIOMOD_Modeling was run
  expected_biomod_run_wd <- file.path(config$sdm_output_dir_intermediate, "biomod2_outputs", 
                                      paste0(current_group_name_from_config, predictor_type_suffix, "_biomod2"), 
                                      species_name_for_saving) 
  
  if (!dir.exists(expected_biomod_run_wd)) {
    b2_hlog_proj("ERROR", paste("Expected BIOMOD_Modeling WD does not exist:", expected_biomod_run_wd, 
                                "Cannot perform projection as BIOMOD_Projection needs to run in that directory.")); return(NULL)
  }
  
  # original_wd_proj <- getwd()
  # setwd(expected_biomod_run_wd)
  b2_hlog_proj("DEBUG", paste("Temp WD for BIOMOD_Projection:", getwd()))
  
  biomod_projection_obj <- NULL
  projected_raster <- NULL
  
  tryCatch({
    biomod_projection_obj <- biomod2::BIOMOD_Projection(
      bm.mod = biomod_model_out, 
      proj.name = projection_run_name, 
      new.env = env_stack_current, 
      models.chosen = model_to_project_path, 
      metric.binary = NULL, # Don't create binary maps at this stage
      compress = TRUE, 
      build.clamping.mask = FALSE, # Usually set to FALSE for current projections
      output.format = ".tif"
    )
    projected_raster_raw <- biomod2::get_predictions(biomod_projection_obj)
    if (!is.null(projected_raster_raw) && inherits(projected_raster_raw, "SpatRaster") && terra::nlyr(projected_raster_raw) > 0) {
      projected_raster <- projected_raster_raw / 1000 # Scale to 0-1
      names(projected_raster) <- paste0("suitability_", selected_algo, "_current")
      b2_hlog_proj("INFO", "Projection created and scaled to 0-1.")
    } else {
      b2_hlog_proj("ERROR", "get_predictions returned invalid raster or no layers.")
    }
  }, error = function(e) {
    b2_hlog_proj("ERROR", paste("Error in BIOMOD_Projection/get_predictions:", e$message))
  }, finally = {
    # setwd(original_wd_proj)
    b2_hlog_proj("DEBUG", paste("Restored WD to:", original_wd_proj))
  })
  
  return(projected_raster)
}

#' Save BIOMOD2 Projection Raster using SDMtune Naming Convention
save_biomod2_projection <- function(prediction_raster, species_name_sanitized_for_files, scenario_name, 
                                    model_type_suffix_for_saving, 
                                    config, logger, species_log_file = NULL) {
  log_prefix_b2_save <- paste(Sys.time(), paste0("[", species_name_sanitized_for_files, "]"))
  b2_hlog_save <- function(level, ...) {
    msg_content <- paste0(..., collapse = " ")
    full_msg <- paste(log_prefix_b2_save, paste0("[SaveBIOMOD2Pred]"), level, "-", msg_content)
    if (!is.null(species_log_file)) cat(full_msg, "\n", file = species_log_file, append = TRUE) else cat(full_msg, "\n")
  }
  if (is.null(prediction_raster) || !inherits(prediction_raster, "SpatRaster")) { b2_hlog_save("ERROR", "Invalid raster."); return(FALSE) }
  
  # Use the existing SDMtune helper for consistent naming and path structure
  pred_file_path <- construct_prediction_filename(species_name_sanitized_for_files, scenario_name, model_type_suffix_for_saving, config)
  if (is.null(pred_file_path)) { b2_hlog_save("ERROR", "Failed to construct BIOMOD2 target prediction path using construct_prediction_filename."); return(FALSE) }
  
  target_dir <- dirname(pred_file_path)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    terra::writeRaster(prediction_raster, filename = pred_file_path, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
    b2_hlog_save("DEBUG", paste("BIOMOD2 Projection raster saved:", basename(pred_file_path), "to", target_dir))
    TRUE
  }, error = function(e) { b2_hlog_save("ERROR", paste("Failed save BIOMOD2 prediction raster:", e$message)); FALSE })
}