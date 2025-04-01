# scripts/helpers/env_processing_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for Environmental Data Loading, Processing, and Analysis (VIF/PCA)
# Includes plotting functions from user's original scripts.
#-------------------------------------------------------------------------------
pacman::p_load(terra, sf, dplyr, readr, ggplot2, corrplot, tools, stringr) # Added stringr

#' Load and Stack Environmental Rasters for a Scenario
#'
#' Loads environmental rasters for a specific scenario, handling current,
#' future, and terrain data based on paths defined in the config. Stacks them
#' into a single SpatRaster. Handles future file structure where time steps
#' are in the filename within the SSP directory.
#'
#' @param scenario_name Character string matching a scenario name in `config$env_scenarios` (e.g., "current", "ssp119_2050").
#' @param config List containing the configuration settings.
#'
#' @return A SpatRaster object containing all layers for the scenario, or NULL on error.
#' @export
load_stack_env_data <- function(scenario_name, config) {
  cat("--- Loading environmental data for scenario:", scenario_name, "---\n")
  
  raster_files <- c()
  
  # --- Get Scenario Specific Files ---
  if (scenario_name == "current") {
    target_dir <- config$scenario_folder_map$current
    if (!dir.exists(target_dir)) {
      warning("Current scenario directory not found: ", target_dir, call. = FALSE)
      return(NULL)
    }
    # List all .tif files directly in the 'current' directory
    scenario_files <- list.files(target_dir, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE) # Changed recursive=FALSE
    scenario_files <- scenario_files[!grepl("\\.aux\\.xml$", scenario_files)]
    raster_files <- c(raster_files, scenario_files)
    cat("  Found", length(scenario_files), "files in", target_dir, "\n")
    
  } else if (grepl("ssp", scenario_name)) {
    # --- Handle Future Scenarios ---
    # 1. Determine the base SSP directory from the map
    target_dir <- config$scenario_folder_map[[scenario_name]] # e.g., data/env/future/ssp119
    
    # 2. Determine the required time step identifier from the name
    time_step_tag <- if (grepl("_2050$", scenario_name)) {
      "_dec50_" # Matches the middle part of filenames like ..._dec50_ltmin.tif
    } else if (grepl("_2100$", scenario_name)) {
      "_dec100_" # Matches the middle part of filenames like ..._dec100_mean.tif
    } else {
      warning("Cannot determine time step (dec50/dec100) from scenario name:", scenario_name, call. = FALSE)
      return(NULL)
    }
    cat("  Targeting time step tag:", time_step_tag, "\n")
    
    if (!dir.exists(target_dir)) {
      warning("Future scenario base directory not found: ", target_dir, call. = FALSE)
      return(NULL)
    }
    
    # 3. List ALL .tif files in the base SSP directory
    all_ssp_files <- list.files(target_dir, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE) # Changed recursive=FALSE
    all_ssp_files <- all_ssp_files[!grepl("\\.aux\\.xml$", all_ssp_files)]
    
    # 4. Filter files based on the time step tag in the filename
    scenario_files <- all_ssp_files[stringr::str_detect(basename(all_ssp_files), fixed(time_step_tag))]
    
    if(length(scenario_files) == 0) {
      warning("No files found containing the tag '", time_step_tag, "' in directory: ", target_dir, call. = FALSE)
      # Don't return NULL yet, still need to add CHL and Terrain
    } else {
      raster_files <- c(raster_files, scenario_files)
      cat("  Found", length(scenario_files), "files matching time step '", time_step_tag, "' in ", target_dir, "\n", sep="")
    }
    
    # 5. Add current Chlorophyll-a layer if it exists
    if (!is.null(config$chl_file_path) && file.exists(config$chl_file_path)) {
      cat("  Adding current Chlorophyll layer:", basename(config$chl_file_path), "\n")
      raster_files <- c(raster_files, config$chl_file_path)
    } else {
      cat("  Warning: Current Chlorophyll layer not found at:", config$chl_file_path, "\n")
    }
    
  } else {
    warning("Unknown scenario name:", scenario_name, call. = FALSE)
    return(NULL)
  }
  
  # --- Get Terrain Files (Same as before) ---
  if (!is.null(config$terrain_folder) && dir.exists(config$terrain_folder)) {
    terrain_files <- list.files(config$terrain_folder, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE) # Changed recursive=FALSE
    terrain_files <- terrain_files[!grepl("\\.aux\\.xml$", terrain_files)]
    terrain_to_include <- sapply(config$terrain_variables_final, function(var_name) {
      grep(paste0("^", var_name, "\\.tif$"), basename(terrain_files), value = TRUE, ignore.case=TRUE) # Match exact basename
    })
    terrain_to_include <- unlist(terrain_to_include)
    terrain_to_include <- terrain_to_include[!sapply(terrain_to_include, is.null) & !duplicated(terrain_to_include)]
    # Get full paths
    terrain_paths_to_include <- file.path(config$terrain_folder, terrain_to_include)
    
    
    if (length(terrain_paths_to_include) > 0) {
      cat("  Adding", length(terrain_paths_to_include), "terrain files:", paste(basename(terrain_paths_to_include), collapse=", "), "\n")
      raster_files <- c(raster_files, terrain_paths_to_include)
    } else {
      cat("  No terrain files matching config$terrain_variables_final found directly in", config$terrain_folder, "\n")
    }
  } else {
    warning("Terrain folder not found or not specified in config:", config$terrain_folder, call. = FALSE)
  }
  
  # --- Stack Rasters ---
  raster_files <- unique(raster_files)
  if (length(raster_files) == 0) {
    warning("No raster files identified to stack for scenario:", scenario_name, call. = FALSE)
    return(NULL)
  }
  
  cat("  Attempting to stack", length(raster_files), "raster files...\n")
  # print(basename(raster_files)) # Optional: print files being stacked
  tryCatch({
    env_stack <- terra::rast(raster_files)
    # Clean layer names: remove time step tags, baseline tags etc. for consistency
    clean_names <- names(env_stack)
    clean_names <- gsub("_dec50_|_dec100_", "_", clean_names) # Remove time step tag
    clean_names <- gsub("_baseline_\\d{4}_\\d{4}", "", clean_names) # Remove baseline year tag
    clean_names <- gsub("_\\d{4}_\\d{4}", "", clean_names) # Remove other year tags if present
    # Use file path sans ext on basename for final cleanup
    clean_names <- tools::file_path_sans_ext(basename(clean_names))
    names(env_stack) <- clean_names
    cat("  Final stacked raster names:", paste(names(env_stack), collapse=", "), "\n")
    # Check for duplicate names after cleaning
    if(any(duplicated(names(env_stack)))) {
      warning("Duplicate layer names found after cleaning for scenario ", scenario_name, ": ",
              paste(names(env_stack)[duplicated(names(env_stack))], collapse=", "),
              ". This may cause issues in analysis.", call.=FALSE)
    }
    return(env_stack)
  }, error = function(e) {
    warning("Error stacking rasters for scenario ", scenario_name, ": ", e$message, call. = FALSE)
    cat("Files attempted to stack:\n")
    print(basename(raster_files))
    return(NULL)
  })
}


# --- Other Helper Functions (preprocess_env_rasters, load_clean_thin_group_occurrences, extract_env_values, plotting functions) ---
# --- remain unchanged from the previous version ---

#' Preprocess Environmental Rasters
preprocess_env_rasters <- function(env_stack, config) {
  cat("--- Preprocessing environmental rasters ---\n")
  processed_stack <- env_stack
  if (config$apply_coral_mask) {
    cat("  Applying coral reef mask/crop...\n")
    if (!is.null(config$coral_shapefile) && file.exists(config$coral_shapefile)) {
      tryCatch({
        coral_areas <- terra::vect(config$coral_shapefile)
        if (terra::crs(coral_areas) != terra::crs(processed_stack)) {
          cat("    Projecting shapefile CRS...\n")
          coral_areas <- terra::project(coral_areas, terra::crs(processed_stack))
        }
        processed_stack <- terra::crop(processed_stack, coral_areas, snap="near")
        processed_stack <- terra::mask(processed_stack, coral_areas)
        cat("  Coral mask/crop applied.\n")
      }, error = function(e) {warning("Error applying coral reef mask/crop: ", e$message, "\n  Skipping mask.", call. = FALSE)})
    } else {warning("Coral shapefile not found or not specified in config. Skipping mask.", call. = FALSE)}
  } else { cat("  Coral masking disabled in config.\n")}
  if (config$apply_depth_filter) {
    cat("  Applying depth filter (", config$depth_min, "m to ", config$depth_max, "m)...\n")
    depth_layer_name <- config$terrain_variables_final[grepl("bathymetry", config$terrain_variables_final, ignore.case=TRUE)][1]
    if (length(depth_layer_name) == 1 && depth_layer_name %in% names(processed_stack)) {
      tryCatch({
        depth_layer <- processed_stack[[depth_layer_name]]
        depth_mask <- depth_layer < config$depth_min | depth_layer > config$depth_max
        processed_stack[depth_mask] <- NA
        min_val <- terra::global(processed_stack[[1]], "min", na.rm=TRUE)$min
        if(is.na(min_val)){warning("All cells masked out by depth filter.", call. = FALSE)} else {cat("  Depth filter applied.\n")}
      }, error = function(e) {warning("Error applying depth filter: ", e$message, "\n  Skipping filter.", call. = FALSE)})
    } else {warning("Depth layer (expected name like 'bathymetry_mean') not found. Cannot apply depth filter.", call. = FALSE)}
  } else { cat("  Depth filtering disabled in config.\n")}
  cat("--- Raster preprocessing finished ---\n")
  return(processed_stack)
}

#' Load, Clean, and Thin Group Occurrence Data
load_clean_thin_group_occurrences <- function(group_occurrence_dir, target_crs, config, env_stack = NULL) {
  cat("--- Loading and cleaning occurrences from:", group_occurrence_dir, "---\n")
  if (!dir.exists(group_occurrence_dir)) {warning("Occurrence directory not found:", group_occurrence_dir, call. = FALSE); return(NULL)}
  occ_files <- list.files(group_occurrence_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(occ_files) == 0) {warning("No occurrence CSV files found in:", group_occurrence_dir, call. = FALSE); return(NULL)}
  tryCatch({
    all_occs_list <- lapply(occ_files, readr::read_csv, col_types = readr::cols(.default = "c"))
    all_occs_raw <- dplyr::bind_rows(all_occs_list)
    cat("  Loaded", nrow(all_occs_raw), "raw records from", length(occ_files), "files.\n")
  }, error = function(e) {warning("Error loading or binding occurrence files: ", e$message, call. = FALSE); return(NULL)})
  if (!all(c("decimalLongitude", "decimalLatitude") %in% names(all_occs_raw))) {warning("Required columns 'decimalLongitude' and 'decimalLatitude' not found.", call.=FALSE); return(NULL)}
  all_occs_clean <- all_occs_raw %>% mutate(decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)), decimalLatitude = suppressWarnings(as.numeric(decimalLatitude))) %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude), decimalLongitude >= -180, decimalLongitude <= 180, decimalLatitude >= -90, decimalLatitude <= 90)
  n_after_clean <- nrow(all_occs_clean); cat("  Retained", n_after_clean, "records after removing NA/invalid coordinates.\n")
  if (n_after_clean == 0) return(NULL)
  tryCatch({occs_sf <- sf::st_as_sf(all_occs_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs)}, error = function(e){warning("Error converting occurrences to sf object: ", e$message, call.=FALSE); return(NULL)})
  target_crs_obj <- sf::st_crs(target_crs)
  if (sf::st_crs(occs_sf) != target_crs_obj) {
    cat("  Transforming occurrence CRS to match raster CRS...\n")
    tryCatch({occs_sf <- sf::st_transform(occs_sf, crs = target_crs_obj)}, error = function(e){warning("Error transforming CRS: ", e$message, call.=FALSE); return(NULL)})
  }
  cat("  Applying spatial thinning (method:", config$thinning_method, ")...\n")
  if (config$thinning_method == "cell") {
    if (is.null(env_stack)) {warning("env_stack is required for 'cell' thinning method. Skipping thinning.", call.=FALSE); occs_thinned <- occs_sf}
    else { tryCatch({
      occs_spatvector <- terra::vect(occs_sf); occs_cells <- terra::extract(env_stack[[1]], occs_spatvector, cells = TRUE)
      valid_cells_indices <- which(!is.na(occs_cells$cell)); unique_cell_indices_in_valid <- which(!duplicated(occs_cells$cell[valid_cells_indices]))
      original_indices_to_keep <- valid_cells_indices[unique_cell_indices_in_valid]; occs_thinned <- occs_sf[original_indices_to_keep, ]
    }, error = function(e){warning("Error during cell-based thinning: ", e$message, "\n  Returning unthinned data.", call.=FALSE); occs_thinned <- occs_sf}) }
  } else { cat("  Thinning method '", config$thinning_method, "' not implemented or thinning disabled. Returning unthinned data.\n", sep=""); occs_thinned <- occs_sf}
  n_after_thin <- nrow(occs_thinned); cat("  Retained", n_after_thin, "records after thinning.\n")
  if (n_after_thin == 0) return(NULL)
  return(occs_thinned)
}

#' Extract Environmental Values at Occurrence Locations
extract_env_values <- function(occurrences_sf, env_stack) {
  cat("--- Extracting environmental values at", nrow(occurrences_sf), "locations ---\n")
  if (is.null(occurrences_sf) || is.null(env_stack) || nrow(occurrences_sf) == 0) {warning("Invalid input for extracting env values.", call. = FALSE); return(NULL)}
  tryCatch({
    if (sf::st_crs(occurrences_sf) != terra::crs(env_stack)) {warning("CRS mismatch between occurrences and rasters during extraction. Transforming occurrences.", call. = FALSE); occurrences_sf <- sf::st_transform(occurrences_sf, terra::crs(env_stack))}
    extracted_values <- terra::extract(env_stack, terra::vect(occurrences_sf))
    if ("ID" %in% names(extracted_values)) {extracted_values <- extracted_values[, -which(names(extracted_values) == "ID"), drop = FALSE]}
    initial_rows <- nrow(extracted_values); extracted_values_clean <- na.omit(extracted_values); final_rows <- nrow(extracted_values_clean)
    cat("  Removed", initial_rows - final_rows, "rows with NA values.\n"); cat("  Returning data frame with", final_rows, "rows and", ncol(extracted_values_clean), "variables.\n")
    if (final_rows == 0) {warning("No data remaining after removing NAs from extracted values.", call.=FALSE); return(NULL)}
    return(as.data.frame(extracted_values_clean))
  }, error = function(e) {warning("Error extracting environmental values: ", e$message, call. = FALSE); return(NULL)})
}

#' Plot VIF results using ggplot (from user's 05 script)
plot_vif_results_original <- function(vif_result, save_path = NULL) {
  if(!is.numeric(vif_result) || is.null(names(vif_result))){warning("plot_vif_results_original expects a named numeric vector (output of car::vif).", call.=FALSE); return(NULL)}
  df <- data.frame(Variable = names(vif_result), VIF = vif_result); num_vars <- nrow(df)
  if (num_vars == 0) return(NULL); colors <- grDevices::colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars)
  df$VIF <- pmax(0, df$VIF)
  p <- ggplot2::ggplot(df, aes(x = reorder(Variable, VIF), y = VIF)) + ggplot2::geom_bar(stat = "identity", aes(fill = Variable), show.legend = FALSE) + ggplot2::scale_fill_manual(values = colors) + ggplot2::labs(title = "VIF Analysis Results (car::vif)", x = "Environment Variables", y = "VIF Value") + ggplot2::scale_y_continuous(limits = c(0, max(5, ceiling(max(df$VIF, na.rm = TRUE)))), breaks = scales::pretty_breaks(n = 5)) + ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(save_path)) {tryCatch({ggplot2::ggsave(save_path, plot = p, width = 8 + 0.1 * num_vars, height = 6, limitsize = FALSE); cat("  VIF plot saved to", save_path, "\n")}, error = function(e) { warning("Failed to save VIF plot: ", e$message, call.=FALSE)})}
  return(p)
}

#' Plot Pearson Correlation using corrplot (from user's 05 script)
plot_correlation_results_original <- function(env_extract, save_path = NULL) {
  if(!is.data.frame(env_extract) && !is.matrix(env_extract) || ncol(env_extract) < 2) {warning("plot_correlation_results_original requires a data frame or matrix with at least 2 columns.", call.=FALSE); return(NULL)}
  env.cor <- stats::cor(env_extract, method = "pearson", use = "pairwise.complete.obs")
  env.p <- NULL; if(requireNamespace("Hmisc", quietly = TRUE)){env.p <- Hmisc::rcorr(as.matrix(env_extract), type="pearson")$P} else {warning("Hmisc package not found, p-values on corrplot may be missing or incomplete.", call.=FALSE)}
  tryCatch({ plot_obj <- function() { corrplot::corrplot(corr = env.cor, method="color", type = "upper", order = "hclust", p.mat = env.p, sig.level = c(.01, .05), insig = "label_sig", pch.cex = 1.5, pch.col = "grey", addCoef.col = "black", number.cex = 0.7, tl.col = "black", tl.srt = 45, diag = FALSE, na.label = "NA", mar=c(0,0,1,0)) }
  if (!is.null(save_path)) { grDevices::png(filename = save_path, width = 8, height = 8, units = "in", res = 300); plot_obj(); grDevices::dev.off(); cat("  Correlation plot saved to", save_path, "\n")} else { plot_obj() }
  return(invisible(NULL))
  }, error = function(e){warning("Failed to create or save correlation plot: ", e$message, call.=FALSE); return(NULL)})
}
#-------------------------------------------------------------------------------