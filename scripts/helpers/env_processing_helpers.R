# scripts/helpers/env_processing_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for Environmental Data Loading, Processing, and Analysis (VIF/PCA)
# Includes plotting functions from user's original scripts.
#-------------------------------------------------------------------------------
pacman::p_load(terra, sf, dplyr, readr, ggplot2, corrplot, Hmisc, tools, stringr) # Added stringr

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
  
  scenario_specific_files <- c() # Files specific to current/future time step
  terrain_paths_to_stack <- c() # Terrain files found
  
  # --- Get Scenario Specific Files (Current or Future) ---
  target_dir <- config$scenario_folder_map[[scenario_name]]
  if (is.null(target_dir) || !dir.exists(target_dir)) {
    warning("Scenario directory not found/mapped for '", scenario_name, "' at path: ", target_dir); return(NULL)
  }
  
  if (scenario_name == "current") {
    current_files <- list.files(target_dir, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE)
    current_files <- current_files[!grepl("\\.aux\\.xml$", current_files)]
    scenario_specific_files <- c(scenario_specific_files, current_files)
    cat("  Found", length(current_files), "files in current directory:", target_dir, "\n")
  } else if (grepl("ssp", scenario_name)) {
    time_step_tag <- if (grepl("_2050$", scenario_name)) "_dec50_" else if (grepl("_2100$", scenario_name)) "_dec100_" else { warning("Invalid future scenario name format."); return(NULL) }
    cat("  Targeting time step tag:", time_step_tag, "in directory:", target_dir, "\n")
    all_ssp_files <- list.files(target_dir, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE)
    all_ssp_files <- all_ssp_files[!grepl("\\.aux\\.xml$", all_ssp_files)]
    # Files matching time tag
    future_time_files <- all_ssp_files[stringr::str_detect(basename(all_ssp_files), fixed(time_step_tag))]
    if(length(future_time_files) > 0) { scenario_specific_files <- c(scenario_specific_files, future_time_files); cat("  Found", length(future_time_files), "files matching time step.\n") } else { warning("No files found containing tag '", time_step_tag, "'") }
    # Find copied/renamed current vars
    copied_vars_to_find <- config$vars_to_copy_to_future
    ssp_name <- sub("_.*", "", scenario_name); time_tag_clean <- gsub("^_|_$", "", time_step_tag)
    expected_copied_names <- sapply(copied_vars_to_find, function(var_basename){
      new_basename_step1 <- stringr::str_replace(var_basename, "_baseline_\\d{4}_\\d{4}", "")
      parts <- stringr::str_split(new_basename_step1, "_")[[1]]
      if (length(parts) >= 3) { paste(parts[1], ssp_name, parts[grepl("depth", parts)], time_tag_clean, parts[length(parts)], sep = "_") } else { paste(new_basename_step1, ssp_name, time_tag_clean, sep = "_") }
    })
    found_copied_files <- all_ssp_files[basename(tools::file_path_sans_ext(all_ssp_files)) %in% expected_copied_names]
    if(length(found_copied_files) > 0) { cat("  Found", length(found_copied_files), "copied current variables:", paste(basename(found_copied_files), collapse=", "), "\n"); scenario_specific_files <- c(scenario_specific_files, found_copied_files) } else { cat("  Note: Did not find expected copied/renamed current variables.\n") }
  } else { warning("Unknown scenario name:", scenario_name); return(NULL) }
  
  # --- Get Terrain Files (Simplified Logic) ---
  cat("  Checking for terrain files...\n")
  if (!is.null(config$terrain_folder) && dir.exists(config$terrain_folder)) {
    # Loop through the final names expected (e.g., "bathymetry_mean", "slope", "rugosity", "distcoast")
    for (terrain_var_final_name in config$terrain_variables_final) {
      expected_terrain_file <- file.path(config$terrain_folder, paste0(terrain_var_final_name, ".tif"))
      if (file.exists(expected_terrain_file)) {
        terrain_paths_to_stack <- c(terrain_paths_to_stack, expected_terrain_file)
        cat("    Found terrain file:", basename(expected_terrain_file), "\n")
      } else {
        cat("    Terrain file not found:", basename(expected_terrain_file), "in", config$terrain_folder, "\n")
      }
    }
  } else { warning("Terrain folder not found or not specified.", call. = FALSE) }
  if(length(terrain_paths_to_stack) == 0) {
    cat("  Warning: No terrain files were found to add.\n")
  }
  
  # --- Combine and Stack ---
  all_raster_files <- unique(c(scenario_specific_files, terrain_paths_to_stack))
  if (length(all_raster_files) == 0) { warning("No raster files identified for scenario:", scenario_name); return(NULL) }
  
  cat("  Attempting to stack", length(all_raster_files), "raster files...\n")
  # print(basename(all_raster_files)) # Debug
  tryCatch({
    env_stack <- terra::rast(all_raster_files)
    layer_names <- tools::file_path_sans_ext(basename(sources(env_stack)))
    # Quick check for obvious duplicates before assigning
    if(any(duplicated(layer_names))){
      warning("Duplicate source file basenames detected before stacking for scenario ", scenario_name, ". Check input file lists and paths.", call.=FALSE)
      # Potentially try unique() again, although it might hide path issues
      # env_stack <- terra::rast(unique(all_raster_files))
      # layer_names <- tools::file_path_sans_ext(basename(sources(env_stack)))
    }
    names(env_stack) <- layer_names
    cat("  Final stacked raster names:", paste(names(env_stack), collapse=", "), "\n")
    if(any(duplicated(names(env_stack)))) { warning("Duplicate layer names exist in the final stack for scenario ", scenario_name, ". This will likely cause errors downstream.") }
    return(env_stack)
  }, error = function(e) { warning("Error stacking rasters: ", e$message); print(basename(all_raster_files)); return(NULL) })
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

#' #' Plot VIF results using ggplot (from user's 05 script)
#' plot_vif_results_original <- function(vif_result, save_path = NULL) {
#'   if(!is.numeric(vif_result) || is.null(names(vif_result))){warning("plot_vif_results_original expects a named numeric vector (output of car::vif).", call.=FALSE); return(NULL)}
#'   df <- data.frame(Variable = names(vif_result), VIF = vif_result); num_vars <- nrow(df)
#'   if (num_vars == 0) return(NULL); colors <- grDevices::colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars)
#'   df$VIF <- pmax(0, df$VIF)
#'   p <- ggplot2::ggplot(df, aes(x = reorder(Variable, VIF), y = VIF)) + ggplot2::geom_bar(stat = "identity", aes(fill = Variable), show.legend = FALSE) + ggplot2::scale_fill_manual(values = colors) + ggplot2::labs(title = "VIF Analysis Results (car::vif)", x = "Environment Variables", y = "VIF Value") + ggplot2::scale_y_continuous(limits = c(0, max(5, ceiling(max(df$VIF, na.rm = TRUE)))), breaks = scales::pretty_breaks(n = 5)) + ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'   if (!is.null(save_path)) {tryCatch({ggplot2::ggsave(save_path, plot = p, width = 8 + 0.1 * num_vars, height = 6, limitsize = FALSE); cat("  VIF plot saved to", save_path, "\n")}, error = function(e) { warning("Failed to save VIF plot: ", e$message, call.=FALSE)})}
#'   return(p)
#' }
#' 
#' #' Plot Pearson Correlation using corrplot (from user's 05 script)
#' plot_correlation_results_original <- function(env_extract, save_path = NULL) {
#'   if(!is.data.frame(env_extract) && !is.matrix(env_extract) || ncol(env_extract) < 2) {warning("plot_correlation_results_original requires a data frame or matrix with at least 2 columns.", call.=FALSE); return(NULL)}
#'   env.cor <- stats::cor(env_extract, method = "pearson", use = "pairwise.complete.obs")
#'   env.p <- NULL; if(requireNamespace("Hmisc", quietly = TRUE)){env.p <- Hmisc::rcorr(as.matrix(env_extract), type="pearson")$P} else {warning("Hmisc package not found, p-values on corrplot may be missing or incomplete.", call.=FALSE)}
#'   tryCatch({ plot_obj <- function() { corrplot::corrplot(corr = env.cor, method="color", type = "upper", order = "hclust", p.mat = env.p, sig.level = c(.01, .05), insig = "label_sig", pch.cex = 1.5, pch.col = "grey", addCoef.col = "black", number.cex = 0.7, tl.col = "black", tl.srt = 45, diag = FALSE, na.label = "NA", mar=c(0,0,1,0)) }
#'   if (!is.null(save_path)) { grDevices::png(filename = save_path, width = 8, height = 8, units = "in", res = 300); plot_obj(); grDevices::dev.off(); cat("  Correlation plot saved to", save_path, "\n")} else { plot_obj() }
#'   return(invisible(NULL))
#'   }, error = function(e){warning("Failed to create or save correlation plot: ", e$message, call.=FALSE); return(NULL)})
#' }

#' Generate Scenario-Specific Variable Names from Core List (Corrected v2)
#'
#' Takes a list of core variable names (typically baseline/simplest form)
#' and a target scenario name, then generates the expected technical filenames/layer names.
#' Handles specific naming convention transformation.
#'
#' @param core_vars Character vector of core variable names (e.g., "thetao_baseline_depthmax_mean", "rugosity").
#' @param scenario_name Character string of the target scenario (e.g., "current", "ssp119_2050", "ssp585_2100").
#' @param config The project configuration list.
#'
#' @return Character vector of variable names expected for the given scenario.
#' @export
generate_scenario_variable_list <- function(core_vars, scenario_name, config) {
  
  scenario_vars <- sapply(core_vars, function(core_var) {
    
    # Handle terrain vars (no scenario tags)
    if (core_var %in% config$terrain_variables_final) {
      return(core_var)
    }
    
    # Handle current scenario
    if (scenario_name == "current") {
      # Expect core_var is already the baseline name
      if (!grepl("_baseline", core_var)) {
        warning("Core variable '", core_var, "' doesn't seem to have '_baseline'. Using as is for 'current'.", call. = FALSE)
      }
      return(core_var)
    }
    
    # Handle future scenarios
    if (grepl("ssp", scenario_name)) {
      ssp_match <- stringr::str_match(scenario_name, "(ssp\\d{3})_(\\d{4})")
      if(is.na(ssp_match[1,1])) {
        warning("Could not parse SSP/Year from scenario name: ", scenario_name, call. = FALSE)
        return(NA) # Indicate failure for this var/scenario combo
      }
      ssp_code <- ssp_match[1, 2] # e.g., "ssp119"
      year_code <- ssp_match[1, 3] # e.g., "2050"
      time_tag_clean <- ifelse(year_code == "2050", "dec50", "dec100")
      
      # Check if it's a variable that was copied/renamed (like chl, par)
      is_copied_var <- FALSE
      copied_stem <- NULL
      for (cpy_stem in config$vars_to_copy_to_future) {
        if (startsWith(core_var, cpy_stem)) {
          is_copied_var <- TRUE
          copied_stem <- cpy_stem # The original baseline stem
          break
        }
      }
      
      if (is_copied_var) {
        # Construct the expected *renamed* file stem
        # e.g., "chl_baseline_depthmax_mean" -> "chl_ssp119_depthmax_dec50_mean"
        new_basename_step1 <- stringr::str_replace(copied_stem, "_baseline_\\d{4}_\\d{4}", "") # Remove years if present
        new_basename_step1 <- stringr::str_replace(new_basename_step1, "_baseline", "") # Remove baseline tag
        parts <- stringr::str_split(new_basename_step1, "_")[[1]]
        if (length(parts) >= 3) { # e.g., chl, depthmax, mean
          var_part <- parts[1]
          depth_part <- parts[grepl("depth", parts)]
          stat_part <- parts[length(parts)]
          scenario_var_name <- paste(var_part, ssp_code, depth_part, time_tag_clean, stat_part, sep = "_")
        } else {
          warning("Unexpected structure for copied var '", copied_stem, "'. Simpler renaming.", call.=FALSE)
          scenario_var_name <- paste(new_basename_step1, ssp_code, time_tag_clean, sep = "_")
        }
      } else {
        # Standard transformation for variables with future data
        # e.g., "thetao_baseline_depthmax_mean" -> "thetao_ssp119_depthmax_dec50_mean"
        scenario_var_name <- gsub("_baseline_\\d{4}_\\d{4}", "", core_var) # Remove years
        scenario_var_name <- gsub("_baseline", "", scenario_var_name) # Remove tag
        # Split by the depth part to insert ssp and time tag correctly
        if (grepl("_depthsurf_", scenario_var_name)) {
          scenario_var_name <- sub("_depthsurf_", paste0("_", ssp_code, "_depthsurf_", time_tag_clean, "_"), scenario_var_name)
        } else if (grepl("_depthmax_", scenario_var_name)) {
          scenario_var_name <- sub("_depthmax_", paste0("_", ssp_code, "_depthmax_", time_tag_clean, "_"), scenario_var_name)
        } else {
          warning("Cannot determine depth part for future name generation in: ", core_var, call. = FALSE)
          return(NA) # Failed conversion
        }
      }
      return(scenario_var_name)
      
    } else {
      warning("Unknown scenario type: ", scenario_name, call. = FALSE)
      return(NA) # Indicate failure
    }
  }, USE.NAMES = FALSE) # Ensure output is unnamed vector
  
  return(na.omit(scenario_vars)) # Remove any NAs from failed conversions
}


#' Plot VIF results using ggplot (Corrected Argument Handling)
#' Uses car::vif output. Applies display name lookup via config object.
#'
#' @param vif_result Named numeric vector (output from car::vif).
#' @param save_path Optional file path to save the plot.
#' @param config The project configuration list (must contain core_var_display_names lookup and get_display_name function).
#' @return A ggplot object, or NULL on error.
#' @export
plot_vif_results_original <- function(vif_result, save_path = NULL, config) {
  if(!is.numeric(vif_result) || is.null(names(vif_result)) || length(vif_result) == 0){
    warning("plot_vif_results_original expects a non-empty named numeric vector.", call.=FALSE); return(NULL)
  }
  if(is.null(config$core_var_display_names) || is.null(config$get_display_name) || !is.function(config$get_display_name)){
    warning("Config object must contain 'core_var_display_names' lookup and 'get_display_name' function.", call.=FALSE)
    # Proceed using technical names as fallback
    display_lookup <- NULL
    get_display_name_func <- function(name, ...) name # Simple fallback
  } else {
    display_lookup <- config$core_var_display_names
    get_display_name_func <- config$get_display_name
  }
  
  original_names <- names(vif_result)
  display_names <- sapply(original_names, get_display_name_func, lookup = display_lookup, USE.NAMES = FALSE)
  
  df <- data.frame(
    OriginalName = original_names,
    DisplayName = display_names,
    VIF = as.numeric(vif_result)
  )
  num_vars <- nrow(df); if (num_vars == 0) return(NULL)
  df$VIF <- pmax(0, df$VIF)
  df$DisplayName <- factor(df$DisplayName, levels = df$DisplayName[order(df$VIF)])
  
  p <- ggplot2::ggplot(df, aes(x = DisplayName, y = VIF)) +
    ggplot2::geom_bar(stat = "identity", aes(fill = OriginalName), show.legend = FALSE) +
    ggplot2::scale_fill_viridis_d(option = "plasma") +
    ggplot2::labs(
      title = "VIF Analysis Results (car::vif)",
      x = "Environmental Variables",
      y = "VIF Value"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, max(config$vif_threshold, ceiling(max(df$VIF, na.rm = TRUE)))),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    ggplot2::geom_hline(yintercept = config$vif_threshold, linetype = "dashed", color = "red") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (!is.null(save_path)) {
    tryCatch({
      ggplot2::ggsave(save_path, plot = p, width = max(8, 4 + 0.2 * num_vars), height = 6, limitsize = FALSE)
      cat("  VIF plot saved to", save_path, "\n")
    }, error = function(e) { warning("Failed to save VIF plot: ", e$message, call.=FALSE)})
  }
  return(p)
}


#' Plot Pearson Correlation using corrplot (Corrected Argument Handling)
#' Uses Pearson correlation results. Applies display name lookup via config object.
#'
#' @param env_extract Data frame or matrix of environmental values.
#' @param save_path Optional file path to save the plot.
#' @param config The project configuration list (must contain core_var_display_names lookup and get_display_name function).
#' @return NULL (invisibly). Plot is saved or printed.
#' @export
plot_correlation_results_original <- function(env_extract, save_path = NULL, config) {
  if(!is.data.frame(env_extract) && !is.matrix(env_extract) || ncol(env_extract) < 2) {
    warning("plot_correlation_results_original requires df/matrix with >= 2 columns.", call.=FALSE); return(invisible(NULL))
  }
  if(is.null(config$core_var_display_names) || is.null(config$get_display_name) || !is.function(config$get_display_name)){
    warning("Config object must contain 'core_var_display_names' lookup and 'get_display_name' function.", call.=FALSE)
    display_lookup <- NULL # Fallback
    get_display_name_func <- function(name, ...) name
  } else {
    display_lookup <- config$core_var_display_names
    get_display_name_func <- config$get_display_name
  }
  
  env.cor <- tryCatch(stats::cor(env_extract, method = "pearson", use = "pairwise.complete.obs"),
                      error = function(e) {warning("Correlation calculation failed: ", e$message, call.=FALSE); NULL})
  if(is.null(env.cor)) return(invisible(NULL))
  
  env.p <- NULL
  if(requireNamespace("Hmisc", quietly = TRUE)){
    numeric_cols <- sapply(env_extract, is.numeric)
    if(sum(numeric_cols) < 2) { warning("Need >= 2 numeric columns for Hmisc::rcorr.", call.=FALSE) }
    else {
      env_extract_num <- as.matrix(env_extract[, numeric_cols, drop = FALSE])
      n_obs <- crossprod(!is.na(env_extract_num))
      if(any(n_obs < 3)) { warning("Insufficient non-NA pairs, p-values from Hmisc may be unreliable.", call.=FALSE) }
      hmisc_result <- tryCatch(Hmisc::rcorr(env_extract_num, type="pearson"), error=function(e){warning("Hmisc::rcorr failed: ",e$message, call.=FALSE); NULL})
      if(!is.null(hmisc_result)) env.p <- hmisc_result$P
    }
  } else { warning("Hmisc package not found, p-values missing.", call.=FALSE) }
  
  original_names <- colnames(env.cor)
  display_names <- sapply(original_names, get_display_name_func, lookup = display_lookup, USE.NAMES = FALSE)
  rownames(env.cor) <- display_names
  colnames(env.cor) <- display_names
  
  if (!is.null(env.p) && all(dim(env.p) == dim(env.cor))) {
    p_original_names <- colnames(env_extract)[numeric_cols]
    p_display_names <- sapply(p_original_names, get_display_name_func, lookup = display_lookup, USE.NAMES = FALSE)
    rownames(env.p) <- p_display_names
    colnames(env.p) <- p_display_names
    # Ensure correct order/subset for p.mat if zero-var cols were removed
    valid_display_names_for_p <- intersect(display_names, p_display_names)
    env.p <- env.p[valid_display_names_for_p, valid_display_names_for_p]
    env.cor_subset <- env.cor[valid_display_names_for_p, valid_display_names_for_p] # Use subset for plotting if p reduced
  } else {
    env.p <- NULL # Disable p-values if mismatch
    env.cor_subset <- env.cor # Plot full matrix
  }
  
  
  # --- Plotting ---
  tryCatch({
    plot_obj <- function() {
      corrplot::corrplot(
        corr = env.cor_subset,      # Use potentially subsetted matrix
        method = "color",
        type = "upper",
        order = "hclust",
        p.mat = env.p,           # Use potentially subsetted/reordered p-matrix
        sig.level = c(.01, .05),
        insig = "label_sig",
        pch.cex = 1.5,
        pch.col = "grey",
        #addCoef.col = "black", # Can make plot crowded
        #number.cex = 0.7,
        tl.col = "black",
        tl.srt = 45,
        diag = FALSE,
        na.label = "NA",
        mar = c(0, 0, 1, 0),
        title = "Pearson Correlation Matrix"
      )
    }
    
    if (!is.null(save_path)) {
      plot_dim <- max(8, 4 + 0.3 * ncol(env.cor_subset))
      grDevices::png(filename = save_path, width = plot_dim, height = plot_dim, units = "in", res = 300)
      plot_obj()
      grDevices::dev.off()
      cat("  Correlation plot saved to", save_path, "\n")
    } else {
      plot_obj()
    }
    return(invisible(NULL))
    
  }, error = function(e){
    warning("Failed to create or save correlation plot: ", e$message, call.=FALSE)
    return(NULL)
  })
}