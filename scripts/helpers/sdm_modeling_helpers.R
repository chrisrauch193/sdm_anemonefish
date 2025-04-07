# scripts/helpers/sdm_modeling_helpers.R
#-------------------------------------------------------------------------------
# Helper functions for running individual species SDMs using SDMtune
# Reverted internal messages to cat/warning for parallel simplicity.
#-------------------------------------------------------------------------------
pacman::p_load(dplyr, sf, terra, dismo, SDMtune, readr, tools, stringr)

# --- load_clean_individual_occ ---
load_clean_individual_occ <- function(species_aphia_id, occurrence_folder, config) {
  occurrence_file <- file.path(occurrence_folder, paste0(species_aphia_id, ".csv"))
  if (!file.exists(occurrence_file)) {
    warning("Occurrence file not found for AphiaID ", species_aphia_id, " in ", occurrence_folder, call. = FALSE); return(NULL)
  }
  tryCatch({ occurrences <- readr::read_csv(occurrence_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  }, error = function(e){ warning("Error reading occurrence file ", occurrence_file, ": ", e$message, call.=FALSE); return(NULL) })
  if (!all(c("decimalLongitude", "decimalLatitude") %in% names(occurrences))) {
    warning("Required columns 'decimalLongitude'/'decimalLatitude' not found in ", occurrence_file, call.=FALSE); return(NULL)
  }
  occ_clean <- occurrences %>% mutate(decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)), decimalLatitude = suppressWarnings(as.numeric(decimalLatitude))) %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude), decimalLongitude >= -180, decimalLongitude <= 180, decimalLatitude >= -90, decimalLatitude <= 90)
  if (nrow(occ_clean) == 0) return(NULL)
  tryCatch({ occ_sf <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = config$occurrence_crs); return(occ_sf)
  }, error = function(e){ warning("Error converting occurrences to sf for AphiaID ", species_aphia_id, ": ", e$message, call.=FALSE); return(NULL) })
}

# --- thin_individual_occ ---
thin_individual_occ <- function(occurrences_sf, predictor_stack, config) {
  # (Keep implementation the same, warnings are okay here)
  if (is.null(occurrences_sf) || nrow(occurrences_sf) == 0 || is.null(predictor_stack)) return(NULL)
  target_crs_obj <- terra::crs(predictor_stack)
  if (sf::st_crs(occurrences_sf) != target_crs_obj) {
    tryCatch({ occurrences_sf <- sf::st_transform(occurrences_sf, crs = target_crs_obj)
    }, error = function(e){ warning("Thinning Error: Failed to transform CRS: ", e$message, call.=FALSE); return(NULL) })
  }
  tryCatch({
    occs_spatvector <- terra::vect(occurrences_sf); occs_cells <- terra::extract(predictor_stack[[1]], occs_spatvector, cells = TRUE)
    valid_cells_indices <- which(!is.na(occs_cells$cell)); if(length(valid_cells_indices) == 0) return(NULL)
    unique_cell_indices_in_valid <- which(!duplicated(occs_cells$cell[valid_cells_indices]))
    original_indices_to_keep <- valid_cells_indices[unique_cell_indices_in_valid]
    occs_thinned <- occurrences_sf[original_indices_to_keep, ]; if(nrow(occs_thinned) == 0) return(NULL)
    return(occs_thinned)
  }, error = function(e){ warning("Error during individual cell-based thinning: ", e$message, call.=FALSE); return(NULL) })
}

# --- generate_sdm_background ---
generate_sdm_background <- function(predictor_stack, n, config) {
  # Using cat for info, warning for issues
  cat("    Generating", n, "background points...\n") # Indented cat for loop context
  if (is.null(predictor_stack)) return(NULL)
  tryCatch({
    mask_rast <- !is.na(predictor_stack[[1]]); bg_points_terra <- terra::spatSample(mask_rast, size = n, method = "random", na.rm = TRUE, xy = TRUE, warn=FALSE)
    n_gen <- nrow(bg_points_terra)
    if(n_gen < n) warning("Only generated ", n_gen, " background points (less than requested ", n, ").", call.=FALSE)
    if(n_gen == 0) { warning("Failed to generate any background points.", call.=FALSE); return(NULL) }
    bg_df <- as.data.frame(bg_points_terra[, c("x", "y")]); return(bg_df)
  }, error = function(e) { warning("Error generating background points: ", e$message, call. = FALSE); return(NULL) })
}

# --- run_sdm_sdmtune_grid ---
run_sdm_sdmtune_grid <- function(occs_thinned_sf, predictor_stack, background_df, config) {
  cat("    Running SDMtune::gridSearch...\n") # Indented cat
  if (is.null(occs_thinned_sf) || nrow(occs_thinned_sf) == 0 || is.null(predictor_stack) || is.null(background_df)) {
    warning("Invalid inputs for running sdmtune.", call. = FALSE); return(NULL)
  }
  occs_coords <- sf::st_coordinates(occs_thinned_sf)
  swd_data <- tryCatch({ sdmtune::prepareSWD(species = "species", p = occs_coords, a = background_df, env = predictor_stack)
  }, error = function(e) { warning("Failed to prepare SWD object: ", e$message, call. = FALSE); return(NULL) })
  if (is.null(swd_data)) return(NULL)
  model <- sdmtune::Maxnet()
  hyper_grid <- config$sdm_tune_grid
  tuned_model <- NULL
  tryCatch({
    tuned_model <- sdmtune::gridSearch(swd_data, hypers = hyper_grid, metric = config$sdm_evaluation_metric,
                                       method = config$sdm_partitions, k = config$sdm_n_folds, save_models = TRUE)
    cat("    SDMtune::gridSearch completed.\n")
  }, error = function(e) { warning("SDMtune::gridSearch failed: ", e$message, call. = FALSE); tuned_model <- NULL })
  return(tuned_model)
}

# --- predict_sdm_sdmtune ---
predict_sdm_sdmtune <- function(sdmtune_results, predictor_stack, config) {
  cat("    Selecting best model and predicting using sdmtune...\n") # Indented cat
  if (is.null(sdmtune_results) || !inherits(sdmtune_results, "SDMtune")) {
    warning("Invalid sdmtune results object provided.", call. = FALSE); return(NULL)
  }
  if (is.null(predictor_stack)) { warning("Predictor stack required.", call. = FALSE); return(NULL) }
  best_model_obj <- tryCatch({
    if(length(sdmtune_results@models) > 0) {
      res_df <- sdmtune::results(sdmtune_results); metric <- config$sdm_evaluation_metric
      metric_cv_name <- paste0(metric, "_cv"); if (!metric_cv_name %in% colnames(res_df)) {
        if("AUC_cv" %in% colnames(res_df)) metric_cv_name <- "AUC_cv" else if ("TSS_cv" %in% colnames(res_df)) metric_cv_name <- "TSS_cv"
        else {warning("Cannot find metric '{metric}_cv' or fallbacks. Using first model.", call.=FALSE); return(sdmtune_results@models[[1]])}
        warning("Metric '{config$sdm_evaluation_metric}_cv' not found, using '{metric_cv_name}' instead.", call.=FALSE)
      }
      best_row_index <- which.max(res_df[[metric_cv_name]]); if(length(best_row_index) == 0 || best_row_index > length(sdmtune_results@models)) {warning("Could not determine best model index. Using first model.", call.=FALSE); best_row_index <- 1}
      cat("      Best model hyperparameters:", paste(names(sdmtune_results@hypers), sdmtune_results@hypers[best_row_index,], collapse=", "), "\n")
      sdmtune_results@models[[best_row_index]]
    } else {warning("No models found within the sdmtune results object.", call.=FALSE); NULL}
  }, error = function(e) {warning("Error accessing best model from sdmtune results: {e$message}", call. = FALSE); NULL})
  if (is.null(best_model_obj)) {warning("Could not retrieve best model object.", call.=FALSE); return(NULL)}
  prediction_raster <- NULL
  tryCatch({
    prediction_raster <- sdmtune::predict(object = best_model_obj, data = predictor_stack, type = "cloglog")
    names(prediction_raster) <- "suitability"; cat("    Prediction raster generated.\n")
  }, error = function(e) {warning("sdmtune::predict failed: {e$message}", call. = FALSE); prediction_raster <- NULL})
  return(prediction_raster)
}
#-------------------------------------------------------------------------------