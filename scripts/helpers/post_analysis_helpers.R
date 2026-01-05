# scripts/helpers/post_analysis_helpers.R

#' Stack Species Prediction Rasters to Create a Richness Layer
#'
#' @param species_df Data frame of species (must have 'scientificName' column).
#' @param base_pred_dir Base directory for current predictions (e.g., config$target_predictions_current_dir).
#' @param scenario_name Name of the scenario (e.g., "current", "ssp119_2050").
#' @param model_suffix Suffix of the model type (e.g., "_pca", "_host_only").
#' @param config Project configuration.
#' @param log_func Logging function.
#' @return A SpatRaster object of summed richness, or NULL.
stack_richness <- function(species_df, base_pred_dir, scenario_name, model_suffix, config, log_func) {
  log_func("INFO", paste("Stacking richness for scenario:", scenario_name, "model_suffix:", model_suffix))
  pred_files <- c()
  for (i in 1:nrow(species_df)) {
    sp_name_sanitized <- gsub(" ", "_", species_df$scientificName[i])
    # Use the centralized helper to construct the expected path
    pred_file <- construct_prediction_filename(sp_name_sanitized, scenario_name, model_suffix, config)
    if (!is.null(pred_file) && file.exists(pred_file)) {
      pred_files <- c(pred_files, pred_file)
    } else {
      log_func("WARN", paste("Missing prediction for", sp_name_sanitized, "scenario", scenario_name, "suffix", model_suffix, "at path:", pred_file %||% "NULL_PATH"))
    }
  }
  if (length(pred_files) > 0) {
    richness_stack <- tryCatch(terra::rast(pred_files), error = function(e) {
      log_func("ERROR", paste("Failed to stack rasters:", e$message)); NULL
    })
    if (is.null(richness_stack)) return(NULL)
    richness_sum <- terra::sum(richness_stack, na.rm = TRUE)
    return(richness_sum)
  } else {
    log_func("WARN", paste("No prediction rasters found to stack for suffix", model_suffix, "scenario", scenario_name))
    return(NULL)
  }
}

#' Calculate Mean Difference in Suitability Between Two Rasters
#'
#' @param present_raster SpatRaster for the present period.
#' @param future_raster SpatRaster for the future period.
#' @param log_func Logging function.
#' @return Numeric mean difference, or NA.
calculate_suitability_change <- function(present_raster, future_raster, log_func) {
  if (!inherits(present_raster, "SpatRaster") || !inherits(future_raster, "SpatRaster")) {
    log_func("ERROR", "Invalid SpatRaster inputs for suitability change."); return(NA)
  }
  if (!terra::compareGeom(present_raster, future_raster, stopOnError = FALSE, res = TRUE)) {
    log_func("WARN", "Rasters have different geometries. Attempting resample for suitability change calc.")
    future_raster <- tryCatch(terra::resample(future_raster, present_raster), error = function(e) {
      log_func("ERROR", paste("Resampling failed for suitability change:", e$message)); NULL
    })
    if (is.null(future_raster)) return(NA)
  }
  diff_raster <- future_raster - present_raster
  mean_diff <- terra::global(diff_raster, "mean", na.rm = TRUE)
  return(as.numeric(mean_diff))
}

# Add other post-analysis helper functions here as needed...