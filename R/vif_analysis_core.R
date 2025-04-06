# R/vif_analysis_core.R

# This script defines the core function for VIF analysis.
# It should be called by scenario-specific driver scripts.

library(raster)
library(usdm) # Though car::vif is used below, usdm might be needed elsewhere
library(terra)
library(dplyr)
library(tools)
library(car) # Specifically for vif

# Source necessary helpers
source("helpers/extract_env.R") # Assumes extract_env.R is in helpers/
source("helpers/plotting_helpers.R") # Assumes plotting_helpers.R is in helpers/

#' Run VIF and Correlation Analysis for SDM Variable Selection
#'
#' Loads environmental and occurrence data, extracts environmental values
#' at occurrence points, performs VIF and correlation analysis, and saves plots.
#'
#' @param env_folder Path to the folder containing environmental raster layers.
#' @param terrain_folder Path to the folder containing terrain raster layers.
#' @param coral_shapefile Path to the coral reef shapefile for potential masking.
#' @param occurrence_folder Path to the folder containing occurrence CSV files.
#' @param initial_selected_vars A character vector of technical variable names
#'   to include in the initial analysis.
#' @param save_location Base directory where results (plots, logs) will be saved.
#' @param output_prefix A string prefix (e.g., "anemone", "anemonefish_future")
#'   used for naming output files and creating a subdirectory within save_location.
#' @param occurrence_crs CRS string (e.g., "EPSG:4326") for occurrence data.
#' @param env_crs CRS string (e.g., "EPSG:4326") for environmental data.
#' @param vif_threshold The VIF threshold value for interpretation (used for plotting line). Default 10.
#'
#' @return Invisibly returns the calculated VIF results vector.
#' @export
#'
run_vif_analysis <- function(
    env_folder,
    terrain_folder,
    coral_shapefile,
    occurrence_folder,
    initial_selected_vars,
    save_location,
    output_prefix,
    occurrence_crs,
    env_crs,
    vif_threshold = 10
) {
  
  cat(paste("--- Starting VIF Analysis for:", output_prefix, "---\n"))
  
  # --- Setup Output Directory ---
  output_dir <- file.path(save_location, output_prefix)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  } else {
    cat("Output directory:", output_dir, "\n")
  }
  
  # --- 1. Load and Prepare Environmental Data ---
  cat("Loading and preparing environmental data...\n")
  
  # Load all available data first
  env_stack_full <- load_and_prepare_env_data(env_folder, terrain_folder, coral_shapefile)
  if (is.null(env_stack_full)) {
    stop("Failed to load and prepare environmental data.")
  }
  
  # Select the desired variables for this analysis run
  env_stack_selected <- env_stack_full[[initial_selected_vars]]
  cat("Selected", nlyr(env_stack_selected), "variables for analysis.\n")
  
  print("load_occurrence_and_extract_env!!!!")
  # --- 2. Load Occurrence Data and Extract Env Values ---
  cat("Loading occurrence data and extracting environmental values...\n")
  env_values_df <- load_occurrence_and_extract_env(
    occurrence_folder = occurrence_folder,
    env_stack = env_stack_selected,
    occurrence_crs = occurrence_crs
  )
  print("ENDING load_occurrence_and_extract_env!!!!")
  
  
  if (is.null(env_values_df) || nrow(env_values_df) == 0) {
    stop("Failed to load occurrence data or extract environmental values. No valid data points found.")
  }
  cat("Extracted values for", nrow(env_values_df), "occurrence points.\n")
  
  # --- 3. Collinearity Analysis ---
  cat("Performing collinearity analysis...\n")
  
  # Remove rows with any NAs first - important for VIF/cor
  env_values_df_clean <- na.omit(env_values_df)
  rows_removed <- nrow(env_values_df) - nrow(env_values_df_clean)
  if (rows_removed > 0) {
    cat("Removed", rows_removed, "rows with NA values.\n")
  }
  if (nrow(env_values_df_clean) < 2) {
    stop("Less than 2 complete observations remain after NA removal. Cannot perform VIF/correlation.")
  }
  if (nrow(env_values_df_clean) < ncol(env_values_df_clean)) {
    warning("Number of observations (", nrow(env_values_df_clean),
            ") is less than the number of variables (", ncol(env_values_df_clean),
            "). VIF results might be unstable.")
  }
  
  # Remove columns with zero variance (can happen after subsetting/NA removal)
  variances <- apply(env_values_df_clean, 2, var, na.rm = TRUE)
  zero_var_cols <- variances == 0 | is.na(variances)
  if (any(zero_var_cols)) {
    removed_cols <- names(variances[zero_var_cols])
    warning("Removing columns with zero or NA variance before VIF/correlation: ",
            paste(removed_cols, collapse=", "))
    env_values_df_analysis <- env_values_df_clean[, !zero_var_cols, drop = FALSE]
    # Update the list of variables used
    final_selected_vars <- colnames(env_values_df_analysis)
    if (ncol(env_values_df_analysis) < 2) {
      stop("Less than 2 variables remaining after removing zero variance columns. Cannot perform VIF/correlation.")
    }
  } else {
    env_values_df_analysis <- env_values_df_clean
    final_selected_vars <- colnames(env_values_df_analysis) # These are the vars actually used
  }
  
  
  # -- Correlation Analysis --
  cat("Calculating and plotting Pearson correlation...\n")
  save_path_cor <- file.path(output_dir, paste0(output_prefix, "_pearson_cor.png"))
  correlation_matrix <- plot_correlation_results(
    env_extract = env_values_df_analysis,
    save_path = save_path_cor
  )
  
  # -- VIF Analysis --
  # Create dummy response variable needed for lm() used by car::vif
  # Using 1 for all rows assumes presence-only data context for VIF
  dummy_response <- rep(1, nrow(env_values_df_analysis))
  vif_data <- cbind(dummy_response, env_values_df_analysis)
  
  # Fit linear model (response doesn't matter for VIF calculation itself)
  # Use formula interface for robustness
  formula_str <- paste("dummy_response ~", paste(colnames(env_values_df_analysis), collapse = " + "))
  lm_result <- lm(as.formula(formula_str), data = vif_data)
  
  # Calculate VIF using car::vif
  cat("Calculating VIF values...\n")
  vif_result <- tryCatch({
    car::vif(lm_result)
  }, error = function(e) {
    warning("VIF calculation failed. Error: ", e$message, ". Skipping VIF plot.")
    return(NULL)
  })
  
  # Plot VIF results
  if (!is.null(vif_result)) {
    cat("Plotting VIF results...\n")
    save_path_vif <- file.path(output_dir, paste0(output_prefix, "_vif_analysis.png"))
    # Ensure only names from vif_result are passed to lookup subset
    vif_var_names <- names(vif_result)
    
    plot_vif_results(
      vif_result = vif_result,
      save_path = save_path_vif
    )
    print("VIF Values:")
    print(vif_result)
  }
  
  
  # --- 4. Save List of Variables Used ---
  # Save the list of variables *actually used* in the final VIF/correlation
  # after potential removals (NA, zero variance)
  variables_file <- file.path(output_dir, paste0(output_prefix, "_variables_used.txt"))
  writeLines(final_selected_vars, variables_file)
  cat("\nList of variables actually used in the analysis saved to:", variables_file, "\n")
  
  cat(paste("--- VIF Analysis for:", output_prefix, "Complete ---\n"))
  
  # Return the VIF results invisibly
  invisible(vif_result)
}