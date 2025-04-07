# scripts/05c_vif_analysis_anemonefish_combined.R
#-------------------------------------------------------------------------------
# INTERACTIVE VIF/Correlation Analysis Script for ANEMONEFISH (Combined Env + Host)
# Analyzes collinearity for the 'current' scenario including host predictor(s).
#-------------------------------------------------------------------------------
cat("--- Running Script 05c: VIF/Correlation Analysis for Anemonefish (Combined Env + Host - Current Only) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, corrplot, ggplot2, car, tools, Hmisc, stringr)
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R"); source(env_helper_path)
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))

# --- 2. Define Group & Scenario ---
group_name <- "anemonefish"
group_occ_dir <- config$anemonefish_occurrence_dir
scenario <- "current" # **** ONLY analyze current scenario for combined model ****
cat("Analyzing Group:", group_name, "(Combined Env + Host) for Scenario:", scenario, "\n")

# --- 3. *** USER DEFINED: CORE ENVIRONMENTAL VARIABLE LIST (Anemonefish Combined) *** ---
# Define the CORE environmental variables (baseline format) to test alongside the host predictor.
# This list should ideally be the FINAL non-collinear set identified from script 05b for the 'current' scenario.
# Comment/uncomment lines here to iterate VIF/correlation *for the combined model*.
#----------------------------------------------------------------
core_env_vars_fish_combined <- c(
  # --- PASTE Your Final Non-Collinear ENV VARS from 05b 'current' analysis ---
  # Example:
  "thetao_baseline_depthmax_mean",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "chl_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
)
cat("--- Core ENVIRONMENTAL variables selected for combined analysis:", paste(core_env_vars_fish_combined, collapse=", "), "---\n")
# --- END USER INTERACTION ---
#-------------------------------------------------------------------------------

# --- 4. Load Current Environmental Data & Subset ---
if(length(core_env_vars_fish_combined) < 1) stop("No environmental variables selected in core_env_vars_fish_combined.")
cat("  Loading and processing selected environmental variables for 'current'...\n")

env_stack_full_raw <- load_stack_env_data(scenario, config)
if(is.null(env_stack_full_raw)) stop("Failed to load current env stack.")

available_layers <- names(env_stack_full_raw)
missing_in_stack <- setdiff(core_env_vars_fish_combined, available_layers)
if(length(missing_in_stack) > 0) {
  stop("Selected env vars missing from stack: ", paste(missing_in_stack, collapse=", "))
}
env_stack_subset <- env_stack_full_raw[[core_env_vars_fish_combined]]
rm(env_stack_full_raw); gc()
env_stack_processed <- preprocess_env_rasters(env_stack_subset, config)
if(is.null(env_stack_processed)) stop("Failed to preprocess env subset.")
rm(env_stack_subset); gc(); target_crs <- terra::crs(env_stack_processed)

# --- 5. Load Anemonefish Occurrences ---
group_occs_sf <- load_clean_thin_group_occurrences(group_occ_dir, target_crs, config, env_stack_processed)
if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) {
  stop("No valid anemonefish occurrences found for VIF/Corr.")
}

# --- 6. Load and Combine Host Prediction Rasters ---
# (Same logic as before: load all anemone predictions for 'current', make max layer)
cat("  Loading ALL host anemone predictions for 'current' scenario...\n")
anemone_list_file <- config$anemone_species_list_file
if (!file.exists(anemone_list_file)) stop("Anemone list missing.")
anemone_species <- readr::read_csv(anemone_list_file, show_col_types = FALSE)
anemone_pred_dir <- config$predictions_dir

host_predictor_list <- list()
for (anem_idx in 1:nrow(anemone_species)) {
  anem_name <- anemone_species$scientificName[anem_idx]
  anem_name_sanitized <- gsub(" ", "_", anem_name)
  anem_pred_file <- file.path(anem_pred_dir, paste0("sdm_prediction_", anem_name_sanitized, "_", scenario, ".tif")) # Only current
  if (file.exists(anem_pred_file)) {
    pred_rast <- tryCatch(terra::rast(anem_pred_file), error = function(e) NULL)
    if (!is.null(pred_rast)) {
      if (!terra::compareGeom(env_stack_processed[[1]], pred_rast, stopOnError=FALSE)) {
        warning("Resampling host predictor ", basename(anem_pred_file), " to match env stack.", call.=FALSE)
        pred_rast <- terra::resample(pred_rast, env_stack_processed[[1]], method="bilinear")
      }
      names(pred_rast) <- paste0("host_", anem_name_sanitized)
      host_predictor_list[[anem_name_sanitized]] <- pred_rast
    } else { warning("Failed to load host pred: ", basename(anem_pred_file)) }
  } else { warning("Host pred file missing: ", basename(anem_pred_file)) }
}
if (length(host_predictor_list) == 0) stop("No host predictions loaded.")

host_stack <- terra::rast(host_predictor_list)
combined_host_suitability <- terra::app(host_stack, fun = "max", na.rm = TRUE)
names(combined_host_suitability) <- "host_suitability_max"
cat("  Created combined MAX host suitability layer.\n")

# Combine with environmental predictors
combined_stack <- c(env_stack_processed, combined_host_suitability)
cat("  Combined predictor stack created with layers:", paste(names(combined_stack), collapse=", "), "\n")

# --- 7. Extract Values using Combined Stack ---
env_values_df <- extract_env_values(group_occs_sf, combined_stack)
if (is.null(env_values_df) || nrow(env_values_df) < 2 || ncol(env_values_df) < 2) {
  stop("Insufficient data extracted using combined stack.")
}

# --- 8. Check Zero Variance ---
variances <- apply(env_values_df, 2, var, na.rm = TRUE)
zero_var_cols <- names(variances[variances == 0 | is.na(variances)])
if (length(zero_var_cols) > 0) {
  cat("  Warning: Removing zero-variance columns:", paste(zero_var_cols, collapse=", "), "\n")
  env_values_df <- env_values_df[, !(names(env_values_df) %in% zero_var_cols), drop = FALSE]
}
if (ncol(env_values_df) < 2) { stop("Less than 2 vars remaining.") }
cat("  Variables remaining for analysis:", paste(names(env_values_df), collapse=", "), "\n")

# --- 9. Run VIF/Correlation Analysis ---
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
vif_plot_path <- file.path(config$log_dir, paste0("vif_plot_", group_name, "_combined_", scenario, "_", run_timestamp, ".png"))
cor_plot_path <- file.path(config$log_dir, paste0("correlation_plot_", group_name, "_combined_", scenario, "_", run_timestamp, ".png"))

# Run VIF
cat("\n  Running VIF analysis (Combined)...\n")
vif_results <- NULL
tryCatch({
  present_dummy <- rep(1, nrow(env_values_df)); lm_data <- cbind(present_dummy, env_values_df)
  valid_colnames <- make.names(colnames(lm_data), unique = TRUE)
  original_colnames <- colnames(lm_data); colnames(lm_data) <- valid_colnames
  predictor_names_valid <- valid_colnames[-1]
  formula_str <- paste("present_dummy ~", paste(predictor_names_valid, collapse = " + "))
  full_model <- stats::lm(as.formula(formula_str), data = lm_data)
  vif_results_raw <- car::vif(full_model)
  vif_names_map <- setNames(original_colnames[-1], valid_colnames[-1])
  names(vif_results_raw) <- vif_names_map[names(vif_results_raw)]
  vif_results <- vif_results_raw
  cat("  VIF values:\n"); print(sort(vif_results, decreasing=TRUE))
  # Create temporary lookup including host for plotting this specific case
  temp_lookup_combined <- config$core_var_display_names
  temp_lookup_combined["host_suitability_max"] <- "Max Host Suitability"
  plot_vif_results_original(vif_results, save_path = vif_plot_path, display_lookup = temp_lookup_combined)
  high_vif_vars <- names(vif_results[vif_results > config$vif_threshold])
  if(length(high_vif_vars) > 0) { cat("  ---> VIF >", config$vif_threshold, ":", paste(high_vif_vars, collapse=", "), "\n") } else { cat("  All VIF <= ", config$vif_threshold, ".\n") }
}, error = function(e) { warning("VIF failed: ", e$message, call.=FALSE) })

# Run Correlation
cat("\n  Running Pearson Correlation (Combined)...\n")
tryCatch({
  temp_lookup_combined <- config$core_var_display_names
  temp_lookup_combined["host_suitability_max"] <- "Max Host Suitability"
  plot_correlation_results_original(env_values_df, save_path = cor_plot_path, display_lookup = temp_lookup_combined)
  cor_matrix <- cor(env_values_df, method = "pearson", use = "pairwise.complete.obs"); cor_matrix_vals <- cor_matrix
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA
  high_cor_pairs_indices <- which(abs(cor_matrix) > config$correlation_threshold, arr.ind = TRUE)
  if(nrow(high_cor_pairs_indices) > 0){
    cat("  ---> Highly correlated pairs (|r| >", config$correlation_threshold, "):\n")
    display_names_in_matrix <- sapply(rownames(cor_matrix_vals), config$get_display_name, lookup = temp_lookup_combined) # Use temp lookup
    for(k in 1:nrow(high_cor_pairs_indices)) {
      row_idx <- high_cor_pairs_indices[k, 1]; col_idx <- high_cor_pairs_indices[k, 2]
      var1_display <- display_names_in_matrix[row_idx]; var2_display <- display_names_in_matrix[col_idx]
      cor_value <- cor_matrix_vals[row_idx, col_idx]
      cat("    - ", var1_display, " & ", var2_display, " (r = ", round(cor_value, 2), ")\n", sep="")
    }
  } else { cat("  No highly correlated pairs found (|r| <= ", config$correlation_threshold, ").\n") }
}, error = function(e) { warning("Correlation failed: ", e$message, call.=FALSE) })

cat("\n  Analysis complete for COMBINED model variables.\n")
cat("  Review outputs:\n")
cat("    - VIF Plot:", vif_plot_path, "\n")
cat("    - Correlation Plot:", cor_plot_path, "\n")
cat("  Modify 'core_env_vars_fish_combined' list and re-run if needed.\n")
cat("  NOTE: This analysis only uses the 'current' scenario host predictions.\n")

cat("\n--- Script 05c finished. ---\n")
#-------------------------------------------------------------------------------