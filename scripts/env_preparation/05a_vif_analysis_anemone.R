# scripts/env_preparation/05a_vif_analysis_anemone.R
#-------------------------------------------------------------------------------
# INTERACTIVE VIF/Correlation Analysis Script for ANEMONE HOSTS
# Analyzes a user-defined CORE set of variables across ALL scenarios.
#-------------------------------------------------------------------------------
cat("--- Running Script 05a: VIF/Correlation Analysis for Anemone Hosts (All Scenarios) ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, corrplot, ggplot2, car, tools, Hmisc, stringr) # Load necessary packages

# Source helpers (ensure generate_scenario_variable_list and plotting functions are loaded)
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R"); source(env_helper_path)
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R")) # For occ loading

# --- 2. Define Group ---
group_name <- "anemone"
group_occ_dir <- config$anemone_occurrence_dir
cat("Analyzing Group:", group_name, "\n")

# --- 3. *** USER DEFINED: CORE VARIABLE LIST (Anemone) *** ---
# Define the CORE variables using BASELINE format (or simplest form for terrain).
# This list will be used to generate names for ALL scenarios.
# Comment/uncomment lines here to iterate VIF/correlation analysis.
# The goal is to arrive at ONE list that works well across scenarios (esp. 'current').
#----------------------------------------------------------------
core_vars_anemone <- c(
  # Temperature
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_range",
  #"thetao_baseline_depthmax_ltmin",
  #"thetao_baseline_depthmax_ltmax",
  #"thetao_baseline_depthsurf_mean",
  # Salinity
  "so_baseline_depthmax_mean",
  # Oxygen
  #"o2_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  #"o2_baseline_depthmax_ltmin",
  #"o2_baseline_depthmax_ltmax",
  # Nitrate
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  #"no3_baseline_depthmax_ltmin",
  #"no3_baseline_depthmax_ltmax",
  # Biology/Light
  "chl_baseline_depthmax_mean",
  #"phyc_baseline_depthmax_mean",
  "par_baseline_depthsurf_mean",
  # Physical
  "sws_baseline_depthsurf_mean",
  #"ph_baseline_depthmax_mean",
  # Terrain (Names should match files in terrain folder)
  "bathymetry_mean",
  "slope",
  # "rugosity",
  "distcoast"
)

core_vars_anemone <- c(
  "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_range",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  "chl_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  # Terrain Vars
  "bathymetry_mean", "distcoast", "rugosity"
)
cat("--- Core variables selected for analysis:", paste(core_vars_anemone, collapse=", "), "---\n")
# --- END USER INTERACTION ---
#-------------------------------------------------------------------------------


# --- 4. Analysis Loop: Scenario -> VIF/Correlation ---
vif_results_all_scenarios <- list() # Store VIF results per scenario

for (scenario in config$env_scenarios) {
  cat("\n================ Processing Scenario:", scenario, " ================\n")
  
  # Generate scenario-specific technical names from the core list
  current_vars_to_test <- generate_scenario_variable_list(core_vars_anemone, scenario, config)
  if(length(current_vars_to_test) < 2) { warning("Need >= 2 variables for scenario '", scenario, "'. Skipping.", call. = FALSE); next }
  cat("  Attempting to analyze:", paste(current_vars_to_test, collapse=", "), "\n")
  
  # Load FULL stack for the scenario
  env_stack_full_raw <- load_stack_env_data(scenario, config)
  if (is.null(env_stack_full_raw)) { warning("Load stack failed for scenario '", scenario, "'. Skipping.", call. = FALSE); next }
  
  # Check if ALL generated variables exist in the loaded stack
  available_layers <- names(env_stack_full_raw)
  missing_in_stack <- setdiff(current_vars_to_test, available_layers)
  if(length(missing_in_stack) > 0) {
    warning("Following generated variables MISSING from stack for scenario '", scenario, "':\n  ",
            paste(missing_in_stack, collapse=", "), "\nCheck naming convention or available files. Skipping scenario.", call. = FALSE)
    print(paste("Available layers in stack:", paste(available_layers, collapse=", ")))
    rm(env_stack_full_raw); gc(); next
  }
  
  # Subset stack to ONLY the variables generated for this scenario
  env_stack_subset <- env_stack_full_raw[[current_vars_to_test]]
  rm(env_stack_full_raw); gc()
  
  # Preprocess the SUBSETTED stack
  env_stack_processed <- preprocess_env_rasters(env_stack_subset, config)
  if (is.null(env_stack_processed)) { warning("Preprocess failed for scenario '", scenario, "'. Skipping.", call. = FALSE); rm(env_stack_subset); gc(); next }
  rm(env_stack_subset); gc(); target_crs <- terra::crs(env_stack_processed)
  
  # Load/Clean/Thin Occurrences (Thinning uses the processed stack for extent/resolution)
  group_occs_sf <- load_clean_thin_group_occurrences(group_occ_dir, target_crs, config, env_stack_processed)
  if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) { warning("No valid occurrences for group '", group_name, "'. Skipping VIF/Corr.", call. = FALSE); next }
  
  # Extract Env Values
  env_values_df <- extract_env_values(group_occs_sf, env_stack_processed)
  if (is.null(env_values_df) || nrow(env_values_df) < 2 || ncol(env_values_df) < 2) { warning("Insufficient env data extracted for scenario '", scenario, "'. Skipping.", call. = FALSE); next }
  
  # Check Zero Variance
  variances <- apply(env_values_df, 2, var, na.rm = TRUE)
  zero_var_cols <- names(variances[variances == 0 | is.na(variances)])
  if (length(zero_var_cols) > 0) {
    cat("  Warning: Removing zero-variance columns post-extraction:", paste(zero_var_cols, collapse=", "), "\n")
    env_values_df <- env_values_df[, !(names(env_values_df) %in% zero_var_cols), drop = FALSE]
  }
  if (ncol(env_values_df) < 2) { warning("Less than 2 variables remaining after NA/zero-var removal for scenario '", scenario, "'. Skipping.", call. = FALSE); next }
  final_vars_analyzed <- names(env_values_df) # Keep track of exactly what was analyzed
  cat("  Variables actually used in VIF/Corr for this scenario:", paste(final_vars_analyzed, collapse=", "), "\n")
  
  # Define Output Paths for this specific run
  run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") # Unique timestamp
  vif_plot_path <- file.path(config$log_dir, paste0("vif_plot_", group_name, "_", scenario, "_", run_timestamp, ".png"))
  cor_plot_path <- file.path(config$log_dir, paste0("correlation_plot_", group_name, "_", scenario, "_", run_timestamp, ".png"))
  
  # Run VIF Analysis
  cat("\n  Running VIF analysis for scenario:", scenario, "...\n")
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
    vif_results_all_scenarios[[scenario]] <- vif_results # Store results
    cat("  VIF values:\n"); print(sort(vif_results, decreasing=TRUE))
    # Pass the config object which contains the get_display_name function and lookup
    plot_vif_results_original(vif_results, save_path = vif_plot_path, config = config)
    high_vif_vars <- names(vif_results[vif_results > config$vif_threshold])
    if(length(high_vif_vars) > 0) { cat("  ---> VIF >", config$vif_threshold, ":", paste(high_vif_vars, collapse=", "), "\n") } else { cat("  All VIF <= ", config$vif_threshold, ".\n") }
  }, error = function(e) { warning("VIF failed for scenario '", scenario, "': ", e$message, call.=FALSE) })
  
  # Run Correlation Analysis
  cat("\n  Running Pearson Correlation for scenario:", scenario, "...\n")
  tryCatch({
    # Pass the config object
    plot_correlation_results_original(env_values_df, save_path = cor_plot_path, config = config)
    cor_matrix <- cor(env_values_df, method = "pearson", use = "pairwise.complete.obs"); cor_matrix_vals <- cor_matrix
    cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA
    high_cor_pairs_indices <- which(abs(cor_matrix) > config$correlation_threshold, arr.ind = TRUE)
    if(nrow(high_cor_pairs_indices) > 0){
      cat("  ---> Highly correlated pairs (|r| >", config$correlation_threshold, "):\n")
      # Use display names via config function
      display_names_in_matrix <- sapply(rownames(cor_matrix_vals), config$get_display_name, lookup = config$core_var_display_names)
      for(k in 1:nrow(high_cor_pairs_indices)) {
        row_idx <- high_cor_pairs_indices[k, 1]; col_idx <- high_cor_pairs_indices[k, 2]
        var1_display <- display_names_in_matrix[row_idx]; var2_display <- display_names_in_matrix[col_idx]
        cor_value <- cor_matrix_vals[row_idx, col_idx]
        cat("    - ", var1_display, " & ", var2_display, " (r = ", round(cor_value, 2), ")\n", sep="")
      }
    } else { cat("  No highly correlated pairs found (|r| <= ", config$correlation_threshold, ").\n") }
  }, error = function(e) { warning("Correlation failed for scenario '", scenario, "': ", e$message, call.=FALSE) })
  
  cat("\n  Analysis complete for scenario '", scenario, "'.\n", sep="")
  rm(env_stack_processed, group_occs_sf, env_values_df, vif_results, cor_matrix, cor_matrix_vals); gc()
} # End scenario loop

cat("\n--- Script 05a finished. Review logs/plots for all scenarios. ---\n")
cat("--- Decide on ONE FINAL set based primarily on 'current' scenario VIF/Corr. ---\n")
cat("--- Check if that final set is acceptable in future scenarios. ---\n")
cat("--- Use that FINAL set for SDM runs in script 06. ---\n")
#-------------------------------------------------------------------------------