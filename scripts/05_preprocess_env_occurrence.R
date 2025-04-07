# scripts/05_preprocess_env_occurrence.R
#-------------------------------------------------------------------------------
# INTERACTIVE VIF/Correlation Analysis Script
# Purpose: Analyze collinearity for a *user-defined* set of variables for a
#          specific group and scenario. Guides manual variable selection for SDMs.
# Removes PCA steps.
#-------------------------------------------------------------------------------
cat("--- Running Script 05: Interactive VIF/Correlation Analysis ---\n")

# --- 1. Setup ---
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, corrplot, ggplot2, car, tools, Hmisc) # Keep car, Hmisc

# Source helpers (ensure get_display_name is available)
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
if (!file.exists(env_helper_path)) stop("Env Helper missing.")
source(env_helper_path) # Make sure plot helpers and get_display_name are here or sourced elsewhere

# --- 2. Define Groups ---
species_groups <- list(
  anemone = list(occurrence_dir = config$anemone_occurrence_dir),
  anemonefish = list(occurrence_dir = config$anemonefish_occurrence_dir)
  # Add more groups if needed
)

# --- 3. *** USER INTERACTION REQUIRED HERE *** ---
# Define the variables to analyze FOR EACH SCENARIO.
# Use the CLEANED layer names (usually filename stems without .tif)
# Comment/uncomment lines within each scenario list to iterate.
#----------------------------------------------------
variables_to_analyze <- list()

# --- CURRENT ---
variables_to_analyze[['current']] <- c(
  # Temperature
  "thetao_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_range",
  "thetao_baseline_depthmax_ltmin",
  "thetao_baseline_depthmax_ltmax",
  # Salinity
  "so_baseline_depthmax_mean",
  # Oxygen
  "o2_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  "o2_baseline_depthmax_ltmin",
  "o2_baseline_depthmax_ltmax",
  # Nitrate
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  "no3_baseline_depthmax_ltmin",
  "no3_baseline_depthmax_ltmax",
  # Biology/Light
  "chl_baseline_depthmax_mean",
  "phyc_baseline_depthmax_mean",
  "par_baseline_depthsurf_mean",
  # Physical
  "sws_baseline_depthsurf_mean",
  "ph_baseline_depthmax_mean",
  # Terrain
  "bathymetry_mean",
  "slope",
  "rugosity",
  "distcoast"
)

# --- SSP119 2050 ---
variables_to_analyze[['ssp119_2050']] <- c(
  # Temperature
  # "thetao_ssp119_depthsurf_dec50_mean", # Add if exists and relevant
  "thetao_ssp119_depthmax_dec50_mean",
  "thetao_ssp119_depthmax_dec50_range",
  "thetao_ssp119_depthmax_dec50_ltmin",
  "thetao_ssp119_depthmax_dec50_ltmax",
  # Salinity
  "so_ssp119_depthmax_dec50_mean",
  # Oxygen
  # "o2_ssp119_depthmax_dec50_mean", # Add if exists and relevant
  "o2_ssp119_depthmax_dec50_range",
  "o2_ssp119_depthmax_dec50_ltmin",
  "o2_ssp119_depthmax_dec50_ltmax",
  # Nitrate
  "no3_ssp119_depthmax_dec50_mean",
  "no3_ssp119_depthmax_dec50_range",
  "no3_ssp119_depthmax_dec50_ltmin",
  "no3_ssp119_depthmax_dec50_ltmax",
  # Biology/Light (Copied)
  "chl_ssp119_depthmax_dec50_mean",
  "phyc_ssp119_depthmax_dec50_mean", # Assuming phyc was also copied/renamed
  "par_ssp119_depthsurf_dec50_mean",
  # Physical
  "sws_ssp119_depthsurf_dec50_mean",
  # "ph_ssp119_depthmax_dec50_mean", # Add if exists and relevant
  # Terrain
  "bathymetry_mean",
  "slope",
  "rugosity",
  "distcoast"
)

# --- SSP119 2100 ---
# Automatically generate from 2050 list if desired, or define manually
variables_to_analyze[['ssp119_2100']] <- gsub("_dec50_", "_dec100_", variables_to_analyze[['ssp119_2050']])
variables_to_analyze[['ssp119_2100']] <- unique(variables_to_analyze[['ssp119_2100']]) # Just in case

# --- SSP585 2050 ---
variables_to_analyze[['ssp585_2050']] <- gsub("ssp119", "ssp585", variables_to_analyze[['ssp119_2050']])
variables_to_analyze[['ssp585_2050']] <- unique(variables_to_analyze[['ssp585_2050']])

# --- SSP585 2100 ---
variables_to_analyze[['ssp585_2100']] <- gsub("ssp119", "ssp585", variables_to_analyze[['ssp119_2100']])
variables_to_analyze[['ssp585_2100']] <- unique(variables_to_analyze[['ssp585_2100']])

# --- END USER INTERACTION SECTION ---
#-------------------------------------------------------------------------------

# --- 4. Main Loop: Scenario -> Group ---
for (scenario in config$env_scenarios) {
  cat("\n=========================================================\n")
  cat("Processing Scenario:", scenario, "\n")
  cat("=========================================================\n")
  
  # --- Get the list of variables to test for this specific scenario ---
  current_vars_to_test <- variables_to_analyze[[scenario]]
  if (is.null(current_vars_to_test) || length(current_vars_to_test) < 2) {
    warning("Skipping scenario '", scenario, "': Variable list is missing or has < 2 variables.", call. = FALSE)
    next
  }
  cat("  Variables selected for analysis in this run:", paste(current_vars_to_test, collapse=", "), "\n")
  
  # --- Load FULL environmental stack first ---
  env_stack_full_raw <- load_stack_env_data(scenario, config)
  if (is.null(env_stack_full_raw)) {
    warning("Failed to load full env stack for scenario: ", scenario, ". Skipping scenario.", call. = FALSE)
    next
  }
  
  # --- Check if ALL selected variables exist in the loaded stack ---
  available_layers <- names(env_stack_full_raw)
  missing_in_stack <- setdiff(current_vars_to_test, available_layers)
  if(length(missing_in_stack) > 0) {
    warning("The following selected variables are MISSING from the loaded stack for scenario '", scenario, "':\n  ",
            paste(missing_in_stack, collapse=", "), "\nPlease check the variable list and filenames. Skipping scenario.", call. = FALSE)
    print(paste("Available layers:", paste(available_layers, collapse=", ")))
    rm(env_stack_full_raw); gc(); next
  }
  
  # --- Subset the stack to ONLY the selected variables ---
  env_stack_subset <- env_stack_full_raw[[current_vars_to_test]]
  rm(env_stack_full_raw); gc() # Free memory
  
  # --- Preprocess the SUBSETTED stack ---
  env_stack_processed <- preprocess_env_rasters(env_stack_subset, config)
  if (is.null(env_stack_processed)) {
    warning("Failed to preprocess env data for scenario: ", scenario, ". Skipping scenario.", call. = FALSE)
    rm(env_stack_subset); gc(); next
  }
  rm(env_stack_subset); gc() # Free memory
  target_crs <- terra::crs(env_stack_processed)
  
  # --- Loop through Groups ---
  for (group_name in names(species_groups)) {
    cat("\n-----------------------------------------------------\n")
    cat("Processing Group:", group_name, "for Scenario:", scenario, "\n")
    cat("-----------------------------------------------------\n")
    
    group_info <- species_groups[[group_name]]
    group_occ_dir <- group_info$occurrence_dir
    
    # Define UNIQUE output file paths for THIS RUN's analysis
    run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") # Add timestamp for uniqueness if iterating
    vif_plot_path <- file.path(config$log_dir, paste0("vif_plot_", group_name, "_", scenario, "_", run_timestamp, ".png"))
    cor_plot_path <- file.path(config$log_dir, paste0("correlation_plot_", group_name, "_", scenario, "_", run_timestamp, ".png"))
    analysis_log_path <- file.path(config$log_dir, paste0("analysis_log_", group_name, "_", scenario, "_", run_timestamp, ".txt"))
    
    # --- Load, Clean, Thin Occurrences ---
    # Use the processed stack (containing only selected vars) for thinning check
    group_occs_sf <- load_clean_thin_group_occurrences(group_occ_dir, target_crs, config, env_stack_processed)
    if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) {
      warning("No valid occurrences found for group: ", group_name, ". Skipping VIF/Corr for this combo.", call. = FALSE)
      next
    }
    
    # --- Extract Environmental Values (using the processed subset stack) ---
    env_values_df <- extract_env_values(group_occs_sf, env_stack_processed)
    if (is.null(env_values_df) || nrow(env_values_df) < 2 || ncol(env_values_df) < 2) {
      warning("Insufficient env data extracted for group: ", group_name, ". Rows:", nrow(env_values_df), " Cols:", ncol(env_values_df), ". Skipping VIF/Corr.", call. = FALSE)
      next
    }
    
    # --- Check for zero-variance columns AGAIN after extraction/NA removal ---
    variances <- apply(env_values_df, 2, var, na.rm = TRUE)
    zero_var_cols <- names(variances[variances == 0 | is.na(variances)]) # Include NA variance
    if (length(zero_var_cols) > 0) {
      cat("  Warning: Removing zero-variance columns after extraction:", paste(zero_var_cols, collapse=", "), "\n")
      env_values_df <- env_values_df[, !(names(env_values_df) %in% zero_var_cols), drop = FALSE]
    }
    if (ncol(env_values_df) < 2) {
      warning("Less than 2 variables remaining after removing zero variance post-extraction. Skipping VIF/Corr.", call. = FALSE)
      next
    }
    cat("  Variables remaining for analysis:", paste(names(env_values_df), collapse=", "), "\n")
    
    # --- Run VIF Analysis (using car::vif) ---
    cat("\n  Running VIF analysis (using lm and car::vif)...\n")
    vif_results <- NULL
    tryCatch({
      present_dummy <- rep(1, nrow(env_values_df))
      lm_data <- cbind(present_dummy, env_values_df)
      
      # Ensure column names are valid for formulas
      valid_colnames <- make.names(colnames(lm_data), unique = TRUE)
      original_colnames <- colnames(lm_data)
      colnames(lm_data) <- valid_colnames
      predictor_names_valid <- valid_colnames[-1] # Exclude dummy response
      
      formula_str <- paste("present_dummy ~", paste(predictor_names_valid, collapse = " + "))
      full_model <- stats::lm(as.formula(formula_str), data = lm_data)
      
      vif_results_raw <- car::vif(full_model)
      
      # Map names back to original if they were changed by make.names
      vif_names_map <- setNames(original_colnames[-1], valid_colnames[-1])
      names(vif_results_raw) <- vif_names_map[names(vif_results_raw)]
      
      vif_results <- vif_results_raw
      
      cat("  VIF values for the selected variables:\n")
      print(sort(vif_results, decreasing=TRUE))
      
      # Use helper to plot (pass the named vector)
      # Ensure your plotting helper uses the display name lookup from config
      plot_vif_results_original(vif_results, save_path = vif_plot_path, display_lookup = config$core_var_display_names) # Pass lookup
      
      # Identify high VIF variables
      high_vif_vars <- names(vif_results[vif_results > config$vif_threshold])
      if(length(high_vif_vars) > 0) {
        cat("  ---> Variables with VIF >", config$vif_threshold, ":", paste(high_vif_vars, collapse=", "), "\n")
      } else {
        cat("  All variables have VIF <= ", config$vif_threshold, ".\n")
      }
      
    }, error = function(e) {
      warning("VIF analysis failed for ", group_name, "-", scenario, ": ", e$message, "\n", call. = FALSE)
    })
    
    # --- Run Pearson Correlation (using corrplot) ---
    cat("\n  Running Pearson Correlation analysis...\n")
    tryCatch({
      # Use helper to plot (pass the dataframe)
      # Ensure your plotting helper uses the display name lookup from config
      plot_correlation_results_original(env_values_df, save_path = cor_plot_path, display_lookup = config$core_var_display_names) # Pass lookup
      
      # Identify highly correlated pairs
      cor_matrix <- cor(env_values_df, method = "pearson", use = "pairwise.complete.obs")
      cor_matrix_vals <- cor_matrix # Keep full matrix for display names later
      cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA # Keep only upper triangle for pair finding
      high_cor_pairs_indices <- which(abs(cor_matrix) > config$correlation_threshold, arr.ind = TRUE)
      
      if(nrow(high_cor_pairs_indices) > 0){
        cat("  ---> Highly correlated pairs (|r| >", config$correlation_threshold, "):\n")
        # Use display names for reporting
        display_names_in_matrix <- sapply(rownames(cor_matrix_vals), get_display_name, lookup = config$core_var_display_names)
        for(k in 1:nrow(high_cor_pairs_indices)) {
          row_idx <- high_cor_pairs_indices[k, 1]
          col_idx <- high_cor_pairs_indices[k, 2]
          var1_display <- display_names_in_matrix[row_idx]
          var2_display <- display_names_in_matrix[col_idx]
          cor_value <- cor_matrix_vals[row_idx, col_idx]
          cat("    - ", var1_display, " & ", var2_display, " (r = ", round(cor_value, 2), ")\n", sep="")
        }
      } else {
        cat("  No highly correlated pairs found (|r| <= ", config$correlation_threshold, ").\n")
      }
    }, error = function(e) {
      warning("Correlation analysis failed for ", group_name, "-", scenario, ": ", e$message, "\n", call. = FALSE)
    })
    
    cat("\n  Analysis complete for this set of variables.\n")
    cat("  Review outputs:\n")
    cat("    - VIF Plot:", vif_plot_path, "\n")
    cat("    - Correlation Plot:", cor_plot_path, "\n")
    cat("    - Console output above for high VIF/correlation.\n")
    cat("  Modify the 'variables_to_analyze' list at the top of script 05 and re-run if needed.\n")
    
    # Clean up memory for this group/scenario combo
    rm(group_occs_sf, env_values_df); gc()
    
  } # End loop groups
  rm(env_stack_processed); gc() # Clean up stack for the scenario
} # End loop scenarios

cat("\n--- Script 05 finished. Manually select final variables for SDM based on these outputs. ---\n")
#-------------------------------------------------------------------------------