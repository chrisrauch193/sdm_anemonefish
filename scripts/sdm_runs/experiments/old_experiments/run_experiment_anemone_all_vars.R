# scripts/sdm_runs/experiments/old_experiments/run_experiment_anemone_all_vars.R
#-------------------------------------------------------------------------------
# EXPERIMENTAL SCRIPT for Anemone SDMs
#
# Purpose: Run the VIF/PCA/SDM pipeline specifically for Anemone species.
#          Starts with a user-defined list of initial variables (using full names).
#          Runs VIF/correlation analysis on this INITIAL set.
#          Requires MANUAL EDITING to define the FINAL variable list for PCA/SDM.
#-------------------------------------------------------------------------------

# --- Clean Workspace (Optional) ---
# rm(list = ls())
# gc()

cat("--- Running EXPERIMENTAL Script: Anemone User-Defined Vars -> PCA -> SDM ---\n")

# --- 1. Load Configuration & Requirements ---
cat("\n--- Loading Configuration and Requirements ---\n")
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }
pacman::p_load(terra, sf, dplyr, readr, corrplot, ggplot2, car, vegan, ggvegan, tools, Hmisc, ENMeval) # Added required packages
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R")
if (!file.exists(env_helper_path)) stop("Env Helper missing.")
if (!file.exists(sdm_helper_path)) stop("SDM Helper missing.")
source(env_helper_path)
source(sdm_helper_path)
cat("--- Configuration and Helpers Loaded ---\n")


# --- 2. Define Experiment Specifics ---
exp_group_name <- "anemone"
exp_suffix <- "_exp_manualvars_fullname" # Suffix for output files/dirs
exp_log_dir <- file.path(config$log_dir, paste0("experiment_anemone", exp_suffix))
exp_pred_dir <- file.path(config$predictions_dir, paste0("experiment_anemone", exp_suffix))
exp_results_dir <- file.path(config$results_dir, paste0("experiment_anemone", exp_suffix))
exp_models_dir <- file.path(config$models_dir, paste0("experiment_anemone", exp_suffix))

# Create experimental output directories
dir.create(exp_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exp_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exp_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exp_models_dir, recursive = TRUE, showWarnings = FALSE)

cat("Experimental outputs will use suffix:", exp_suffix, "\n")
cat("Log/Plot directory:", exp_log_dir, "\n")
cat("Prediction directory:", exp_pred_dir, "\n")

# --- 3. Verify Necessary Inputs ---
cat("\n--- Verifying Inputs ---\n")
anemone_list_file <- config$anemone_species_list_file
if (!file.exists(anemone_list_file)) stop("Anemone list missing.")
anemone_species <- readr::read_csv(anemone_list_file, show_col_types = FALSE)
cat("  Anemone species list loaded.\n")
anemone_occ_dir <- config$anemone_occurrence_dir
if (!check_any_output_exists(anemone_occ_dir, "\\.csv$")) stop("Anemone occurrences missing.")
cat("  Anemone occurrence directory check passed.\n")
if (!dir.exists(config$scenario_folder_map$current) || !dir.exists(config$terrain_folder)) stop("Base env data missing.")
cat("  Base environmental directory check passed.\n")
cat("--- Input Verification Complete ---\n")


#-------------------------------------------------------------------------------
# !!! --- USER DEFINED: INITIAL FULL VARIABLE LIST (Using Cleaned Names) --- !!!
#-------------------------------------------------------------------------------
# Define the initial list using CLEANED layer names expected AFTER loading
# via the load_stack_env_data helper (which uses filename stems).
# *** VERIFY THESE NAMES AGAINST HELPER OUTPUT FOR YOUR DATA ***

initial_full_variable_list <- list()

# --- CURRENT Scenario Variables ---
initial_full_variable_list[['current']] <- c(
  "no3_baseline_depthmax_ltmax",
  "no3_baseline_depthmax_ltmin",
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  "o2_baseline_depthmax_ltmax",
  "o2_baseline_depthmax_ltmin",
  "o2_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  "phyc_baseline_depthmax_mean",
  "so_baseline_depthmax_mean",
  "sws_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_ltmax",
  "thetao_baseline_depthmax_ltmin",
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_range",
  "thetao_baseline_depthsurf_mean",
  "chl_baseline_depthmax_mean",
  "par_baseline_depthsurf_mean",
  # Terrain Vars
  "bathymetry_mean", "slope", "rugosity", "distcoast"
)

# # --- SSP119 dec50 (2050) ---
# # List based on your debug output for ssp119 dec50 files + terrain
# initial_full_variable_list[['ssp119_2050']] <- c(
#   "no3_ssp119_depthmax_dec50_ltmax",
#   "no3_ssp119_depthmax_dec50_ltmin",
#   "no3_ssp119_depthmax_dec50_mean",
#   "no3_ssp119_depthmax_dec50_range",
#   "o2_ssp119_depthmax_dec50_ltmax",
#   "o2_ssp119_depthmax_dec50_ltmin",
#   "o2_ssp119_depthmax_dec50_mean",
#   "o2_ssp119_depthmax_dec50_range",
#   "phyc_ssp119_depthmax_dec50_mean",
#   "so_ssp119_depthmax_dec50_mean",
#   "sws_ssp119_depthsurf_dec50_mean", # Assuming depthsurf
#   "thetao_ssp119_depthmax_dec50_ltmax",
#   "thetao_ssp119_depthmax_dec50_ltmin",
#   "thetao_ssp119_depthmax_dec50_mean",
#   "thetao_ssp119_depthmax_dec50_range",
#   "thetao_ssp119_depthsurf_dec50_mean", # Assuming depthsurf
#   # Copied & Renamed Baseline Vars
#   "chl_ssp119_depthmax_dec50_mean",
#   "par_ssp119_depthsurf_dec50_mean", # Assuming depthsurf
#   # Terrain Vars
#   "bathymetry_mean", "slope", "rugosity", "distcoast"
# )


# --- FIRST CUT Scenario Variables ---
initial_full_variable_list[['current']] <- c(
  "par_baseline_depthsurf_mean",
  "sws_baseline_depthsurf_mean",
  # "thetao_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_range",
  "so_baseline_depthmax_mean",
  "no3_baseline_depthmax_mean",
  "no3_baseline_depthmax_range",
  # "no3_baseline_depthmax_ltmax",
  # "no3_baseline_depthmax_ltmin",
  "chl_baseline_depthmax_mean",
  # "o2_baseline_depthmax_mean",
  "o2_baseline_depthmax_range",
  # "o2_baseline_depthmax_ltmax",
  # "o2_baseline_depthmax_ltmin",
  # "phyc_baseline_depthmax_mean",
  # "ph_baseline_depthmax_mean",
  # Terrain Vars
  "bathymetry_mean", "distcoast", "rugosity"
)

# PROPOSED ANEMONEFISH Starting List
initial_full_variable_list[['current']] <- c(
  "par_baseline_depthsurf_mean",    # Vision/orientation
  "sws_baseline_depthsurf_mean",    # Lower priority - larval transport?
  "thetao_baseline_depthmax_mean",  # Juvenile/Adult temp
  "thetao_baseline_depthmax_range", # Temp variability/stress
  "so_baseline_depthmax_mean",      # Osmoregulation
  "no3_baseline_depthmax_mean",     # Larval food web proxy
  "no3_baseline_depthmax_range",  # Consider removing (lower priority for fish)
  "chl_baseline_depthmax_mean",     # Larval food web proxy
  "o2_baseline_depthmax_ltmin",     # ADDED/SWAPPED - Physiological limit
  # Terrain Vars
  "bathymetry_mean",                # Depth limit
  "distcoast",                      # Coastal affinity
  "rugosity"                        # Habitat structure
)

# --- SSP119 dec50 (2050) ---
# List based on your debug output for ssp119 dec50 files + terrain
initial_full_variable_list[['ssp119_2050']] <- c(
  "par_ssp119_depthsurf_dec50_mean",
  "sws_ssp119_depthsurf_dec50_mean",
  "thetao_ssp119_depthmax_dec50_mean",
  "thetao_ssp119_depthmax_dec50_range",
  "so_ssp119_depthmax_dec50_mean",
  "no3_ssp119_depthmax_dec50_mean",
  "no3_ssp119_depthmax_dec50_range",
  "chl_ssp119_depthmax_dec50_mean",
  "o2_ssp119_depthmax_dec50_range",
  # Terrain Vars
  "bathymetry_mean", "distcoast", "rugosity"
)

# --- SSP119 dec100 (2100) ---
# Construct names based on the dec50 list, just changing the tag
initial_full_variable_list[['ssp119_2100']] <- gsub("_dec50_", "_dec100_", initial_full_variable_list[['ssp119_2050']])
initial_full_variable_list[['ssp119_2100']] <- unique(initial_full_variable_list[['ssp119_2100']]) # Ensure unique

# --- SSP585 dec50 (2050) ---
# Construct names based on the ssp119_2050 list, changing the SSP
initial_full_variable_list[['ssp585_2050']] <- gsub("ssp119", "ssp585", initial_full_variable_list[['ssp119_2050']])
initial_full_variable_list[['ssp585_2050']] <- unique(initial_full_variable_list[['ssp585_2050']]) # Ensure unique

# --- SSP585 dec100 (2100) ---
# Construct names based on the ssp119_2100 list, changing the SSP
initial_full_variable_list[['ssp585_2100']] <- gsub("ssp119", "ssp585", initial_full_variable_list[['ssp119_2100']])
initial_full_variable_list[['ssp585_2100']] <- unique(initial_full_variable_list[['ssp585_2100']]) # Ensure unique


# Verify constructed lists (optional print)
# print(initial_full_variable_list[['ssp119_2100']])
# print(initial_full_variable_list[['ssp585_2050']])
# print(initial_full_variable_list[['ssp585_2100']])

cat("--- Initial variable lists defined for experiment ---\n")
#-------------------------------------------------------------------------------

# --- (Rest of the script: VIF/Corr loop, Manual selection, PCA/SDM loop remains the same) ---

# --- 4. VIF/Correlation Analysis Loop (No Step, Uses Initial List) ---
cat("\n--- PART 1: VIF & Correlation Analysis (Anemones, User Initial List) ---\n")
vif_results_full <- list(); correlation_plots <- list()
run_vif_correlation <- TRUE
if (run_vif_correlation) {
  for (scenario in names(initial_full_variable_list)) {
    cat("\n-----------------\n Analyzing Scenario:", scenario, "\n-----------------\n")
    current_initial_vars <- initial_full_variable_list[[scenario]]
    if(is.null(current_initial_vars) || length(current_initial_vars) < 2) { warning("Skipping: Initial list too short."); next }
    cat("  Using initial variables:", paste(current_initial_vars, collapse=", "), "\n")
    env_stack_raw <- load_stack_env_data(scenario, config); if (is.null(env_stack_raw)) {warning("Failed load env"); next}
    # Ensure ALL initial vars are actually present in the loaded stack
    if (!all(current_initial_vars %in% names(env_stack_raw))) {
      missing_vars <- current_initial_vars[!current_initial_vars %in% names(env_stack_raw)]
      warning("Cannot find variables: '", paste(missing_vars, collapse=","), "' in loaded stack for scenario ", scenario, ". Skipping.", call.=FALSE)
      print("Available layers:"); print(names(env_stack_raw))
      rm(env_stack_raw); gc(); next
    }
    env_stack_subset <- env_stack_raw[[current_initial_vars]]; env_stack_processed <- preprocess_env_rasters(env_stack_subset, config); if (is.null(env_stack_processed)) {warning("Failed preprocess env"); next}
    rm(env_stack_raw, env_stack_subset); gc(); target_crs <- terra::crs(env_stack_processed)
    group_occs_sf <- load_clean_thin_group_occurrences(anemone_occ_dir, target_crs, config, env_stack_processed); if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) {warning("No occurrences"); next}
    env_values_df <- extract_env_values(group_occs_sf, env_stack_processed); if (is.null(env_values_df) || ncol(env_values_df) < 2) {warning("Insufficient extracted data"); next}
    variances <- apply(env_values_df, 2, var, na.rm = TRUE); zero_var_cols <- names(variances[variances == 0]); if (length(zero_var_cols) > 0) { cat("  Warning: Removing zero-var cols:", paste(zero_var_cols, collapse=", "), "\n"); env_values_df <- env_values_df[, !(names(env_values_df) %in% zero_var_cols), drop = FALSE] }; if (ncol(env_values_df) < 2) {warning("Less than 2 vars remaining"); next}
    
    # VIF Analysis (on initial vars)
    cat("  Running VIF analysis on initial set...\n"); vif_on_initial <- NULL
    tryCatch({
      present_dummy <- rep(1, nrow(env_values_df)); lm_data <- cbind(present_dummy, env_values_df)
      valid_names <- make.names(names(env_values_df), unique=TRUE); colnames(lm_data) <- c("present_dummy", valid_names)
      formula_str <- paste("present_dummy ~", paste(valid_names, collapse = " + "))
      full_model_initial <- stats::lm(as.formula(formula_str), data = lm_data)
      vif_on_initial <- car::vif(full_model_initial); names(vif_on_initial) <- names(env_values_df) # Map back names
      cat("  VIF values for initial variable set:\n"); print(sort(vif_on_initial, decreasing=TRUE))
      vif_results_full[[scenario]] <- vif_on_initial
      vif_plot_path <- file.path(exp_log_dir, paste0("vif_plot_initial_", exp_group_name, "_", scenario, exp_suffix, ".png"))
      plot_vif_results_original(vif_on_initial, save_path = vif_plot_path)
    }, error = function(e) { warning("VIF analysis failed: ", e$message) })
    
    # Pearson Correlation (on initial vars)
    cat("  Running Correlation Analysis on initial set...\n")
    cor_plot_path <- file.path(exp_log_dir, paste0("correlation_plot_initial_", exp_group_name, "_", scenario, exp_suffix, ".png"))
    plot_correlation_results_original(env_values_df, save_path = cor_plot_path)
    correlation_plots[[scenario]] <- cor_plot_path
    rm(env_stack_processed, group_occs_sf, env_values_df); gc()
  } # End scenario loop
} else { cat ("Skipping VIF/Correlation analysis.\n") }
cat("\n--- PART 1 FINISHED: Review INITIAL VIF/Correlation plots in", exp_log_dir, "---\n")

#-------------------------------------------------------------------------------
# !!! --- MANUAL STEP: DEFINE FINAL VARIABLES FOR PCA --- !!!
#-------------------------------------------------------------------------------
# Based on INITIAL VIF/Corr plots, define the FINAL list for PCA below.
# Use the EXACT layer names from the 'Initial variables for VIF/Corr' printout above.

# stop("SCRIPT HALTED: Review INITIAL VIF/Corr plots in the log directory then define 'final_selected_vars_experiment' list below and comment out this stop() line.")

final_selected_vars_experiment <- list()

# --- Fill this list MANUALLY using EXACT layer names ---
# Example (replace with your actual choices after review):
final_selected_vars_experiment[['current']] <- c(
  "chl_baseline_2000_2018_depthmax_mean",
  "no3_baseline_2000_2018_depthmax_ltmin",
  "o2_baseline_2000_2018_depthmax_range",
  "par_mean_baseline_2000_2020_depthsurf",
  "so_baseline_2000_2019_depthmax_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)

final_selected_vars_experiment[['ssp119_2050']] <- c(
  "chl_baseline_2000_2018_depthmax_mean", # Copied from current
  "no3_ssp119_depthmax_dec50_ltmin",
  "o2_ssp119_depthmax_dec50_range",
  "par_mean_baseline_2000_2020_depthsurf", # Copied from current
  "so_ssp119_depthmax_dec50_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)

final_selected_vars_experiment[['ssp119_2100']] <- c(
  "chl_baseline_2000_2018_depthmax_mean", # Copied from current
  "no3_ssp119_depthmax_dec100_ltmin",
  "o2_ssp119_depthmax_dec100_range",
  "par_mean_baseline_2000_2020_depthsurf", # Copied from current
  "so_ssp119_depthmax_dec100_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)

final_selected_vars_experiment[['ssp585_2050']] <- c(
  "chl_baseline_2000_2018_depthmax_mean", # Copied from current
  "no3_ssp585_depthmax_dec50_ltmin",
  "o2_ssp585_depthmax_dec50_range",
  "par_mean_baseline_2000_2020_depthsurf", # Copied from current
  "so_ssp585_depthmax_dec50_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)

final_selected_vars_experiment[['ssp585_2100']] <- c(
  "chl_baseline_2000_2018_depthmax_mean", # Copied from current
  "no3_ssp585_depthmax_dec100_ltmin",
  "o2_ssp585_depthmax_dec100_range",
  "par_mean_baseline_2000_2020_depthsurf", # Copied from current
  "so_ssp585_depthmax_dec100_mean",
  "bathymetry_mean",
  "distcoast",
  "rugosity"
)
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# --- Safety Check ---
needed_scenarios <- config$env_scenarios
missing_selections <- setdiff(needed_scenarios, names(final_selected_vars_experiment))
if(length(missing_selections) > 0) {
  stop("FATAL: Manual variable selections missing for: ", paste(missing_selections, collapse=", "))
}
cat("\n--- MANUAL STEP COMPLETE: Using final variable selections for PCA/SDM ---\n")
#-------------------------------------------------------------------------------


# --- 5. PCA & Projection Loop ---
cat("\n--- PART 2: PCA Analysis & Projection (Anemones, Final Selected Vars) ---\n")
run_pca_projection <- TRUE
if(run_pca_projection) {
  # Loop using names from the final selection list
  for (scenario in names(final_selected_vars_experiment)) {
    cat("\n-----------------------------------------------------\n")
    cat(" Running PCA & Projection for Scenario:", scenario, "\n")
    cat("-----------------------------------------------------\n")
    final_selected_vars <- final_selected_vars_experiment[[scenario]]
    if(is.null(final_selected_vars) || length(final_selected_vars) < 2) { warning("Skipping PCA: < 2 vars selected."); next }
    cat("  Using variables for PCA:", paste(final_selected_vars, collapse=", "), "\n")
    # Reload data - ensure it contains the selected variables
    env_stack_raw <- load_stack_env_data(scenario, config); if (is.null(env_stack_raw)) {warning("Failed load env"); next}
    # Check if ALL selected vars are in the loaded stack BEFORE subsetting/processing
    if (!all(final_selected_vars %in% names(env_stack_raw))) {
      missing_vars <- final_selected_vars[!final_selected_vars %in% names(env_stack_raw)]
      warning("FATAL: Manually selected variables '", paste(missing_vars, collapse=","), "' not found in loaded stack for scenario ", scenario, ". Check initial list and helper output names. Skipping.", call.=FALSE)
      print("Available layers:"); print(names(env_stack_raw))
      rm(env_stack_raw); gc(); next
    }
    env_stack_subset <- env_stack_raw[[final_selected_vars]] # Subset to ONLY selected vars BEFORE processing
    env_stack_processed <- preprocess_env_rasters(env_stack_subset, config); if (is.null(env_stack_processed)) {warning("Failed preprocess env"); next}
    rm(env_stack_raw, env_stack_subset); gc(); target_crs <- terra::crs(env_stack_processed)
    group_occs_sf <- load_clean_thin_group_occurrences(anemone_occ_dir, target_crs, config, env_stack_processed); if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) {warning("No occurrences"); next}
    env_values_df <- extract_env_values(group_occs_sf, env_stack_processed); if (is.null(env_values_df) || ncol(env_values_df) < 2) {warning("Insufficient extracted data"); next}
    # Use the already selected df (env_values_df now only contains selected vars)
    env_values_df_selected <- env_values_df
    if (nrow(env_values_df_selected) < 2) {warning("Insufficient rows after selection/NA removal for PCA"); next}
    
    # Run PCA (Original Method)
    cat("  Running PCA (vegan::rda)...\n"); pca_model_obj <- NULL
    tryCatch({
      env_std <- scale(env_values_df_selected, center = TRUE, scale = TRUE); pca_model_obj <- vegan::rda(env_std)
      # Plotting
      pca_biplot_path <- file.path(exp_log_dir, paste0("pca_biplot_", exp_group_name, "_", scenario, exp_suffix, ".png"))
      pca_contribution_path <- file.path(exp_log_dir, paste0("pca_contribution_", exp_group_name, "_", scenario, exp_suffix, ".png"))
      # (Plotting code remains the same...)
      pca_summary <- summary(pca_model_obj); prop_explained <- pca_summary$cont$importance[2, 1:2]; site_scores <- as.data.frame(vegan::scores(pca_model_obj, display = "sites", choices = 1:2)); species_scores <- as.data.frame(vegan::scores(pca_model_obj, display = "species", choices = 1:2)); scaling_factor <- min((max(site_scores$PC1)-min(site_scores$PC1))/(max(species_scores$PC1)-min(species_scores$PC1)), (max(site_scores$PC2)-min(site_scores$PC2))/(max(species_scores$PC2)-min(species_scores$PC2)), na.rm=TRUE) * 0.9; species_scores_scaled <- species_scores * scaling_factor
      pca_plot <- ggplot() + geom_point(data = site_scores, aes(x = PC1, y = PC2), color = "salmon", size = 1) + geom_segment(data = species_scores_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.25, "cm")), color = "black") + geom_text(data = species_scores_scaled, aes(x = PC1 * 1.1, y = PC2 * 1.1, label = rownames(species_scores_scaled)), hjust = 0.5, vjust = 0.5, size = 3, color = "black", check_overlap = TRUE) + labs(x = paste("PCA Axis 1 (", round(prop_explained[1] * 100, 2), "%)", sep = ""), y = paste("PCA Axis 2 (", round(prop_explained[2] * 100, 2), "%)", sep = ""), title = paste("PCA Biplot (Exp) -", exp_group_name, scenario)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      ggsave(filename = pca_biplot_path, plot = pca_plot, width = 8, height = 7)
      var_loadings <- vegan::scores(pca_model_obj, display = "species", choices = 1:config$n_pca_components); var_loadings_df <- as.data.frame(var_loadings); var_loadings_df$Variable <- rownames(var_loadings_df); var_loadings_long <- tidyr::pivot_longer(var_loadings_df, cols = starts_with("PC"), names_to = "Principal_Component", values_to = "Contribution"); var_loadings_long$Principal_Component <- factor(var_loadings_long$Principal_Component, levels = paste0("PC", 1:config$n_pca_components)); num_vars <- length(unique(var_loadings_long$Variable)); colors <- grDevices::colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars)
      contribution_plot <- ggplot(var_loadings_long, aes(x = Principal_Component, y = Contribution, fill = Variable)) + geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = colors) + labs(title = paste("Variable Contribution (Exp) -", exp_group_name, scenario), x = "Principal Component", y = "Contribution") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
      ggsave(filename = pca_contribution_path, plot = contribution_plot, width = 10, height = 6)
      cat("  PCA plots saved.\n"); experimental_pca_models[[scenario]] <- pca_model_obj
    }, error = function(e) { warning("PCA failed: ", e$message); pca_model_obj <- NULL })
    
    # Project PCA onto Rasters
    if (!is.null(pca_model_obj)) {
      cat("  Projecting PCA onto rasters...\n")
      tryCatch({
        env_stack_selected_for_pca <- env_stack_processed # Already subsetted
        pca_prcomp <- stats::prcomp(env_values_df_selected, center = TRUE, scale. = TRUE)
        pca_rasters <- terra::predict(env_stack_selected_for_pca, pca_prcomp, index = 1:config$n_pca_components); names(pca_rasters) <- paste0("PC", 1:config$n_pca_components)
        pca_raster_file <- file.path(exp_log_dir, paste0("pca_rasters_", exp_group_name, "_", scenario, exp_suffix, ".tif"))
        terra::writeRaster(pca_rasters, filename = pca_raster_file, overwrite = TRUE, gdal=c("COMPRESS=LZW"))
        cat("  Experimental PCA raster saved:", pca_raster_file, "\n"); experimental_pca_paths[[scenario]] <- pca_raster_file
      }, error = function(e) { warning("PCA projection failed: ", e$message) })
    }
    rm(env_stack_processed, group_occs_sf, env_values_df, env_values_df_selected, pca_model_obj); gc()
  } # End scenario loop PCA
} else { cat("Skipping PCA/Projection.\n"); temp_pca_paths_path <- file.path(exp_log_dir, "experimental_pca_paths.rds"); if(file.exists(temp_pca_paths_path)) {experimental_pca_paths <- readRDS(temp_pca_paths_path); cat("Loaded previous exp PCA paths.\n")} else {warning("Cannot load exp PCA paths.")} }
saveRDS(experimental_pca_paths, file.path(exp_log_dir, "experimental_pca_paths.rds"))
saveRDS(experimental_pca_models, file.path(exp_log_dir, "experimental_pca_models.rds"))
cat("\n--- PART 2 FINISHED ---\n")


#-------------------------------------------------------------------------------
# --- 6. Run Standard SDMs Loop (Anemones Only, Using Experimental PCA Layers) ---
#-------------------------------------------------------------------------------
cat("\n--- PART 3: Running Standard SDMs (Anemones Only, Experimental Predictors) ---\n")
if(length(experimental_pca_paths) == 0) stop("Exp PCA paths missing.")
run_sdms <- TRUE
if(run_sdms) {
  # Loop using names from the generated PCA paths list
  needed_scenarios_for_sdm <- names(experimental_pca_paths)
  if(length(needed_scenarios_for_sdm) == 0) stop("No experimental PCA paths found to run SDMs.")
  
  total_sdms_run <- 0; total_sdms_skipped <- 0
  for (i in 1:nrow(anemone_species)) {
    species_name_sanitized <- gsub(" ", "_", anemone_species$scientificName[i]); species_aphia_id <- anemone_species$AphiaID[i]
    cat("\n------------------\nProcessing Species:", anemone_species$scientificName[i], "\n------------------\n")
    occ_sf_clean <- load_clean_individual_occ(species_aphia_id, anemone_occ_dir, config)
    if(is.null(occ_sf_clean) || nrow(occ_sf_clean) < config$min_occurrences_sdm) { warning("Skipping: Not enough occurrences."); total_sdms_skipped <- total_sdms_skipped + length(needed_scenarios_for_sdm); next }
    for (scenario in needed_scenarios_for_sdm) {
      cat("  --- Scenario:", scenario, "---\n")
      pred_file <- file.path(exp_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", scenario, exp_suffix, ".tif"))
      results_file <- file.path(exp_results_dir, paste0("sdm_results_", species_name_sanitized, "_", scenario, exp_suffix, ".rds"))
      eval_file <- file.path(exp_results_dir, paste0("sdm_eval_", species_name_sanitized, "_", scenario, exp_suffix, ".csv"))
      model_obj_file <- file.path(exp_models_dir, paste0("sdm_model_", species_name_sanitized, "_", scenario, exp_suffix, ".rds"))
      pca_stack_path <- experimental_pca_paths[[scenario]] # Path comes from the list now
      if (is.null(pca_stack_path) || !file.exists(pca_stack_path)) { warning("Exp PCA stack missing."); total_sdms_skipped <- total_sdms_skipped + 1; next }
      if (file.exists(pred_file)) { cat("    Prediction file exists. Skipping.\n"); total_sdms_skipped <- total_sdms_skipped + 1; next }
      predictor_stack <- tryCatch(terra::rast(pca_stack_path), error = function(e) NULL)
      if(is.null(predictor_stack)) { warning("Failed load PCA stack."); total_sdms_skipped <- total_sdms_skipped + 1; next }
      occ_sf_thinned <- thin_individual_occ(occ_sf_clean, predictor_stack, config)
      if(is.null(occ_sf_thinned) || nrow(occ_sf_thinned) < config$min_occurrences_sdm) { warning("Skipping: Not enough occs after thinning."); total_sdms_skipped <- total_sdms_skipped + 1; next }
      background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config); if (is.null(background_points)) {warning("Failed background gen."); total_sdms_skipped <- total_sdms_skipped + 1; next}
      enmeval_results <- run_sdm_enmeval(occ_sf_thinned, predictor_stack, background_points, config); if (is.null(enmeval_results)) {warning("ENMeval failed."); total_sdms_skipped <- total_sdms_skipped + 1; next}
      tryCatch(saveRDS(enmeval_results, file = results_file), error=function(e){warning("Failed save ENMeval results.")}); tryCatch(readr::write_csv(enmeval_results@results, file = eval_file), error=function(e){warning("Failed save Eval table.")})
      prediction_raster <- predict_sdm_best(enmeval_results, predictor_stack, config); if (is.null(prediction_raster)) {warning("Prediction failed."); total_sdms_skipped <- total_sdms_skipped + 1; next}
      tryCatch({ terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES")); cat("    Prediction raster saved.\n"); total_sdms_run <- total_sdms_run + 1 }, error=function(e){warning("Failed save prediction raster."); total_sdms_skipped <- total_sdms_skipped + 1})
      tryCatch({eval_results_df<-enmeval_results@results; metric<-config$sdm_evaluation_metric; if(!metric %in% names(eval_results_df)) metric<-"AICc"; select_best<-if(metric=="AICc") which.min else which.max; valid_rows<-!is.na(eval_results_df[[metric]]); if(sum(valid_rows)>0){best_model_index<-select_best(eval_results_df[[metric]][valid_rows]); original_indices<-which(valid_rows); best_model_row_index<-original_indices[best_model_index]; best_model_obj<-enmeval_results@models[[best_model_row_index]]; saveRDS(best_model_obj, file=model_obj_file); cat("    Best model object saved.\n")}}, error=function(e){warning("Failed to save best model object: ", e$message)})
      rm(predictor_stack, occ_sf_thinned, background_points, enmeval_results, prediction_raster); gc()
    } # End scenario loop
    rm(occ_sf_clean); gc()
  } # End species loop
  cat("\n--- Exp Anemone SDM runs finished ---\n"); cat("Total SDMs run:", total_sdms_run, "\n"); cat("Total SDMs skipped:", total_sdms_skipped, "\n")
} else { cat("Skipping SDM runs.\n") }

cat("\n--- EXPERIMENTAL Script Finished ---\n")
#-------------------------------------------------------------------------------






# Start final 8
# mean surface temperature
# mean current velocity (surface wind speed)
# mean salinity
# mean temperature
# mean nitrate concentration
# nitrate concentration range
# mean chlorophyll concentration
# dissolved oxygen concentration range
# mean phytoplankton concentration.
# bath >= 50m

# General start
# mean depthsurf light availability (PAR)
# mean depthsurf wind speed
# mean depthsurf temperature
# mean depthmax temperature
# ltmin depthmax temperature
# ltmax depthmax temperature
# mean depthmax salinity
# mean depthmax nitrate concentration
# ltmin depthmax nitrate concentration
# ltmax depthmax nitrate concentration
# mean depthmax chlorophyll concentration
# mean depthmax dissolved oxygen concentration
# ltmin depthmax dissolved oxygen concentration
# ltmax depthmax dissolved oxygen concentration
# mean depthmax phytoplankton concentration.
# distcoast
# rugosity
# bath >= 50m


# Initial Selection (TODO: Check colinearity)
# mean depthmax light availability (PAR)
# mean surface wind speed
# mean surface temperature
# mean depthmax temperature
# ltmin depthmax temperature
# ltmax depthmax temperature
# mean depthmax salinity
# ltmin depthmax nitrate concentration
# ltmax depthmax nitrate concentration
# mean depthmax chlorophyll concentration
# ltmin depthmax dissolved oxygen concentration
# ltmax depthmax dissolved oxygen concentration
# mean depthmax phytoplankton concentration.
# distcoast
# rugosity
# bath >= 50m

# Potential additions for anemonefish only
# current speed at different depths (for dispersal, maybe add in dispersal land barriers)






#FINAL FINAL ANEMONE HOST
initial_full_variable_list[['current']] <- c(
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

# PROPOSED ANEMONEFISH Starting List
initial_full_variable_list[['current']] <- c(
  "par_baseline_depthsurf_mean",    # Vision/orientation
  "sws_baseline_depthsurf_mean",    # Lower priority - larval transport?
  "thetao_baseline_depthsurf_mean", # ADDED - Larval temp
  "thetao_baseline_depthmax_mean",  # Juvenile/Adult temp
  "thetao_baseline_depthmax_range", # Temp variability/stress
  "so_baseline_depthmax_mean",      # Osmoregulation
  # "no3_baseline_depthmax_mean",     # Larval food web proxy
  # "no3_baseline_depthmax_range",  # Consider removing (lower priority for fish)
  "chl_baseline_depthmax_mean",     # Larval food web proxy
  "o2_baseline_depthmax_ltmin",     # ADDED/SWAPPED - Physiological limit
  # "ph_baseline_depthmax_mean",           # ADDED - CRITICAL for behavior/physiology
  # Terrain Vars
  "bathymetry_mean",                # Depth limit
  "distcoast",                      # Coastal affinity
  "rugosity"                        # Habitat structure
)

