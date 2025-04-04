# scripts/run_experiment_anemone_all_vars.R
#-------------------------------------------------------------------------------
# EXPERIMENTAL SCRIPT for Anemone SDMs
#
# Purpose: Run the VIF/PCA/SDM pipeline specifically for Anemone species,
#          starting VIF/correlation analysis with ALL available environmental
#          variables for each scenario. Allows manual selection of variables
#          for PCA within this script for experimentation.
#
# !!! IMPORTANT !!!
# 1. Assumes scripts 01-04 have been run successfully via 00_run_all_scripts.R
#    (i.e., packages installed, species lists exist, occurrences downloaded,
#     environmental data downloaded and processed).
# 2. Outputs from this script (plots, PCA rasters, SDM results) will be saved
#    with an "_exp_allvars" suffix or in specific experimental subdirectories
#    to avoid overwriting main workflow results.
# 3. Requires MANUAL EDITING of this script after the first run section
#    to define the 'final_selected_vars_experiment' list before running
#    the PCA and SDM sections.
#-------------------------------------------------------------------------------

# --- Clean Workspace (Optional) ---
# rm(list = ls())
# gc()

cat("--- Running EXPERIMENTAL Script: Anemone Full Env Vars -> Manual Select -> PCA -> SDM ---\n")

# --- 1. Load Configuration & Requirements ---
cat("\n--- Loading Configuration and Requirements ---\n")
if (!exists("config")) {
  source("scripts/config.R")
  if (!exists("config")) stop("Failed to load config object from scripts/config.R")
}
# Ensure required packages are loaded
pacman::p_load(terra, sf, dplyr, readr, corrplot, ggplot2, car, vegan, ggvegan, tools, Hmisc, ENMeval)

# Source helper functions
env_helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R")
if (!file.exists(env_helper_path)) stop("Env Helper file not found: ", env_helper_path)
if (!file.exists(sdm_helper_path)) stop("SDM Helper file not found: ", sdm_helper_path)
source(env_helper_path)
source(sdm_helper_path)
cat("--- Configuration and Helpers Loaded ---\n")


# --- 2. Define Experiment Specifics ---
exp_group_name <- "anemone"
exp_suffix <- "_exp_allvars" # Suffix for output files/dirs
exp_log_dir <- file.path(config$log_dir, paste0("experiment_anemone", exp_suffix))
exp_pred_dir <- file.path(config$predictions_dir, paste0("experiment_anemone", exp_suffix))
exp_results_dir <- file.path(config$results_dir, paste0("experiment_anemone", exp_suffix))
exp_models_dir <- file.path(config$models_dir, paste0("experiment_anemone", exp_suffix))

# Create experimental output directories
dir.create(exp_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exp_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exp_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exp_models_dir, recursive = TRUE, showWarnings = FALSE)

cat("Experimental outputs will be saved with suffix:", exp_suffix, "\n")
cat("Log/Plot directory:", exp_log_dir, "\n")
cat("Prediction directory:", exp_pred_dir, "\n")

# --- 3. Verify Necessary Inputs from Main Workflow ---
cat("\n--- Verifying Inputs ---\n")
# Check species list
anemone_list_file <- config$anemone_species_list_file
if (!file.exists(anemone_list_file)) stop("Anemone species list not found: ", anemone_list_file)
anemone_species <- readr::read_csv(anemone_list_file, show_col_types = FALSE)
cat("  Anemone species list loaded.\n")

# Check occurrence directory
anemone_occ_dir <- config$anemone_occurrence_dir
if (!check_any_output_exists(anemone_occ_dir, "\\.csv$")) {
  stop("Anemone occurrence files seem missing from: ", anemone_occ_dir)
}
cat("  Anemone occurrence directory check passed.\n")

# Check environmental data directories (basic check)
if (!dir.exists(config$scenario_folder_map$current) || !dir.exists(config$terrain_folder)) {
  stop("Base environmental data directories (current/terrain) not found in: ", config$env_data_dir)
}
cat("  Base environmental directory check passed.\n")
cat("--- Input Verification Complete ---\n")


# --- 4. VIF/Correlation Analysis Loop (Run this part first) ---
cat("\n--- PART 1: VIF & Correlation Analysis (Anemones, All Initial Vars) ---\n")
# This stores the variables suggested by VIF/step before manual filtering
suggested_vars_from_step <- list()
# This will store paths to experimental PCA rasters eventually
experimental_pca_paths <- list()
experimental_pca_models <- list() # Store vegan::rda models

run_vif_correlation <- TRUE # Set to FALSE if you have already run this and just want to define variables and run PCA/SDM

if (run_vif_correlation) {
  for (scenario in config$env_scenarios) {
    cat("\n-----------------------------------------------------\n")
    cat(" Analyzing Scenario:", scenario, "for VIF/Correlation\n")
    cat("-----------------------------------------------------\n")
    
    # Load/Preprocess Env Data (using helpers)
    env_stack_raw <- load_stack_env_data(scenario, config)
    if (is.null(env_stack_raw)) {warning("Failed load env for ", scenario); next}
    env_stack_processed <- preprocess_env_rasters(env_stack_raw, config)
    if (is.null(env_stack_processed)) {warning("Failed preprocess env for ", scenario); next}
    rm(env_stack_raw); gc()
    target_crs <- terra::crs(env_stack_processed)
    
    # Load/Clean/Thin Group Occurrences (Anemones)
    group_occs_sf <- load_clean_thin_group_occurrences(anemone_occ_dir, target_crs, config, env_stack_processed)
    if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) {warning("No valid occurrences for ", scenario); next}
    
    # Extract Env Values (using ALL layers initially)
    env_values_df <- extract_env_values(group_occs_sf, env_stack_processed)
    if (is.null(env_values_df) || ncol(env_values_df) < 2) {warning("Insufficient env data extracted for ", scenario); next}
    
    # Check/Remove Zero Variance
    variances <- apply(env_values_df, 2, var, na.rm = TRUE)
    zero_var_cols <- names(variances[variances == 0])
    if (length(zero_var_cols) > 0) {
      cat("  Warning: Removing zero-variance columns:", paste(zero_var_cols, collapse=", "), "\n")
      env_values_df <- env_values_df[, !(names(env_values_df) %in% zero_var_cols), drop = FALSE]
    }
    if (ncol(env_values_df) < 2) {warning("Less than 2 vars remaining for ", scenario); next}
    
    cat("  Initial variables for VIF/Corr:", paste(names(env_values_df), collapse=", "), "\n")
    
    # Run VIF/step (Original Method)
    cat("  Running VIF/step...\n")
    step_selected_vars <- NULL
    tryCatch({
      present_dummy <- rep(1, nrow(env_values_df))
      lm_data <- cbind(present_dummy, env_values_df)
      full_model <- stats::lm(present_dummy ~ ., data = lm_data)
      stepped_model <- stats::step(full_model, direction = "both", trace = 0)
      step_selected_vars_lm <- names(coefficients(stepped_model))[-1]
      step_selected_vars_lm <- gsub("^`|`$", "", step_selected_vars_lm) # Clean names
      
      if(length(step_selected_vars_lm) >= 2) {
        vif_on_stepped <- car::vif(stepped_model)
        # Save VIF plot with experimental suffix
        vif_plot_path <- file.path(exp_log_dir, paste0("vif_plot_", exp_group_name, "_", scenario, exp_suffix, ".png"))
        plot_vif_results_original(vif_on_stepped, save_path = vif_plot_path)
        step_selected_vars <- step_selected_vars_lm
      } else {
        warning("step() selected less than 2 variables for ", scenario)
        step_selected_vars <- step_selected_vars_lm # Still store the limited result
      }
    }, error = function(e) {
      warning("VIF/step failed for ", scenario, ": ", e$message, ". Using all variables.")
      step_selected_vars <- names(env_values_df) # Fallback
    })
    
    # Run Correlation Plot (Original Method)
    if(!is.null(step_selected_vars) && length(step_selected_vars) >= 2) {
      cat("  Running Correlation Analysis on step-selected vars...\n")
      cor_plot_path <- file.path(exp_log_dir, paste0("correlation_plot_", exp_group_name, "_", scenario, exp_suffix, ".png"))
      plot_correlation_results_original(env_values_df[, step_selected_vars, drop=FALSE], save_path = cor_plot_path)
      # Save suggested vars list
      suggested_vars_path <- file.path(exp_log_dir, paste0("suggested_vars_after_step_", exp_group_name, "_", scenario, exp_suffix, ".txt"))
      writeLines(step_selected_vars, suggested_vars_path)
      suggested_vars_from_step[[scenario]] <- step_selected_vars
      cat("  Suggested variables and plots saved to:", exp_log_dir, "\n")
    } else {
      cat("  Skipping correlation plot as < 2 variables selected by step.\n")
      suggested_vars_from_step[[scenario]] <- step_selected_vars # Store limited/null result
    }
    rm(env_stack_processed, group_occs_sf, env_values_df); gc()
  } # End scenario loop for VIF/Corr
} else {
  cat ("Skipping VIF/Correlation analysis as run_vif_correlation is FALSE.\n")
  # Attempt to load previously suggested variables if skipping
  # This requires the user to have saved them appropriately before
  # temp_suggested_vars_path <- file.path(exp_log_dir, "suggested_vars_from_step_list.rds") # Example save path
  # if(file.exists(temp_suggested_vars_path)) {
  #     suggested_vars_from_step <- readRDS(temp_suggested_vars_path)
  #     cat("Loaded previously suggested variables.\n")
  # } else {
  #    warning("Cannot load suggested variables, list will be empty.")
  # }
}
# Optional: Save the list of suggested vars for reference
# saveRDS(suggested_vars_from_step, file.path(exp_log_dir, "suggested_vars_from_step_list.rds"))

cat("\n--- PART 1 FINISHED: Review VIF/Correlation plots in", exp_log_dir, "---\n")

#-------------------------------------------------------------------------------
# !!! --- MANUAL STEP: DEFINE FINAL VARIABLES FOR PCA --- !!!
#-------------------------------------------------------------------------------
# Based on the plots and suggested variable lists generated above (saved in
# exp_log_dir), define the final list of variables you want to use for PCA
# for EACH scenario in this experiment.
# Replace the example lists below with your actual choices.

final_selected_vars_experiment <- list()

# Example structure (Populate this for EACH scenario you ran VIF/Corr for):
final_selected_vars_experiment[['current']] <- c(
  'thetao_baseline_2000_2019_depthmax', # Replace with your chosen variables
  'so_baseline_2000_2019_depthmax',
  'distcoast',
  'rugosity',
  'bathymetry_mean'
  # ... add others based on your review
)
final_selected_vars_experiment[['ssp119_2050']] <- c(
  'thetao_mean', # Use cleaned base names matching stack output
  'so_mean',
  'distcoast',
  'rugosity',
  'bathymetry_mean'
  # ... add others based on your review
)
final_selected_vars_experiment[['ssp119_2100']] <- c(
  'thetao_mean',
  'so_mean',
  'distcoast',
  'rugosity',
  'bathymetry_mean'
  # ... add others based on your review
)
final_selected_vars_experiment[['ssp585_2050']] <- c(
  'thetao_mean',
  'so_mean',
  'distcoast',
  'rugosity',
  'bathymetry_mean'
  # ... add others based on your review
)
final_selected_vars_experiment[['ssp585_2100']] <- c(
  'thetao_mean',
  'so_mean',
  'distcoast',
  'rugosity',
  'bathymetry_mean'
  # ... add others based on your review
)
# Add entries for all scenarios you intend to model in this experiment!

# --- Safety Check: Ensure selections are made for needed scenarios ---
needed_scenarios <- config$env_scenarios # Could potentially filter this if needed
missing_selections <- setdiff(needed_scenarios, names(final_selected_vars_experiment))
if(length(missing_selections) > 0) {
  stop("FATAL: Manual variable selections are missing for the following scenarios in 'final_selected_vars_experiment': ",
       paste(missing_selections, collapse=", "), ". Please define them above and re-run.")
}

cat("\n--- MANUAL STEP COMPLETE: Using defined variable selections for PCA/SDM ---\n")
#-------------------------------------------------------------------------------


# --- 5. PCA & Projection Loop (Run this part after defining variables) ---
cat("\n--- PART 2: PCA Analysis & Projection (Anemones, Selected Vars) ---\n")
run_pca_projection <- TRUE # Set to FALSE to skip if already done

if(run_pca_projection) {
  for (scenario in needed_scenarios) { # Only loop scenarios with selections
    cat("\n-----------------------------------------------------\n")
    cat(" Running PCA & Projection for Scenario:", scenario, "\n")
    cat("-----------------------------------------------------\n")
    
    # Get the final selected variables for this scenario
    final_selected_vars <- final_selected_vars_experiment[[scenario]]
    if(is.null(final_selected_vars) || length(final_selected_vars) < 2) {
      warning("Skipping PCA for ", scenario, ": Fewer than 2 variables selected.")
      next
    }
    cat("  Using variables for PCA:", paste(final_selected_vars, collapse=", "), "\n")
    
    # --- Reload data (necessary as it was cleared in previous loop or if run separately) ---
    env_stack_raw <- load_stack_env_data(scenario, config)
    if (is.null(env_stack_raw)) {warning("Failed load env for ", scenario); next}
    env_stack_processed <- preprocess_env_rasters(env_stack_raw, config)
    if (is.null(env_stack_processed)) {warning("Failed preprocess env for ", scenario); next}
    rm(env_stack_raw); gc()
    target_crs <- terra::crs(env_stack_processed)
    group_occs_sf <- load_clean_thin_group_occurrences(anemone_occ_dir, target_crs, config, env_stack_processed)
    if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) {warning("No valid occurrences for ", scenario); next}
    env_values_df <- extract_env_values(group_occs_sf, env_stack_processed)
    if (is.null(env_values_df) || ncol(env_values_df) < 2) {warning("Insufficient env data extracted for ", scenario); next}
    # Check selected variables exist in the loaded stack AND extracted values
    if(!all(final_selected_vars %in% names(env_stack_processed)) || !all(final_selected_vars %in% names(env_values_df))) {
      warning("Selected PCA variables mismatch loaded data for scenario ", scenario, ". Skipping PCA.", call.=FALSE)
      next
    }
    env_values_df_selected <- env_values_df[, final_selected_vars, drop = FALSE]
    if (nrow(env_values_df_selected) < 2) {warning("Insufficient rows after selection for PCA in ", scenario); next}
    
    
    # --- Run PCA (Original Method) ---
    cat("  Running PCA (vegan::rda)...\n")
    pca_model_obj <- NULL
    tryCatch({
      env_std <- scale(env_values_df_selected, center = TRUE, scale = TRUE)
      pca_model_obj <- vegan::rda(env_std)
      
      # Save Plots with experimental suffix
      pca_biplot_path <- file.path(exp_log_dir, paste0("pca_biplot_", exp_group_name, "_", scenario, exp_suffix, ".png"))
      pca_contribution_path <- file.path(exp_log_dir, paste0("pca_contribution_", exp_group_name, "_", scenario, exp_suffix, ".png"))
      # Plotting code adapted from helper/original script 06
      pca_summary <- summary(pca_model_obj); prop_explained <- pca_summary$cont$importance[2, 1:2]
      site_scores <- as.data.frame(vegan::scores(pca_model_obj, display = "sites", choices = 1:2))
      species_scores <- as.data.frame(vegan::scores(pca_model_obj, display = "species", choices = 1:2))
      scaling_factor <- min((max(site_scores$PC1)-min(site_scores$PC1))/(max(species_scores$PC1)-min(species_scores$PC1)), (max(site_scores$PC2)-min(site_scores$PC2))/(max(species_scores$PC2)-min(species_scores$PC2)), na.rm=TRUE) * 0.9
      species_scores_scaled <- species_scores * scaling_factor
      pca_plot <- ggplot() + geom_point(data = site_scores, aes(x = PC1, y = PC2), color = "salmon", size = 1) + geom_segment(data = species_scores_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.25, "cm")), color = "black") + geom_text(data = species_scores_scaled, aes(x = PC1 * 1.1, y = PC2 * 1.1, label = rownames(species_scores_scaled)), hjust = 0.5, vjust = 0.5, size = 3, color = "black", check_overlap = TRUE) + labs(x = paste("PCA Axis 1 (", round(prop_explained[1] * 100, 2), "%)", sep = ""), y = paste("PCA Axis 2 (", round(prop_explained[2] * 100, 2), "%)", sep = ""), title = paste("PCA Biplot (Exp) -", exp_group_name, scenario)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      ggsave(filename = pca_biplot_path, plot = pca_plot, width = 8, height = 7)
      
      var_loadings <- vegan::scores(pca_model_obj, display = "species", choices = 1:config$n_pca_components); var_loadings_df <- as.data.frame(var_loadings); var_loadings_df$Variable <- rownames(var_loadings_df)
      var_loadings_long <- tidyr::pivot_longer(var_loadings_df, cols = starts_with("PC"), names_to = "Principal_Component", values_to = "Contribution"); var_loadings_long$Principal_Component <- factor(var_loadings_long$Principal_Component, levels = paste0("PC", 1:config$n_pca_components))
      num_vars <- length(unique(var_loadings_long$Variable)); colors <- grDevices::colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars)
      contribution_plot <- ggplot(var_loadings_long, aes(x = Principal_Component, y = Contribution, fill = Variable)) + geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = colors) + labs(title = paste("Variable Contribution (Exp) -", exp_group_name, scenario), x = "Principal Component", y = "Contribution") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
      ggsave(filename = pca_contribution_path, plot = contribution_plot, width = 10, height = 6)
      cat("  PCA plots saved to:", exp_log_dir, "\n")
      
      # Store model
      if (is.null(experimental_pca_models[[scenario]])) experimental_pca_models[[scenario]] <- list()
      experimental_pca_models[[scenario]] <- pca_model_obj
      
    }, error = function(e) { warning("PCA failed for ", scenario, ": ", e$message); pca_model_obj <- NULL })
    
    
    # --- Project PCA onto Rasters ---
    if (!is.null(pca_model_obj)) {
      cat("  Projecting PCA onto rasters...\n")
      tryCatch({
        env_stack_selected_for_pca <- env_stack_processed[[final_selected_vars]]
        # Use prcomp for prediction
        pca_prcomp <- stats::prcomp(env_values_df_selected, center = TRUE, scale. = TRUE)
        pca_rasters <- terra::predict(env_stack_selected_for_pca, pca_prcomp, index = 1:config$n_pca_components)
        names(pca_rasters) <- paste0("PC", 1:config$n_pca_components)
        
        # Save with experimental suffix
        pca_raster_file <- file.path(exp_log_dir, paste0("pca_rasters_", exp_group_name, "_", scenario, exp_suffix, ".tif"))
        terra::writeRaster(pca_rasters, filename = pca_raster_file, overwrite = TRUE, gdal=c("COMPRESS=LZW"))
        cat("  Experimental PCA raster layers saved to:", pca_raster_file, "\n")
        # Store path
        experimental_pca_paths[[scenario]] <- pca_raster_file
        
      }, error = function(e) { warning("PCA projection failed for ", scenario, ": ", e$message) })
    }
    rm(env_stack_processed, group_occs_sf, env_values_df, env_values_df_selected, pca_model_obj); gc()
  } # End scenario loop for PCA
} else {
  cat("Skipping PCA/Projection as run_pca_projection is FALSE.\n")
  # Attempt to load previously created paths if skipping
  temp_pca_paths_path <- file.path(exp_log_dir, "experimental_pca_paths.rds")
  if(file.exists(temp_pca_paths_path)) {
    experimental_pca_paths <- readRDS(temp_pca_paths_path)
    cat("Loaded previously generated experimental PCA paths.\n")
  } else {
    warning("Cannot load experimental PCA paths, SDM step may fail.")
  }
}
# Save experimental PCA paths and models
saveRDS(experimental_pca_paths, file.path(exp_log_dir, "experimental_pca_paths.rds"))
saveRDS(experimental_pca_models, file.path(exp_log_dir, "experimental_pca_models.rds"))

cat("\n--- PART 2 FINISHED: PCA and Projection complete for selected variables ---\n")

#-------------------------------------------------------------------------------
# --- 6. Run Standard SDMs Loop (Anemones Only, Using Experimental PCA Layers) ---
#-------------------------------------------------------------------------------
cat("\n--- PART 3: Running Standard SDMs (Anemones Only, Experimental Predictors) ---\n")

# Check if experimental PCA paths were generated/loaded
if(length(experimental_pca_paths) == 0) {
  stop("Experimental PCA raster paths are missing. Cannot run SDMs. Ensure Part 2 completed successfully.")
}

run_sdms <- TRUE # Set to FALSE to skip SDM runs

if(run_sdms) {
  total_sdms_run <- 0
  total_sdms_skipped <- 0
  
  for (i in 1:nrow(anemone_species)) {
    
    species_name_sanitized <- gsub(" ", "_", anemone_species$scientificName[i])
    species_aphia_id <- anemone_species$AphiaID[i]
    
    cat("\n-----------------------------------------------------\n")
    cat(" Processing Species:", anemone_species$scientificName[i], "\n")
    cat("-----------------------------------------------------\n")
    
    # Load/clean individual occurrences
    occ_sf_clean <- load_clean_individual_occ(species_aphia_id, anemone_occ_dir, config)
    if(is.null(occ_sf_clean) || nrow(occ_sf_clean) < config$min_occurrences_sdm) {
      warning("Skipping: Not enough occurrences (", nrow(occ_sf_clean), ").", call.=FALSE)
      total_sdms_skipped <- total_sdms_skipped + length(needed_scenarios)
      next
    }
    
    for (scenario in needed_scenarios) { # Loop only through scenarios with PCA rasters
      cat("  --- Scenario:", scenario, "---\n")
      
      # Define EXPERIMENTAL Output Paths
      pred_file <- file.path(exp_pred_dir, paste0("sdm_prediction_", species_name_sanitized, "_", scenario, exp_suffix, ".tif"))
      results_file <- file.path(exp_results_dir, paste0("sdm_results_", species_name_sanitized, "_", scenario, exp_suffix, ".rds"))
      eval_file <- file.path(exp_results_dir, paste0("sdm_eval_", species_name_sanitized, "_", scenario, exp_suffix, ".csv"))
      model_obj_file <- file.path(exp_models_dir, paste0("sdm_model_", species_name_sanitized, "_", scenario, exp_suffix, ".rds"))
      
      # Check skip conditions
      pca_stack_path <- experimental_pca_paths[[scenario]]
      if (is.null(pca_stack_path) || !file.exists(pca_stack_path)) {
        warning("Experimental PCA stack missing for scenario ", scenario, ". Skipping SDM.", call.=FALSE)
        total_sdms_skipped <- total_sdms_skipped + 1; next
      }
      if (file.exists(pred_file)) { # Simple check if already run for this experiment
        cat("    Prediction file exists. Skipping SDM run.\n")
        total_sdms_skipped <- total_sdms_skipped + 1; next
      }
      
      # Load EXPERIMENTAL PCA predictor stack
      predictor_stack <- tryCatch(terra::rast(pca_stack_path), error = function(e) NULL)
      if(is.null(predictor_stack)) { warning("Failed load PCA stack for ", scenario); total_sdms_skipped <- total_sdms_skipped + 1; next }
      
      # Thin occurrences
      occ_sf_thinned <- thin_individual_occ(occ_sf_clean, predictor_stack, config)
      if(is.null(occ_sf_thinned) || nrow(occ_sf_thinned) < config$min_occurrences_sdm) {
        warning("Skipping: Not enough occurrences after thinning (", nrow(occ_sf_thinned), ").", call.=FALSE)
        total_sdms_skipped <- total_sdms_skipped + 1; next
      }
      
      # Generate Background Points
      background_points <- generate_sdm_background(predictor_stack, config$background_points_n, config)
      if (is.null(background_points)) {warning("Failed background gen for ", scenario); total_sdms_skipped <- total_sdms_skipped + 1; next}
      
      # Run ENMeval
      enmeval_results <- run_sdm_enmeval(occ_sf_thinned, predictor_stack, background_points, config)
      if (is.null(enmeval_results)) {warning("ENMeval failed for ", scenario); total_sdms_skipped <- total_sdms_skipped + 1; next}
      tryCatch(saveRDS(enmeval_results, file = results_file), error=function(e){warning("Failed to save ENMeval results.")})
      tryCatch(readr::write_csv(enmeval_results@results, file = eval_file), error=function(e){warning("Failed to save Eval table.")})
      
      # Predict Best Model
      prediction_raster <- predict_sdm_best(enmeval_results, predictor_stack, config)
      if (is.null(prediction_raster)) {warning("Prediction failed for ", scenario); total_sdms_skipped <- total_sdms_skipped + 1; next}
      
      # Save Prediction Raster
      tryCatch({ terra::writeRaster(prediction_raster, filename = pred_file, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES")); cat("    Prediction raster saved.\n"); total_sdms_run <- total_sdms_run + 1 }, error=function(e){warning("Failed save prediction raster."); total_sdms_skipped <- total_sdms_skipped + 1})
      
      # Save Best Model Object
      tryCatch({eval_results_df<-enmeval_results@results; metric<-config$sdm_evaluation_metric; if(!metric %in% names(eval_results_df)) metric<-"AICc"; select_best<-if(metric=="AICc") which.min else which.max; valid_rows<-!is.na(eval_results_df[[metric]]); if(sum(valid_rows)>0){best_model_index<-select_best(eval_results_df[[metric]][valid_rows]); original_indices<-which(valid_rows); best_model_row_index<-original_indices[best_model_index]; best_model_obj<-enmeval_results@models[[best_model_row_index]]; saveRDS(best_model_obj, file=model_obj_file); cat("    Best model object saved.\n")}}, error=function(e){warning("Failed to save best model object: ", e$message)})
      
      rm(predictor_stack, occ_sf_thinned, background_points, enmeval_results, prediction_raster); gc()
    } # End scenario loop
    rm(occ_sf_clean); gc()
  } # End species loop
  
  cat("\n--- Experimental Anemone SDM runs finished ---\n")
  cat("Total SDMs successfully run:", total_sdms_run, "\n")
  cat("Total SDMs skipped:", total_sdms_skipped, "\n")
} else {
  cat("Skipping SDM runs as run_sdms is FALSE.\n")
}

cat("\n--- EXPERIMENTAL Script Finished ---\n")
#-------------------------------------------------------------------------------