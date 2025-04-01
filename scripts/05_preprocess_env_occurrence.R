# scripts/05_preprocess_env_occurrence.R
#-------------------------------------------------------------------------------
# Preprocess Environmental Data and Occurrences for VIF/PCA Analysis
# Uses methods from user's original scripts (car::vif, step, corrplot, vegan::rda)
# Loops through scenarios and groups, requires manual selection input.
#-------------------------------------------------------------------------------
cat("--- Running Script 05: Preprocess Env Data & Occurrences (Original Methods) ---\n")

# Ensure config is loaded
if (!exists("config")) {
  source("scripts/config.R")
  if (!exists("config")) {
    stop("Failed to load config object from scripts/config.R")
  }
}

# Load necessary packages
pacman::p_load(terra, sf, dplyr, readr, corrplot, ggplot2, car, vegan, ggvegan, tools, Hmisc) # Added car, Hmisc

# Source helper functions (including plotting)
helper_path <- file.path(config$helpers_dir, "env_processing_helpers.R")
if (!file.exists(helper_path)) stop("Helper file not found: ", helper_path)
source(helper_path)

# --- Define Groups ---
species_groups <- list(
  anemone = list(occurrence_dir = config$anemone_occurrence_dir),
  anemonefish = list(occurrence_dir = config$anemonefish_occurrence_dir)
)

# --- Initialize lists to store results across scenarios/groups ---
selected_variables_for_pca <- list() # Stores FINAL user selections for PCA
pca_raster_paths <- list()           # Stores paths to generated PCA rasters
pca_models <- list()                 # Stores fitted PCA models (rda objects)

# --- Check if previous results exist and load if not forcing rerun ---
if (!config$force_rerun$preprocess_env_occurrence) {
  if (file.exists(config$selected_vars_rds_path)) {
    cat("Loading existing selected variables from:", config$selected_vars_rds_path, "\n")
    selected_variables_for_pca <- readRDS(config$selected_vars_rds_path)
  }
  if (file.exists(config$pca_raster_paths_rds_path)) {
    cat("Loading existing PCA raster paths from:", config$pca_raster_paths_rds_path, "\n")
    pca_raster_paths <- readRDS(config$pca_raster_paths_rds_path)
  }
  if (file.exists(config$pca_models_rds_path)) {
    cat("Loading existing PCA models from:", config$pca_models_rds_path, "\n")
    pca_models <- readRDS(config$pca_models_rds_path)
  }
}

# --- Main Loop: Scenario -> Group ---
for (scenario in config$env_scenarios) {
  cat("\n=========================================================\n")
  cat("Processing Scenario:", scenario, "\n")
  cat("=========================================================\n")
  
  # 1. Load and Preprocess Environmental Data for the Scenario
  env_stack_raw <- load_stack_env_data(scenario, config)
  if (is.null(env_stack_raw)) {
    warning("Failed to load env data for scenario: ", scenario, ". Skipping scenario.", call. = FALSE)
    next
  }
  env_stack_processed <- preprocess_env_rasters(env_stack_raw, config)
  if (is.null(env_stack_processed)) {
    warning("Failed to preprocess env data for scenario: ", scenario, ". Skipping scenario.", call. = FALSE)
    next
  }
  if(!identical(env_stack_raw, env_stack_processed)) rm(env_stack_raw); gc()
  target_crs <- terra::crs(env_stack_processed)
  
  for (group_name in names(species_groups)) {
    cat("\n-----------------------------------------------------\n")
    cat("Processing Group:", group_name, "for Scenario:", scenario, "\n")
    cat("-----------------------------------------------------\n")
    
    group_info <- species_groups[[group_name]]
    group_occ_dir <- group_info$occurrence_dir
    
    # Define output file paths for VIF/Corr plots and suggested vars
    vif_plot_path <- file.path(config$log_dir, paste0("vif_plot_", group_name, "_", scenario, ".png"))
    cor_plot_path <- file.path(config$log_dir, paste0("correlation_plot_", group_name, "_", scenario, ".png"))
    suggested_vars_path <- file.path(config$log_dir, paste0("suggested_vars_after_step_", group_name, "_", scenario, ".txt")) # Renamed
    # PCA related paths
    pca_biplot_path <- file.path(config$log_dir, paste0("pca_biplot_", group_name, "_", scenario, ".png"))
    pca_contribution_path <- file.path(config$log_dir, paste0("pca_contribution_", group_name, "_", scenario, ".png"))
    pca_raster_file <- file.path(config$log_dir, paste0("pca_rasters_", group_name, "_", scenario, ".tif"))
    
    # Check if results exist and skip if not forcing
    skip_this_combo <- FALSE
    if (!config$force_rerun$preprocess_env_occurrence &&
        !is.null(pca_raster_paths[[group_name]][[scenario]]) &&
        file.exists(pca_raster_paths[[group_name]][[scenario]])) {
      cat("PCA raster path already recorded for", group_name, "-", scenario, ". Skipping VIF/PCA steps.\n")
      skip_this_combo <- TRUE
    }
    
    # --- Steps 2-6: Occurrences, Extraction, VIF/Pearson (Original Methods) ---
    env_values_df <- NULL # Initialize
    step_selected_vars <- NULL # Initialize vars selected by step()
    if (!skip_this_combo) {
      # 2. Load, Clean, Thin Group Occurrences
      group_occs_sf <- load_clean_thin_group_occurrences(group_occ_dir, target_crs, config, env_stack_processed)
      if (is.null(group_occs_sf) || nrow(group_occs_sf) == 0) {
        warning("No valid occurrences found for group: ", group_name, ". Skipping VIF/PCA.", call. = FALSE)
        next
      }
      
      # 3. Extract Environmental Values
      env_values_df <- extract_env_values(group_occs_sf, env_stack_processed)
      if (is.null(env_values_df) || nrow(env_values_df) < 2 || ncol(env_values_df) < 2) {
        warning("Insufficient env data extracted for group: ", group_name, ". Rows:", nrow(env_values_df), " Cols:", ncol(env_values_df), ". Skipping VIF/PCA.", call. = FALSE)
        next
      }
      
      # Check for zero-variance columns
      variances <- apply(env_values_df, 2, var, na.rm = TRUE)
      zero_var_cols <- names(variances[variances == 0])
      if (length(zero_var_cols) > 0) {
        cat("  Warning: Removing zero-variance columns:", paste(zero_var_cols, collapse=", "), "\n")
        env_values_df <- env_values_df[, !(names(env_values_df) %in% zero_var_cols), drop = FALSE]
      }
      if (ncol(env_values_df) < 2) {
        warning("Less than 2 variables remaining after removing zero variance. Skipping VIF/PCA.", call. = FALSE)
        next
      }
      
      # 4. Run VIF Analysis (Using car::vif and step as in original script 05)
      cat("\n  Running VIF analysis (using lm, car::vif, step)...\n")
      tryCatch({
        # Create dummy response variable (as in original script)
        present_dummy <- rep(1, nrow(env_values_df))
        lm_data <- cbind(present_dummy, env_values_df)
        
        # Initial full model
        full_model <- stats::lm(present_dummy ~ ., data = lm_data)
        
        # Optional: Calculate initial VIF values on full model
        initial_vif <- car::vif(full_model)
        cat("  Initial VIF values (max =", round(max(initial_vif, na.rm=TRUE),1), "):\n")
        # print(sort(initial_vif, decreasing=TRUE))
        
        # Use step() for variable selection based on AIC (as in original script)
        # Note: step() is primarily for model selection (AIC/BIC), not strictly multicollinearity.
        # It might remove variables useful for PCA. Using vifstep is often preferred for collinearity.
        # But sticking to original method here.
        cat("  Running step() for variable selection...\n")
        stepped_model <- stats::step(full_model, direction = "both", trace = 0) # trace=0 for less output
        step_selected_vars_lm <- names(coefficients(stepped_model))[-1] # Get predictor names (remove intercept)
        # Clean names if lm added backticks ` `
        step_selected_vars_lm <- gsub("^`|`$", "", step_selected_vars_lm)
        
        
        if (length(step_selected_vars_lm) < 2) {
          warning("step() selected less than 2 variables. VIF/Correlation cannot be calculated reliably.", call.=FALSE)
          step_selected_vars <- step_selected_vars_lm # Store the limited selection
        } else {
          # Calculate VIF *on the variables selected by step*
          vif_on_stepped <- car::vif(stepped_model)
          cat("  VIF values for variables selected by step():\n")
          print(sort(vif_on_stepped, decreasing=TRUE))
          
          # Use the helper function to plot these VIF values
          plot_vif_results_original(vif_on_stepped, save_path = vif_plot_path)
          
          # Check if any VIF values in the stepped model are still > threshold
          if(any(vif_on_stepped > config$vif_threshold)) {
            cat("  Warning: Some variables selected by step() still have VIF >", config$vif_threshold, "\n")
          }
          step_selected_vars <- step_selected_vars_lm # Store the selection
        }
        
      }, error = function(e) {
        warning("VIF/step analysis failed for ", group_name, "-", scenario, ": ", e$message, "\n  Proceeding with all variables, but caution advised.", call. = FALSE)
        step_selected_vars <- names(env_values_df) # Use all if failed
      })
      
      # Ensure we have at least 2 variables for correlation/PCA
      if(is.null(step_selected_vars) || length(step_selected_vars) < 2) {
        warning("Fewer than 2 variables selected after VIF/step. Skipping Correlation/PCA.", call. = FALSE)
        next
      }
      cat("  Variables selected by step():", paste(step_selected_vars, collapse=", "), "\n")
      
      
      # 5. Run Pearson Correlation (Using original corrplot approach)
      cat("\n  Running Pearson Correlation analysis on step-selected variables...\n")
      # Use the helper function which replicates the original corrplot logic
      plot_correlation_results_original(env_values_df[, step_selected_vars, drop=FALSE], save_path = cor_plot_path)
      
      # Identify highly correlated pairs (for user review)
      cor_matrix <- cor(env_values_df[, step_selected_vars, drop = FALSE], method = "pearson", use = "pairwise.complete.obs")
      cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA
      high_cor_pairs <- which(abs(cor_matrix) > config$correlation_threshold, arr.ind = TRUE)
      if(nrow(high_cor_pairs) > 0){
        cat("  Highly correlated pairs (|r| >", config$correlation_threshold, ") among step-selected variables (for user review):\n")
        for(i in 1:nrow(high_cor_pairs)) {
          cat("    -", rownames(cor_matrix)[high_cor_pairs[i,1]], "&", colnames(cor_matrix)[high_cor_pairs[i,2]],
              "(r =", round(cor_matrix[high_cor_pairs[i,1], high_cor_pairs[i,2]], 2), ")\n")
        }
      } else {
        cat("  No highly correlated pairs found among step-selected variables.\n")
      }
      
      # 6. Save suggested variables (output from step)
      writeLines(step_selected_vars, suggested_vars_path)
      cat("  Step-selected variables saved for review:", suggested_vars_path, "\n")
      cat("  Correlation plot saved for review:", cor_plot_path, "\n")
      cat("  VIF plot saved for review:", vif_plot_path, "\n")
      
    } # End if(!skip_this_combo) for VIF/Corr
    
    
    # --- MANUAL INTERVENTION POINT ---
    cat("\n--- MANUAL STEP REQUIRED for:", group_name, "-", scenario, "---\n")
    cat("1. Review the VIF plot saved in:\n  ", vif_plot_path, "\n")
    cat("2. Review the correlation plot saved in:\n  ", cor_plot_path, "\n")
    cat("3. Review the variables selected by step() saved in:\n  ", suggested_vars_path, "\n")
    cat("4. Decide on the FINAL list of variables for PCA based on:\n")
    cat("    - Variables selected by step().\n")
    cat("    - Removing one variable from any highly correlated pairs (|r| >", config$correlation_threshold, ").\n")
    cat("    - Ecological relevance.\n")
    cat("5. UPDATE the 'selected_variables_for_pca' list within this script OR save/edit '", config$selected_vars_rds_path, "' and reload.\n", sep="")
    cat("   Example R code:\n")
    cat("   selected_variables_for_pca[['", group_name, "']][['", scenario, "']] <- c('var1', 'var3', 'var5', ...)\n", sep="")
    cat("6. Re-run this script (05). It will detect the selection and proceed with PCA/Projection.\n")
    cat("------------------------------------------------------\n")
    
    # --- Check if final selection exists ---
    final_selected_vars <- NULL
    if (!is.null(selected_variables_for_pca[[group_name]][[scenario]])) {
      final_selected_vars <- selected_variables_for_pca[[group_name]][[scenario]]
      cat("Found existing final selection for", group_name, "-", scenario, ":", paste(final_selected_vars, collapse=", "), "\n")
    } else {
      cat("No final variable selection found for", group_name, "-", scenario, ". Stopping processing for this combo.\n")
      # Save current state before stopping/continuing
      saveRDS(selected_variables_for_pca, config$selected_vars_rds_path)
      saveRDS(pca_raster_paths, config$pca_raster_paths_rds_path)
      saveRDS(pca_models, config$pca_models_rds_path)
      next # Skip to next group/scenario
    }
    
    # --- Ensure data is loaded if VIF/Corr was skipped ---
    if(is.null(env_values_df)) {
      cat("Reloading occurrences and extracting env values for PCA...\n")
      group_occs_sf <- load_clean_thin_group_occurrences(group_occ_dir, target_crs, config, env_stack_processed)
      if(is.null(group_occs_sf)) next
      env_values_df <- extract_env_values(group_occs_sf, env_stack_processed)
      if(is.null(env_values_df)) next
    }
    # Check selected variables exist in extracted data
    if(!all(final_selected_vars %in% names(env_values_df))){
      missing_vars <- final_selected_vars[!final_selected_vars %in% names(env_values_df)]
      warning("Some selected variables for PCA are not present in the extracted data for ", group_name, "-", scenario, ": ", paste(missing_vars, collapse=", "), ". Skipping PCA.", call.=FALSE)
      next
    }
    env_values_df_selected <- env_values_df[, final_selected_vars, drop = FALSE]
    if (is.null(env_values_df_selected) || nrow(env_values_df_selected) < 2 || ncol(env_values_df_selected) < 2) {
      warning("Insufficient data for PCA after final selection. Rows:", nrow(env_values_df_selected), " Cols:", ncol(env_values_df_selected), ". Skipping PCA.", call. = FALSE)
      next
    }
    
    # --- Steps 7-9: PCA (Original Method) and Projection ---
    cat("\n  Running PCA (vegan::rda) on final selected variables...\n")
    pca_model_obj <- NULL
    tryCatch({
      # 7. Standardize and Run PCA (using vegan::rda as in original script 06)
      env_std <- scale(env_values_df_selected, center = TRUE, scale = TRUE)
      pca_model_obj <- vegan::rda(env_std) # Use rda as original
      
      # --- Plot PCA Biplot (using ggplot approach from original script 06) ---
      pca_summary <- summary(pca_model_obj)
      prop_explained <- pca_summary$cont$importance[2, 1:2] # Prop explained by PC1, PC2
      site_scores <- as.data.frame(vegan::scores(pca_model_obj, display = "sites", choices = 1:2))
      species_scores <- as.data.frame(vegan::scores(pca_model_obj, display = "species", choices = 1:2))
      # Optional scaling (as in original script 06)
      scaling_factor <- min(
        (max(site_scores$PC1)-min(site_scores$PC1))/(max(species_scores$PC1)-min(species_scores$PC1)),
        (max(site_scores$PC2)-min(site_scores$PC2))/(max(species_scores$PC2)-min(species_scores$PC2))
      ) * 0.9
      species_scores_scaled <- species_scores * scaling_factor
      
      pca_plot <- ggplot() +
        geom_point(data = site_scores, aes(x = PC1, y = PC2), color = "salmon", size = 1) +
        geom_segment(data = species_scores_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.25, "cm")), color = "black") +
        geom_text(data = species_scores_scaled, aes(x = PC1 * 1.1, y = PC2 * 1.1, label = rownames(species_scores_scaled)), # Adjust label position
                  hjust = 0.5, vjust = 0.5, size = 3, color = "black", check_overlap = TRUE) +
        labs(x = paste("PCA Axis 1 (", round(prop_explained[1] * 100, 2), "%)", sep = ""),
             y = paste("PCA Axis 2 (", round(prop_explained[2] * 100, 2), "%)", sep = ""),
             title = paste("PCA Biplot -", group_name, scenario)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      ggsave(filename = pca_biplot_path, plot = pca_plot, width = 8, height = 7)
      cat("  PCA biplot saved to:", pca_biplot_path, "\n")
      
      # --- Plot Variable Contributions (using ggplot approach from original script 06) ---
      var_loadings <- vegan::scores(pca_model_obj, display = "species", choices = 1:config$n_pca_components)
      var_loadings_df <- as.data.frame(var_loadings)
      var_loadings_df$Variable <- rownames(var_loadings_df)
      var_loadings_long <- tidyr::pivot_longer(var_loadings_df, cols = starts_with("PC"),
                                               names_to = "Principal_Component", values_to = "Contribution")
      var_loadings_long$Principal_Component <- factor(var_loadings_long$Principal_Component, levels = paste0("PC", 1:config$n_pca_components))
      num_vars <- length(unique(var_loadings_long$Variable))
      colors <- grDevices::colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars)
      
      contribution_plot <- ggplot(var_loadings_long, aes(x = Principal_Component, y = Contribution, fill = Variable)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = colors) +
        labs(title = paste("Variable Contribution to PCs -", group_name, scenario),
             x = "Principal Component", y = "Contribution") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
      ggsave(filename = pca_contribution_path, plot = contribution_plot, width = 10, height = 6)
      cat("  PCA contribution plot saved to:", pca_contribution_path, "\n")
      
      # Store the vegan::rda model object
      if (is.null(pca_models[[group_name]])) pca_models[[group_name]] <- list()
      pca_models[[group_name]][[scenario]] <- pca_model_obj
      cat("  PCA model (vegan::rda) stored.\n")
      
    }, error = function(e) {
      warning("PCA analysis (vegan::rda) failed for ", group_name, "-", scenario, ": ", e$message, call. = FALSE)
      pca_model_obj <- NULL
    })
    
    
    # 8. Project PCA onto Rasters
    if (!is.null(pca_model_obj)) {
      cat("\n  Projecting PCA onto selected rasters...\n")
      tryCatch({
        if(!all(final_selected_vars %in% names(env_stack_processed))) {
          stop("Selected variables used for PCA do not match layers in the processed raster stack.")
        }
        env_stack_selected_for_pca <- env_stack_processed[[final_selected_vars]]
        
        # Use stats::prcomp for prediction compatibility with terra::predict
        pca_prcomp <- stats::prcomp(env_values_df_selected, center = TRUE, scale. = TRUE)
        pca_rasters <- terra::predict(env_stack_selected_for_pca, pca_prcomp, index = 1:config$n_pca_components)
        names(pca_rasters) <- paste0("PC", 1:config$n_pca_components)
        
        terra::writeRaster(pca_rasters, filename = pca_raster_file, overwrite = TRUE, gdal=c("COMPRESS=LZW"))
        cat("  PCA raster layers saved to:", pca_raster_file, "\n")
        
        if (is.null(pca_raster_paths[[group_name]])) pca_raster_paths[[group_name]] <- list()
        pca_raster_paths[[group_name]][[scenario]] <- pca_raster_file
        
      }, error = function(e) {
        warning("PCA projection failed for ", group_name, "-", scenario, ": ", e$message, call. = FALSE)
      })
    } # End if PCA model exists
    
    rm(group_occs_sf, env_values_df, env_values_df_selected, pca_model_obj); gc()
    
  } # End loop groups
  rm(env_stack_processed); gc()
} # End loop scenarios


# --- Save the final results ---
cat("\n--- Saving final selections and paths ---\n")
tryCatch({
  saveRDS(selected_variables_for_pca, config$selected_vars_rds_path)
  cat("Final selected variables saved to:", config$selected_vars_rds_path, "\n")
}, error=function(e){warning("Error saving selected_variables_for_pca: ", e$message)})

tryCatch({
  saveRDS(pca_raster_paths, config$pca_raster_paths_rds_path)
  cat("PCA raster paths saved to:", config$pca_raster_paths_rds_path, "\n")
}, error=function(e){warning("Error saving pca_raster_paths: ", e$message)})

tryCatch({
  saveRDS(pca_models, config$pca_models_rds_path)
  cat("PCA models (vegan::rda objects) saved to:", config$pca_models_rds_path, "\n")
}, error=function(e){warning("Error saving pca_models: ", e$message)})


cat("\n--- Script 05 finished. Check log directory for plots and suggested variables. Update selections if needed and re-run. ---\n")
#-------------------------------------------------------------------------------