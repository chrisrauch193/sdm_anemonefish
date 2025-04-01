# scripts/07_run_anemonefish_biotic_sdms.R
#-------------------------------------------------------------------------------
# Run Biotic SDMs for Anemonefish Species
#
# Workflow:
# 1. Load inputs: fish list, association list, paths to standard anemone SDM
#    predictions, paths to fish PCA rasters.
# 2. Loop through each anemonefish species.
# 3. Find its associated host anemone(s).
# 4. Loop through each environmental scenario.
# 5. Load relevant predictors:
#    - Anemone prediction raster(s) for the scenario.
#    - Anemonefish PCA raster stack for the scenario.
# 6. Run "Biotic-Only" SDM:
#    - Predictors: Anemone suitability map(s).
#    - Thin fish occurrences, generate background points based on anemone maps.
#    - Run ENMeval, predict, save results.
# 7. Run "Combined" SDM:
#    - Predictors: Anemone suitability map(s) + Fish PCA rasters.
#    - Thin fish occurrences, generate background points based on combined layers.
#    - Run ENMeval, predict, save results.
#-------------------------------------------------------------------------------
cat("--- Running Script 07: Run Anemonefish Biotic SDMs ---\n")

# Ensure config is loaded
if (!exists("config")) {
  source("scripts/config.R")
  if (!exists("config")) {
    stop("Failed to load config object from scripts/config.R")
  }
}

# Load necessary packages
pacman::p_load(terra, sf, dplyr, readr, ENMeval, tools, stringr)

# Source helper functions
sdm_helper_path <- file.path(config$helpers_dir, "sdm_modeling_helpers.R")
if (!file.exists(sdm_helper_path)) stop("SDM Helper file not found: ", sdm_helper_path)
source(sdm_helper_path)

# --- Load Inputs ---

# 1. Anemonefish Species List
tryCatch({
  anemonefish_species <- readr::read_csv(config$anemonefish_species_list_file, show_col_types = FALSE)
}, error = function(e) {
  stop("Failed to load anemonefish species list: ", config$anemonefish_species_list_file, ". Error: ", e$message)
})

# 2. Anemone-Anemonefish Association List
if (!file.exists(config$anemone_fish_association_file)) {
  stop("Anemone-Anemonefish association file not found: ", config$anemone_fish_association_file)
}
tryCatch({
  associations <- readr::read_csv(config$anemone_fish_association_file, show_col_types = FALSE)
  # Ensure required columns exist (adjust names if different in your file)
  required_cols <- c("AnemonefishScientificName", "AssociatedAnemoneScientificName")
  if (!all(required_cols %in% names(associations))) {
    stop("Association file must contain columns: ", paste(required_cols, collapse=", "))
  }
}, error = function(e) {
  stop("Failed to load association file: ", config$anemone_fish_association_file, ". Error: ", e$message)
})


# 3. Anemonefish PCA Raster Paths (from script 05)
if (!file.exists(config$pca_raster_paths_rds_path)) {
  stop("PCA raster paths file not found: ", config$pca_raster_paths_rds_path, ". Run script 05 first.")
}
pca_raster_paths_all <- readRDS(config$pca_raster_paths_rds_path)
# Check if anemonefish paths exist
if(is.null(pca_raster_paths_all[['anemonefish']])) {
  stop("Paths for anemonefish PCA rasters not found in ", config$pca_raster_paths_rds_path)
}
fish_pca_raster_paths <- pca_raster_paths_all[['anemonefish']]


# 4. Standard Anemone Prediction Raster Directory
anemone_pred_dir <- config$predictions_dir # Assuming standard predictions are here


# 5. Create output subdirectories for biotic models
pred_dir_biotic_only <- file.path(config$predictions_dir, "biotic_only")
pred_dir_combined <- file.path(config$predictions_dir, "biotic_combined")
results_dir_biotic_only <- file.path(config$results_dir, "biotic_only")
results_dir_combined <- file.path(config$results_dir, "biotic_combined")
models_dir_biotic_only <- file.path(config$models_dir, "biotic_only")
models_dir_combined <- file.path(config$models_dir, "biotic_combined")

dir.create(pred_dir_biotic_only, recursive = TRUE, showWarnings = FALSE)
dir.create(pred_dir_combined, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir_biotic_only, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir_combined, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir_biotic_only, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir_combined, recursive = TRUE, showWarnings = FALSE)


# --- Main Loop: Anemonefish Species -> Scenario ---

total_biotic_only_run <- 0
total_combined_run <- 0
total_biotic_skipped <- 0

for (i in 1:nrow(anemonefish_species)) {
  
  fish_sciname <- anemonefish_species$scientificName[i]
  fish_sciname_sanitized <- gsub(" ", "_", fish_sciname)
  fish_aphia_id <- anemonefish_species$AphiaID[i]
  
  cat("\n=========================================================\n")
  cat("Processing Fish:", fish_sciname, "(", i, "/", nrow(anemonefish_species), ")\n")
  
  # Find associated host anemones for this fish
  host_anemone_names <- associations %>%
    filter(AnemonefishScientificName == fish_sciname) %>%
    pull(AssociatedAnemoneScientificName) %>%
    unique()
  
  if (length(host_anemone_names) == 0) {
    warning("No associated host anemones found for ", fish_sciname, " in association file. Skipping biotic SDMs.", call.=FALSE)
    total_biotic_skipped <- total_biotic_skipped + (length(config$env_scenarios) * 2) # Skip both model types for all scenarios
    next # Skip to next fish species
  }
  cat("  Associated Host(s):", paste(host_anemone_names, collapse=", "), "\n")
  
  # Load and clean fish occurrences ONCE per species
  fish_occ_sf_clean <- load_clean_individual_occ(fish_aphia_id, config$anemonefish_occurrence_dir, config)
  
  if(is.null(fish_occ_sf_clean) || nrow(fish_occ_sf_clean) < config$min_occurrences_sdm) {
    warning("Skipping fish ", fish_sciname, ": Not enough occurrences after cleaning (", nrow(fish_occ_sf_clean), "<", config$min_occurrences_sdm, ").", call.=FALSE)
    total_biotic_skipped <- total_biotic_skipped + (length(config$env_scenarios) * 2)
    next # Skip to next fish species
  }
  cat("  Loaded and cleaned", nrow(fish_occ_sf_clean), "initial occurrences for fish.\n")
  
  
  # Inner Loop: Scenarios
  for (scenario in config$env_scenarios) {
    cat("\n-----------------------------------------------------\n")
    cat("  Processing Scenario:", scenario, "for Fish:", fish_sciname, "\n")
    cat("-----------------------------------------------------\n")
    
    # --- Load Scenario-Specific Predictors ---
    
    # 1. Load Host Anemone Prediction Raster(s)
    anemone_predictor_list <- list()
    missing_anemone_pred <- FALSE
    for (anemone_name in host_anemone_names) {
      anemone_name_sanitized <- gsub(" ", "_", anemone_name)
      anemone_pred_file <- file.path(anemone_pred_dir, paste0("sdm_prediction_", anemone_name_sanitized, "_", scenario, ".tif"))
      
      if (!file.exists(anemone_pred_file)) {
        warning("Host anemone prediction raster not found for: ", anemone_name, " scenario: ", scenario, "\n  File expected: ", anemone_pred_file, call.=FALSE)
        missing_anemone_pred <- TRUE
        # Decide how to handle: skip this host? skip the scenario? skip the fish?
        # For now, we'll skip the biotic models for this scenario if ANY host is missing.
        break # Exit the inner loop over anemones
      }
      anemone_pred_rast <- tryCatch(terra::rast(anemone_pred_file), error = function(e) NULL)
      if (is.null(anemone_pred_rast)) {
        warning("Failed to load host anemone prediction: ", anemone_pred_file, call.=FALSE)
        missing_anemone_pred <- TRUE
        break
      }
      # Standardize layer name
      names(anemone_pred_rast) <- paste0("anem_", anemone_name_sanitized)
      anemone_predictor_list[[anemone_name_sanitized]] <- anemone_pred_rast
    }
    
    if (missing_anemone_pred) {
      warning("Skipping biotic models for scenario ", scenario, " due to missing host anemone prediction(s).", call.=FALSE)
      total_biotic_skipped <- total_biotic_skipped + 2 # Skip both model types
      next # Skip to next scenario
    }
    
    # Stack if multiple hosts
    anemone_predictor_stack <- tryCatch(terra::rast(anemone_predictor_list), error=function(e)NULL)
    if(is.null(anemone_predictor_stack)){
      warning("Failed to stack anemone predictors for scenario ", scenario, call.=FALSE)
      total_biotic_skipped <- total_biotic_skipped + 2
      next
    }
    cat("    Loaded", length(names(anemone_predictor_stack)), "host anemone predictor(s).\n")
    
    
    # 2. Load Anemonefish PCA Raster Stack
    fish_pca_stack_path <- fish_pca_raster_paths[[scenario]]
    if (is.null(fish_pca_stack_path) || !file.exists(fish_pca_stack_path)) {
      warning("Anemonefish PCA raster stack not found for scenario: ", scenario, " at path: ", fish_pca_stack_path, ". Skipping Combined SDM.", call. = FALSE)
      fish_pca_stack <- NULL # Mark as unavailable
    } else {
      fish_pca_stack <- tryCatch(terra::rast(fish_pca_stack_path), error=function(e)NULL)
      if(is.null(fish_pca_stack)){
        warning("Failed to load fish PCA stack: ", fish_pca_stack_path, call.=FALSE)
      } else {
        cat("    Loaded fish PCA predictor stack:", basename(fish_pca_stack_path), "\n")
      }
    }
    
    
    # --- Run Biotic-Only SDM ---
    cat("    --- Running Biotic-Only Model ---\n")
    
    # Define output paths
    pred_file_bo <- file.path(pred_dir_biotic_only, paste0("sdm_prediction_biotic_only_", fish_sciname_sanitized, "_", scenario, ".tif"))
    results_file_bo <- file.path(results_dir_biotic_only, paste0("sdm_results_biotic_only_", fish_sciname_sanitized, "_", scenario, ".rds"))
    eval_file_bo <- file.path(results_dir_biotic_only, paste0("sdm_eval_biotic_only_", fish_sciname_sanitized, "_", scenario, ".csv"))
    model_file_bo <- file.path(models_dir_biotic_only, paste0("sdm_model_biotic_only_", fish_sciname_sanitized, "_", scenario, ".rds"))
    
    if (!config$force_rerun$run_biotic_sdms && file.exists(pred_file_bo)) {
      cat("      Biotic-only prediction exists. Skipping.\n")
      total_biotic_skipped <- total_biotic_skipped + 1
    } else {
      # Thin occurrences using ANEMONE predictor stack
      fish_occ_thinned_bo <- thin_individual_occ(fish_occ_sf_clean, anemone_predictor_stack, config)
      if (is.null(fish_occ_thinned_bo) || nrow(fish_occ_thinned_bo) < config$min_occurrences_sdm) {
        warning("Skipping Biotic-Only: Not enough occurrences after thinning (", nrow(fish_occ_thinned_bo), ").", call.=FALSE)
        total_biotic_skipped <- total_biotic_skipped + 1
      } else {
        cat("      Thinned occurrences for Biotic-Only:", nrow(fish_occ_thinned_bo), "points.\n")
        # Generate background points using ANEMONE predictor stack
        background_points_bo <- generate_sdm_background(anemone_predictor_stack, config$background_points_n, config)
        if (is.null(background_points_bo)) {
          warning("Failed to generate background points for Biotic-Only. Skipping.", call. = FALSE)
          total_biotic_skipped <- total_biotic_skipped + 1
        } else {
          # Run ENMeval
          enmeval_results_bo <- run_sdm_enmeval(fish_occ_thinned_bo, anemone_predictor_stack, background_points_bo, config)
          if (!is.null(enmeval_results_bo)) {
            tryCatch(saveRDS(enmeval_results_bo, file = results_file_bo), error=function(e){warning("Failed to save ENMeval results (Biotic-Only).")})
            tryCatch(readr::write_csv(enmeval_results_bo@results, file = eval_file_bo), error=function(e){warning("Failed to save Eval table (Biotic-Only).")})
            cat("      ENMeval results saved (Biotic-Only).\n")
            
            # Predict
            prediction_raster_bo <- predict_sdm_best(enmeval_results_bo, anemone_predictor_stack, config)
            if (!is.null(prediction_raster_bo)) {
              tryCatch({
                terra::writeRaster(prediction_raster_bo, filename = pred_file_bo, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
                cat("      Prediction raster saved (Biotic-Only):", basename(pred_file_bo), "\n")
                total_biotic_only_run <- total_biotic_only_run + 1
              }, error=function(e){warning("Failed to save prediction raster (Biotic-Only).")})
              
              # Save best model object
              tryCatch({
                eval_results_df<-enmeval_results_bo@results; metric<-config$sdm_evaluation_metric; if(!metric %in% names(eval_results_df)) metric<-"AICc"; select_best<-if(metric=="AICc") which.min else which.max; valid_rows<-!is.na(eval_results_df[[metric]]); if(sum(valid_rows)>0){best_model_index<-select_best(eval_results_df[[metric]][valid_rows]); original_indices<-which(valid_rows); best_model_row_index<-original_indices[best_model_index]; best_model_obj<-enmeval_results_bo@models[[best_model_row_index]]; saveRDS(best_model_obj, file=model_file_bo); cat("      Best model object saved (Biotic-Only).\n")} # Condense selection logic
              }, error=function(e){warning("Failed to save best model object (Biotic-Only): ", e$message)})
              
            } else { total_biotic_skipped <- total_biotic_skipped + 1 }
            rm(prediction_raster_bo)
          } else { total_biotic_skipped <- total_biotic_skipped + 1 } # ENMeval failed
          rm(enmeval_results_bo)
        } # Background points generated
        rm(background_points_bo)
      } # Enough points after thinning
      rm(fish_occ_thinned_bo)
    } # Check if file exists
    
    
    # --- Run Combined Biotic + Abiotic SDM ---
    cat("    --- Running Combined Biotic + Abiotic Model ---\n")
    
    # Define output paths
    pred_file_comb <- file.path(pred_dir_combined, paste0("sdm_prediction_combined_", fish_sciname_sanitized, "_", scenario, ".tif"))
    results_file_comb <- file.path(results_dir_combined, paste0("sdm_results_combined_", fish_sciname_sanitized, "_", scenario, ".rds"))
    eval_file_comb <- file.path(results_dir_combined, paste0("sdm_eval_combined_", fish_sciname_sanitized, "_", scenario, ".csv"))
    model_file_comb <- file.path(models_dir_combined, paste0("sdm_model_combined_", fish_sciname_sanitized, "_", scenario, ".rds"))
    
    if(is.null(fish_pca_stack)){
      cat("      Skipping Combined model: Fish PCA stack unavailable for this scenario.\n")
      total_biotic_skipped <- total_biotic_skipped + 1
    } else if (!config$force_rerun$run_biotic_sdms && file.exists(pred_file_comb)) {
      cat("      Combined prediction exists. Skipping.\n")
      total_biotic_skipped <- total_biotic_skipped + 1
    } else {
      # Combine predictors
      cat("      Combining anemone and fish PCA predictors...\n")
      combined_predictor_stack <- tryCatch({
        # Ensure CRS match before stacking (should match if processed correctly)
        if(terra::crs(anemone_predictor_stack) != terra::crs(fish_pca_stack)){
          warning("CRS mismatch between anemone predictors and fish PCA stack. Attempting to project fish PCA stack.", call.=FALSE)
          fish_pca_stack <- terra::project(fish_pca_stack, terra::crs(anemone_predictor_stack))
        }
        # Ensure extents align (resample fish PCA to anemone predictors extent/res)
        if (!terra::compareGeom(anemone_predictor_stack, fish_pca_stack, stopOnError=FALSE, messages=TRUE)) {
          cat("      Resampling fish PCA stack to match anemone predictors...\n")
          # Use first layer of anemone stack as template
          fish_pca_stack <- terra::resample(fish_pca_stack, anemone_predictor_stack[[1]], method="bilinear")
        }
        c(anemone_predictor_stack, fish_pca_stack) # Combine
      }, error = function(e){
        warning("Failed to combine predictor stacks: ", e$message, call.=FALSE)
        NULL
      })
      
      if(is.null(combined_predictor_stack)){
        warning("Skipping Combined model due to error combining predictors.", call.=FALSE)
        total_biotic_skipped <- total_biotic_skipped + 1
      } else {
        cat("      Combined stack created with layers:", paste(names(combined_predictor_stack), collapse=", "), "\n")
        # Thin occurrences using COMBINED predictor stack
        fish_occ_thinned_comb <- thin_individual_occ(fish_occ_sf_clean, combined_predictor_stack, config)
        if (is.null(fish_occ_thinned_comb) || nrow(fish_occ_thinned_comb) < config$min_occurrences_sdm) {
          warning("Skipping Combined: Not enough occurrences after thinning (", nrow(fish_occ_thinned_comb), ").", call.=FALSE)
          total_biotic_skipped <- total_biotic_skipped + 1
        } else {
          cat("      Thinned occurrences for Combined:", nrow(fish_occ_thinned_comb), "points.\n")
          # Generate background points using COMBINED predictor stack
          background_points_comb <- generate_sdm_background(combined_predictor_stack, config$background_points_n, config)
          if (is.null(background_points_comb)) {
            warning("Failed to generate background points for Combined. Skipping.", call. = FALSE)
            total_biotic_skipped <- total_biotic_skipped + 1
          } else {
            # Run ENMeval
            enmeval_results_comb <- run_sdm_enmeval(fish_occ_thinned_comb, combined_predictor_stack, background_points_comb, config)
            if (!is.null(enmeval_results_comb)) {
              tryCatch(saveRDS(enmeval_results_comb, file = results_file_comb), error=function(e){warning("Failed to save ENMeval results (Combined).")})
              tryCatch(readr::write_csv(enmeval_results_comb@results, file = eval_file_comb), error=function(e){warning("Failed to save Eval table (Combined).")})
              cat("      ENMeval results saved (Combined).\n")
              
              # Predict
              prediction_raster_comb <- predict_sdm_best(enmeval_results_comb, combined_predictor_stack, config)
              if (!is.null(prediction_raster_comb)) {
                tryCatch({
                  terra::writeRaster(prediction_raster_comb, filename = pred_file_comb, overwrite = TRUE, gdal=c("COMPRESS=LZW", "TFW=YES"))
                  cat("      Prediction raster saved (Combined):", basename(pred_file_comb), "\n")
                  total_combined_run <- total_combined_run + 1
                }, error=function(e){warning("Failed to save prediction raster (Combined).")})
                
                # Save best model object
                tryCatch({
                  eval_results_df<-enmeval_results_comb@results; metric<-config$sdm_evaluation_metric; if(!metric %in% names(eval_results_df)) metric<-"AICc"; select_best<-if(metric=="AICc") which.min else which.max; valid_rows<-!is.na(eval_results_df[[metric]]); if(sum(valid_rows)>0){best_model_index<-select_best(eval_results_df[[metric]][valid_rows]); original_indices<-which(valid_rows); best_model_row_index<-original_indices[best_model_index]; best_model_obj<-enmeval_results_comb@models[[best_model_row_index]]; saveRDS(best_model_obj, file=model_file_comb); cat("      Best model object saved (Combined).\n")}
                }, error=function(e){warning("Failed to save best model object (Combined): ", e$message)})
                
              } else { total_biotic_skipped <- total_biotic_skipped + 1 }
              rm(prediction_raster_comb)
            } else { total_biotic_skipped <- total_biotic_skipped + 1 } # ENMeval failed
            rm(enmeval_results_comb)
          } # Background points generated
          rm(background_points_comb)
        } # Enough points after thinning
        rm(fish_occ_thinned_comb)
      } # Predictors combined ok
      rm(combined_predictor_stack)
    } # Check if file exists / PCA available
    
    
    # Clean up scenario-specific rasters
    rm(anemone_predictor_stack, fish_pca_stack); gc()
    
  } # End scenario loop
  
  # Clean up species data
  rm(fish_occ_sf_clean); gc()
  
} # End fish species loop


cat("\n=========================================================\n")
cat("Biotic SDM runs finished.\n")
cat("Total Biotic-Only models successfully run:", total_biotic_only_run, "\n")
cat("Total Combined models successfully run:", total_combined_run, "\n")
cat("Total Biotic/Combined SDMs skipped:", total_biotic_skipped, "\n")
cat("=========================================================\n")

cat("\n--- Script 07 finished. ---\n")
#-------------------------------------------------------------------------------