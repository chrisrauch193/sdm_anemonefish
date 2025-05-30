# scripts/sdm_runs/experiments/sdmtune__combined_tweatking_run.R
#-------------------------------------------------------------------------------
# Combined Tweaking Script: Runs SDMtune, then BIOMOD2 with multiple algorithms.
# - MAXNET uses SDMtune hyperparameters if available.
# - Other algorithms use 'bigboss' defaults.
# - Builds an ensemble model from best individual models.
# For a single species and current scenario.
#-------------------------------------------------------------------------------

cat("--- Running Script sdm_combined_tweaking_run.R (Multi-Algo Ensemble) ---\n")

# --- 1. Setup: Load Config FIRST ---
if (file.exists("scripts/config.R")) {
  source("scripts/config.R")
  if (!exists("config") || !is.list(config)) {
    stop("FATAL: 'config' list object not found or invalid after sourcing scripts/config.R")
  }
} else {
  stop("FATAL: Configuration file 'scripts/config.R' not found.")
}

# --- 2. Load Required Packages ---
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr, 
               biomod2, maxnet, blockCV, 
               presenceabsence, randomForest, gbm, mda, gam, earth, xgboost 
)

# --- 3. Source Helper Functions ---
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))

# --- 4. Define Group Specifics & Predictor Type ---
cat(paste(Sys.time(), "[INFO] --- Starting Script sdm_combined_tweaking_run.R (Multi-Algo Ensemble) ---\n"))

group_name <- "anemone" 
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors 
predictor_type_suffix_sdmtune <- ifelse(use_pca, "_pca", "_vif") 

cat(paste(Sys.time(), "[INFO] --- Processing Group:", group_name, "---\n"))
cat(paste(Sys.time(), "[INFO] --- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---\n"))

# --- 5. Load Predictor Information (for the 'current' scenario) ---
current_scenario_name <- "current"
env_predictor_stack_current <- NULL 
predictor_paths_or_list_sdmtune <- NULL 

if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (!file.exists(pca_paths_rds)) stop("PCA raster paths RDS file not found.")
  predictor_paths_or_list_sdmtune <- readRDS(pca_paths_rds) 
  current_pca_path <- predictor_paths_or_list_sdmtune[[current_scenario_name]]
  if (is.null(current_pca_path) || !file.exists(current_pca_path)) stop("PCA stack for 'current' scenario not found.")
  env_predictor_stack_current <- tryCatch(terra::rast(current_pca_path), error = function(e) {
    cat(paste(Sys.time(), "[ERROR] Failed to load PCA stack:", e$message, "\n")); NULL
  })
} else {
  core_vif_vars <- config$final_vars_vif_anemone 
  predictor_paths_or_list_sdmtune <- core_vif_vars 
  scenario_vif_vars <- generate_scenario_variable_list(core_vif_vars, current_scenario_name, config)
  if(length(scenario_vif_vars) < 1) stop("No VIF variables defined for current scenario.")
  env_predictor_stack_current <- load_selected_env_data(current_scenario_name, scenario_vif_vars, config)
}
if(is.null(env_predictor_stack_current)) stop("Failed to load 'current' environmental predictors.")
cat(paste(Sys.time(), "[INFO] Loaded 'current' env stack for SDMtune & BIOMOD2. Layers:", paste(names(env_predictor_stack_current), collapse=", "), "\n"))

# --- 6. Create Intermediate Output Dirs ---
sdmtune_intermediate_models_dir_path <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix_sdmtune))
sdmtune_intermediate_results_dir_path <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix_sdmtune))
dir.create(sdmtune_intermediate_models_dir_path, recursive = TRUE, showWarnings = FALSE)
dir.create(sdmtune_intermediate_results_dir_path, recursive = TRUE, showWarnings = FALSE)
dir.create(config$species_log_dir, recursive = TRUE, showWarnings = FALSE) 
cat(paste(Sys.time(), "[INFO] SDMtune intermediate dirs checked/created for:", group_name, predictor_type_suffix_sdmtune, "\n"))

biomod2_intermediate_base_dir <- file.path(config$sdm_output_dir_intermediate, "biomod2_outputs")
biomod2_folder_suffix <- paste0(predictor_type_suffix_sdmtune, "_biomod2") 
biomod2_species_group_dir_path <- file.path(biomod2_intermediate_base_dir, paste0(group_name, biomod2_folder_suffix))
dir.create(biomod2_species_group_dir_path, recursive = TRUE, showWarnings = FALSE)
cat(paste(Sys.time(), "[INFO] BIOMOD2 intermediate output base for this run type:", biomod2_species_group_dir_path, "\n"))

# --- 7. Load Species List & Select One Species for Tweaking ---
species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
if(nrow(species_df) == 0) stop("Species list is empty.")
species_row_to_process <- species_df[1, ] 
species_name_actual <- species_row_to_process$scientificName
species_name_sanitized_for_files <- gsub(" ", "_", species_name_actual) 
species_name_for_biomod_resp <- gsub(" ", ".", species_name_sanitized_for_files) 
species_aphia_id <- species_row_to_process$AphiaID
config$AphiaID_for_seed <- species_aphia_id 

cat(paste(Sys.time(), "[INFO] --- Processing Species:", species_name_actual, "(AphiaID:", species_aphia_id, ") ---\n"))

# --- 8. Define Combined SDMtune + BIOMOD2 Processing Function ---
process_species_combined_sdm_and_biomod2 <- function(
    sp_row, cfg, 
    env_pred_paths_or_list_sdmtune_arg,
    global_current_env_stack, 
    grp_name, pred_suffix_sdmtune, use_pca_flag, occ_dir_sp, 
    tuning_scenario_sdmtune = "current"
) {
  
  sp_name <- sp_row$scientificName
  sp_name_sanitized_files <- gsub(" ", "_", sp_name)
  sp_name_biomod <- gsub(" ", ".", sp_name_sanitized_files)
  sp_aphia <- sp_row$AphiaID
  
  log_prefix <- paste(Sys.time(), paste0("[", sp_name, "]"))
  cat(log_prefix, "INFO --- Starting COMBINED SDMtune & BIOMOD2 (Multi-Algo Ensemble) processing ---\n")
  
  sdmtune_model_file_path <- file.path(cfg$models_dir_intermediate, paste0(grp_name, pred_suffix_sdmtune), paste0("sdm_model_", sp_name_sanitized_files, pred_suffix_sdmtune, ".rds"))
  sdmtune_tuning_rds_path <- file.path(cfg$results_dir_intermediate, paste0(grp_name, pred_suffix_sdmtune), paste0("sdm_tuning_", sp_name_sanitized_files, pred_suffix_sdmtune, "_object.rds"))
  
  cat(log_prefix, "DEBUG Loading and cleaning occurrences...\n")
  cfg_occ_load <- cfg; cfg_occ_load$thinning_method <- NULL 
  occ_data_raw_list <- load_clean_individual_occ_coords(sp_aphia, occ_dir_sp, cfg_occ_load, logger=NULL, species_log_file=NULL)
  if (is.null(occ_data_raw_list) || is.null(occ_data_raw_list$coords) || occ_data_raw_list$count < cfg$min_occurrences_sdm) {
    cat(log_prefix, "ERROR Insufficient occurrences after basic cleaning:", occ_data_raw_list$count %||% 0, "\n"); return(list(status="skipped_occurrences_initial"))
  }
  occs_coords_cleaned <- occ_data_raw_list$coords 
  cat(log_prefix, "INFO Occurrences after basic cleaning:", nrow(occs_coords_cleaned), "\n")
  
  occs_coords_thinned_sac <- occs_coords_cleaned
  occurrence_count_final <- nrow(occs_coords_thinned_sac)
  if (cfg$apply_sac_thinning) {
    cat(log_prefix, "DEBUG Applying SAC thinning...\n")
    sac_result <- thin_occurrences_by_sac(occs_coords_cleaned, global_current_env_stack, cfg, logger=NULL, species_log_file=NULL)
    if (!is.null(sac_result) && !is.null(sac_result$coords_thinned)) {
      occs_coords_thinned_sac <- sac_result$coords_thinned
      occurrence_count_final <- sac_result$n_thinned
      cat(log_prefix, "INFO Occs after SAC thinning:", occurrence_count_final, "(Dist:", sac_result$thinning_distance_km %||% "NA", "km)\n")
      if (occurrence_count_final < cfg$min_occurrences_sdm) { cat(log_prefix, "ERROR Insufficient occs after SAC.\n"); return(list(status="skipped_occurrences_sac"))}
    } else { cat(log_prefix, "WARN SAC thinning failed or no change, using pre-SAC coords.\n") }
  } else { cat(log_prefix, "INFO SAC thinning disabled.\n") }
  
  final_occs_coords <- occs_coords_thinned_sac
  final_occs_sf <- sf::st_as_sf(as.data.frame(final_occs_coords), coords=c("longitude","latitude"), crs=cfg$occurrence_crs)
  final_occs_sf$AphiaID <- sp_aphia
  if(sf::st_crs(final_occs_sf) != terra::crs(global_current_env_stack)){ final_occs_sf <- sf::st_transform(final_occs_sf, crs=terra::crs(global_current_env_stack)) }
  
  # ===========================================================================
  # --- SDMtune Workflow ---
  # ===========================================================================
  cat(log_prefix, "INFO --- Starting SDMtune Workflow ---\n")
  sdmtune_model_obj <- NULL; sdmtune_tuning_obj <- NULL; sdmtune_best_hypers_obj <- NULL
  sdmtune_swd <- NULL; sdmtune_bg_pts <- NULL; sdmtune_spec_stack <- NULL
  
  if (!cfg$force_rerun$run_standard_sdms && file.exists(sdmtune_model_file_path)) {
    cat(log_prefix, "INFO [SDMtune] Loading existing model:", basename(sdmtune_model_file_path), "\n")
    sdmtune_model_obj <- readRDS(sdmtune_model_file_path)
    if (file.exists(sdmtune_tuning_rds_path)) sdmtune_tuning_obj <- readRDS(sdmtune_tuning_rds_path)
    if (!is.null(sdmtune_tuning_obj) && inherits(sdmtune_tuning_obj, "SDMtune") && !is.null(attr(sdmtune_tuning_obj, "best_hypers"))) {
      sdmtune_best_hypers_obj <- attr(sdmtune_tuning_obj, "best_hypers")
    }
    
    bg_ret_sdm <- generate_sdm_background_obis(final_occs_sf, global_current_env_stack, cfg, NULL, NULL, sp_aphia)
    if(!is.null(bg_ret_sdm)) { 
      sdmtune_bg_pts <- bg_ret_sdm$background_points
      sdmtune_spec_stack <- bg_ret_sdm$species_specific_stack
      if(!is.null(sdmtune_bg_pts) && !is.null(sdmtune_spec_stack)) {
        sdmtune_swd <- tryCatch(SDMtune::prepareSWD(species=sp_name, p=final_occs_coords, a=sdmtune_bg_pts, env=sdmtune_spec_stack, verbose=FALSE), error=function(e)NULL)
      }
    }
  } else {
    cat(log_prefix, "INFO [SDMtune] Running new tuning and training...\n")
    bg_ret_sdm <- generate_sdm_background_obis(final_occs_sf, global_current_env_stack, cfg, NULL, NULL, sp_aphia)
    if (is.null(bg_ret_sdm) || is.null(bg_ret_sdm$background_points) || is.null(bg_ret_sdm$species_specific_stack)) {
      cat(log_prefix, "ERROR [SDMtune] Background/stack generation failed.\n"); return(list(status="error_sdmtune_bg"))
    }
    sdmtune_bg_pts <- bg_ret_sdm$background_points
    sdmtune_spec_stack <- bg_ret_sdm$species_specific_stack
    
    sdmtune_swd <- tryCatch(SDMtune::prepareSWD(species=sp_name, p=final_occs_coords, a=sdmtune_bg_pts, env=sdmtune_spec_stack, verbose=FALSE), error=function(e){cat(log_prefix, "ERROR [SDMtune] SWD prep failed.\n"); NULL})
    if(is.null(sdmtune_swd)) return(list(status="error_sdmtune_swd"))
    
    sdmtune_spatial_folds <- create_spatial_cv_folds_simplified(sdmtune_swd, sdmtune_spec_stack, cfg, NULL, NULL)
    if(is.null(sdmtune_spatial_folds)) { cat(log_prefix, "ERROR [SDMtune] Spatial folds failed.\n"); return(list(status="error_sdmtune_folds"))}
    
    sdmtune_tuning_obj <- run_sdm_tuning_scv(final_occs_coords, sdmtune_spec_stack, sdmtune_bg_pts, cfg, NULL, sp_name, NULL)
    if(is.null(sdmtune_tuning_obj) || !inherits(sdmtune_tuning_obj, "SDMtune") || is.null(attr(sdmtune_tuning_obj, "best_hypers"))) { 
      cat(log_prefix, "ERROR [SDMtune] Tuning failed or did not produce best_hypers attribute.\n"); return(list(status="error_sdmtune_tuning"))
    }
    sdmtune_best_hypers_obj <- attr(sdmtune_tuning_obj, "best_hypers")
    save_tuning_results(sdmtune_tuning_obj, sp_name_sanitized_files, pred_suffix_sdmtune, cfg, NULL, NULL)
    
    sdmtune_model_obj <- train_final_sdm(final_occs_coords, sdmtune_spec_stack, sdmtune_bg_pts, sdmtune_best_hypers_obj, cfg, NULL, sp_name, NULL)
    if(is.null(sdmtune_model_obj)) { cat(log_prefix, "ERROR [SDMtune] Training failed.\n"); return(list(status="error_sdmtune_train"))}
    save_final_model(sdmtune_model_obj, sp_name_sanitized_files, pred_suffix_sdmtune, grp_name, cfg, NULL, NULL)
  }
  
  if (!is.null(sdmtune_model_obj) && !is.null(sdmtune_spec_stack)) {
    log_final_model_metrics(sdmtune_model_obj, sdmtune_swd, sdmtune_spec_stack, sdmtune_tuning_obj, sp_name_sanitized_files, grp_name, pred_suffix_sdmtune, cfg, NULL, NULL)
    if(!is.null(sdmtune_swd)) calculate_and_save_vi(sdmtune_model_obj, sdmtune_swd, sp_name_sanitized_files, grp_name, pred_suffix_sdmtune, cfg, NULL, NULL)
  }
  cat(log_prefix, "INFO --- Finished SDMtune Workflow ---\n")
  
  # ===========================================================================
  # --- BIOMOD2 Workflow (Multi-Algo Ensemble) ---
  # ===========================================================================
  cat(log_prefix, "INFO --- Starting BIOMOD2 Workflow (Multi-Algo Ensemble) ---\n")
  myBiomodEM <- NULL # Initialize Ensemble Model object
  
  if (is.null(sdmtune_bg_pts) || is.null(sdmtune_spec_stack)) {
    cat(log_prefix, "ERROR [BIOMOD2] Missing background points or species stack from SDMtune. Skipping BIOMOD2.\n")
  } else {
    myBiomodData <- format_data_for_biomod2(sp_name_biomod, final_occs_coords, sdmtune_bg_pts, sdmtune_spec_stack, NULL)
    if (is.null(myBiomodData)) { cat(log_prefix, "ERROR [BIOMOD2] Data formatting failed.\n") }
    else {
      cat(log_prefix, "INFO [BIOMOD2] Data formatted.\n")
      myBiomodCVTable <- create_biomod2_block_cv_table(myBiomodData, sdmtune_spec_stack, cfg, NULL)
      if (is.null(myBiomodCVTable)) cat(log_prefix, "WARN [BIOMOD2] blockCV table failed. BIOMOD_Modeling will use internal random CV.\n")
      
      all_algos_to_run <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'RF', 'MAXNET', 'XGBOOST', 'MARS')
      # all_algos_to_run <- c('RF', 'MAXNET')
      
      user_options_for_maxnet <- list()
      maxnet_args_for_biomod_run <- list() 
      
      if (!is.null(sdmtune_best_hypers_obj) && inherits(sdmtune_best_hypers_obj, "data.frame") && nrow(sdmtune_best_hypers_obj) == 1) {
        cat(log_prefix, "INFO [BIOMOD2] Preparing MAXNET parameters from SDMtune best_hypers.\n")
        if ("reg" %in% names(sdmtune_best_hypers_obj)) {
          maxnet_args_for_biomod_run$regmult <- as.numeric(sdmtune_best_hypers_obj$reg)
        }
        if ("fc" %in% names(sdmtune_best_hypers_obj)) {
          maxnet_args_for_biomod_run$classes <- as.character(sdmtune_best_hypers_obj$fc)
        }
        if (length(maxnet_args_for_biomod_run) > 0) {
          user_options_for_maxnet[['MAXNET.binary.maxnet.maxnet']] <- list('_allData_allRun' = maxnet_args_for_biomod_run)
          cat(log_prefix, "INFO [BIOMOD2] Tuned MAXNET parameters prepared:", paste(names(maxnet_args_for_biomod_run), unlist(maxnet_args_for_biomod_run), collapse=", "), "\n")
        } else {
          cat(log_prefix, "WARN [BIOMOD2] No usable hyperparameters (reg, fc) found for MAXNET.\n")
        }
      } else {
        cat(log_prefix, "WARN [BIOMOD2] SDMtune best_hypers not available for MAXNET. It will use 'bigboss' defaults.\n")
      }
      
      cat(log_prefix, "INFO [BIOMOD2] Models to attempt:", paste(all_algos_to_run, collapse=", "), "\n")
      cfg$group_name_for_biomod_output <- grp_name 
      
      biomod_options_obj <- biomod2::bm_ModelingOptions(
        data.type = 'binary',
        models = all_algos_to_run,
        strategy = 'user.defined',
        user.val = user_options_for_maxnet, 
        user.base = 'bigboss', 
        bm.format = myBiomodData 
      )
      
      myBiomodModelOut <- run_biomod2_models_with_blockcv(
        myBiomodData, biomod_options_obj, all_algos_to_run, myBiomodCVTable,
        sp_name_sanitized_files, pred_suffix_sdmtune, 
        grp_name, cfg, NULL
      )
      
      # TODO: Skip single algo projections for now
      # --- Individual Model Projections (Current Scenario) ---
      if (!is.null(myBiomodModelOut) && length(myBiomodModelOut@models.computed) > 0) {
        cat(log_prefix, "INFO [BIOMOD2] Individual modeling completed. Projecting current scenario for all computed models...\n")
        cat("INFO SKIPPING SINGLE PROJECTIONS")
        # for (algo_fullname_computed in myBiomodModelOut@models.computed) {
        #   algo_shortname_computed <- tail(strsplit(algo_fullname_computed, "_")[[1]], 1)
        #   if (algo_shortname_computed %in% all_algos_to_run) {
        #     cat(log_prefix, "DEBUG [BIOMOD2] Projecting individual model:", algo_shortname_computed, "\n")
        #     proj <- project_biomod2_models_current(myBiomodModelOut, sdmtune_spec_stack,
        #                                            sp_name_sanitized_files, pred_suffix_sdmtune,
        #                                            algo_shortname_computed, cfg, NULL)
        #     if(!is.null(proj)) {
        #       projection_suffix_indiv <- ""
        #       if (algo_shortname_computed == "MAXNET") {
        #         projection_suffix_indiv <- if (length(maxnet_args_for_biomod_run) > 0 && "MAXNET.binary.maxnet.maxnet" %in% names(user_options_for_maxnet)) {
        #           paste0(pred_suffix_sdmtune, "_biomod2_MAXNET_tuned")
        #         } else {
        #           paste0(pred_suffix_sdmtune, "_biomod2_MAXNET_default")
        #         }
        #       } else {
        #         projection_suffix_indiv <- paste0(pred_suffix_sdmtune, "_biomod2_", algo_shortname_computed, "_default")
        #       }
        #       save_biomod2_projection(proj, sp_name_sanitized_files, tuning_scenario_sdmtune, projection_suffix_indiv, cfg, NULL, NULL)
        #     } else {
        #       cat(log_prefix, "WARN [BIOMOD2] Projection failed for individual model", algo_shortname_computed, "\n")
        #     }
        #   }
        # }
        
        # --- Ensemble Modeling ---
        cat(log_prefix, "INFO [BIOMOD2] Building ensemble model...\n")
        
        
        
        # current_wd_b2 <- getwd()
        # setwd(config$species_biomod_run_dir)
        get_built_models(myBiomodModelOut , full.name = NULL, PA = NULL, run = NULL, algo = c('RF', 'MAXNET'))
        
        cat("INFO PRINTING ALL MODEL RESULTS\n")
        all_evaluations <- get_evaluations(myBiomodModelOut)
        print(all_evaluations)
        cat("INFO END PRINTING ALL MODEL RESULTS END\n")
        
        myBiomodEM <- tryCatch({
          biomod2::BIOMOD_EnsembleModeling(
            bm.mod = myBiomodModelOut,
            models.chosen = 'all', # Let metric.select and metric.select.thresh filter
            em.by = 'all', # Use full models for ensembling
            em.algo = c('EMmean'), # As per paper's implication of equal weighting
            metric.select = c('TSS', 'ROC'),
            metric.select.thresh = c(0.7, 0.8), # TSS > 0.7 AND ROC > 0.8
            metric.eval = c('TSS', 'ROC'), # Metrics to evaluate the ensemble model
            var.import = 0, # No VI for ensemble for now
            seed.val = cfg$AphiaID_for_seed %||% 4242 
          )
        }, error = function(e) {
          cat(log_prefix, "ERROR [BIOMOD2] BIOMOD_EnsembleModeling failed:", e$message, "\n")
          return(NULL)
        })
        
        cat("INFO PRINTING ENSEMBLE MODEL RESULTS\n")
        ensemble_evaluations <- get_evaluations(myBiomodModelOut)
        print(ensemble_evaluations)
        cat("INFO END PRINTING ENSEMBLE MODEL RESULTS END\n")
        
        if (!is.null(myBiomodEM) && length(myBiomodEM@em.models_kept) > 0) {
          cat(log_prefix, "INFO [BIOMOD2] Ensemble model built. Models kept:", paste(myBiomodEM@em.models_kept, collapse=", "), "\n")
          
          # --- Ensemble Model Projection ---
          cat(log_prefix, "INFO [BIOMOD2] Projecting ensemble model to current scenario...\n")
          myBiomodEMProj <- tryCatch({
            biomod2::BIOMOD_EnsembleForecasting(
              bm.em = myBiomodEM,
              proj.name = paste0("Ensemble_Current_", format(Sys.time(), "%Y%m%d%H%M%S")),
              new.env = global_current_env_stack, # Project to full current extent
              models.chosen = get_built_models(myBiomodEM), # Project only the successfully built ensemble(s)
              metric.binary = NULL, # Keep continuous
              compress = TRUE,
              output.format = ".tif"
            )
          }, error = function(e) {
            cat(log_prefix, "ERROR [BIOMOD2] BIOMOD_EnsembleForecasting failed:", e$message, "\n")
            return(NULL)
          })
          
          if (!is.null(myBiomodEMProj)) {
            ens_proj_raster_raw <- biomod2::get_predictions(myBiomodEMProj)
            if (!is.null(ens_proj_raster_raw) && inherits(ens_proj_raster_raw, "SpatRaster") && terra::nlyr(ens_proj_raster_raw) > 0) {
              ens_proj_raster <- ens_proj_raster_raw / 1000 # Scale to 0-1
              names(ens_proj_raster) <- paste0("suitability_EnsembleMean_filtered_current")
              
              ensemble_projection_suffix <- paste0(pred_suffix_sdmtune, "_biomod2_EnsembleMean_filtered")
              save_biomod2_projection(ens_proj_raster, sp_name_sanitized_files, tuning_scenario_sdmtune, ensemble_projection_suffix, cfg, NULL, NULL)
              
              # --- Log Ensemble Metrics ---
              ensemble_evals_df <- biomod2::get_evaluations(myBiomodEM, as.data.frame = TRUE)
              if (!is.null(ensemble_evals_df) && nrow(ensemble_evals_df) > 0) {
                cat(log_prefix, "INFO [BIOMOD2] Ensemble Model Evaluations:\n")
                print(ensemble_evals_df)
                # TODO: Could add a more structured logging for ensemble metrics here if needed
              } else {
                cat(log_prefix, "WARN [BIOMOD2] Could not retrieve evaluations for the ensemble model.\n")
              }
            } else {
              cat(log_prefix, "ERROR [BIOMOD2] get_predictions for ensemble returned invalid raster.\n")
            }
          }
        } else {
          cat(log_prefix, "WARN [BIOMOD2] Ensemble modeling resulted in no models being kept or failed.\n")
        }
      } else { 
        cat(log_prefix, "ERROR [BIOMOD2] Individual modeling failed or no models computed. Cannot proceed to ensembling.\n")
      }
    }
  }
  cat(log_prefix, "INFO --- Finished BIOMOD2 Workflow (Multi-Algo Ensemble) ---\n")
  # setwd(current_wd_b2)
  
  cat(log_prefix, "INFO [SDMtune] SDMtune predictions are not run by default in this combined script. Focus is on BIOMOD2 outputs.\n")
  
  cat(log_prefix, "INFO --- Combined processing finished ---\n")
  # Return the main BIOMOD_Modeling output and the EnsembleModeling output for inspection
  return(list(status = "success_combined_multi_algo_ensemble", 
              species = sp_name, 
              occ_count = occurrence_count_final, 
              message = "SDMtune and BIOMOD2 (Multi-Algo Ensemble) steps completed.",
              biomod_model_out = myBiomodModelOut, # Return individual models object
              biomod_ensemble_model_out = myBiomodEM # Return ensemble object
  ))
} # End process_species_combined_sdm_and_biomod2


# --- 9. Run for the Selected Species ---
config$use_parallel <- FALSE; config$num_cores <- 1
if(exists("future")) future::plan(future::sequential)

cat(paste(Sys.time(), "[INFO] --- Calling COMBINED SDMtune & BIOMOD2 (Multi-Algo Ensemble) processing for species:", species_name_actual, "---\n"))
combined_run_result <- process_species_combined_sdm_and_biomod2(
  sp_row = species_row_to_process,
  cfg = config,
  env_pred_paths_or_list_sdmtune_arg = predictor_paths_or_list_sdmtune, 
  global_current_env_stack = env_predictor_stack_current, 
  grp_name = group_name,
  pred_suffix_sdmtune = predictor_type_suffix_sdmtune, 
  use_pca_flag = use_pca,
  occ_dir_sp = occurrence_dir
)

# --- 10. Process Results (Simplified for Tweaking) ---
cat(paste(Sys.time(), "[INFO]--- Combined SDMtune & BIOMOD2 (Multi-Algo Ensemble) Tweaking Run Result ---\n"))
if (!is.null(combined_run_result)) {
  cat(paste(Sys.time(), "[STATUS:", combined_run_result$status, "]", combined_run_result$species, "-", combined_run_result$message %||% "Completed.", "\n"))
  if (!is.null(combined_run_result$biomod_model_out)) {
    cat("  Individual models computed:", paste(combined_run_result$biomod_model_out@models.computed, collapse=", "),"\n")
  }
  if (!is.null(combined_run_result$biomod_ensemble_model_out)) {
    cat("  Ensemble models built:", paste(get_built_models(combined_run_result$biomod_ensemble_model_out), collapse=", "),"\n")
    cat("  Ensemble evaluations:\n")
    print(get_evaluations(combined_run_result$biomod_ensemble_model_out))
  }
} else {
  cat(paste(Sys.time(), "[ERROR] Combined processing returned NULL overall.\n"))
}

gc(full=TRUE)
cat(paste(Sys.time(), "[INFO]--- Script sdm_combined_tweaking_run.R (Multi-Algo Ensemble) finished. ---\n"))
#-------------------------------------------------------------------------------