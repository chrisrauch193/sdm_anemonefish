# scripts/sdm_runs/experiments/sdmtune_all_algos_biomod2_combined_run.R
#-------------------------------------------------------------------------------
# Combined Tweaking Script: 
# 1. Runs SDMtune for Maxnet, RF, ANN, BRT.
# 2. Runs BIOMOD2 with multiple algorithms:
#    - SDMtune-tuned algorithms use their best hyperparameters.
#    - Other BIOMOD2 algorithms use 'bigboss' defaults.
# 3. Builds an ensemble model from best-performing BIOMOD2 individual models.
# For a single species and current scenario.
#-------------------------------------------------------------------------------

cat("--- Running Script sdmtune_all_algos_biomod2_combined_run.R ---\n")

# --- 1. Setup: Load Config FIRST ---
# ... (same as before) ...
if (file.exists("scripts/config.R")) {
  source("scripts/config.R")
  if (!exists("config") || !is.list(config)) {
    stop("FATAL: 'config' list object not found or invalid after sourcing scripts/config.R")
  }
  # Ensure new hypergrids are loaded from config
  if(is.null(config$sdm_tune_grids) || !is.list(config$sdm_tune_grids)) {
    stop("FATAL: 'config$sdm_tune_grids' not found or invalid. Ensure it's defined in config.R")
  }
} else {
  stop("FATAL: Configuration file 'scripts/config.R' not found.")
}


# --- 2. Load Required Packages ---
# ... (same as before) ...
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr, 
               biomod2, maxnet, blockCV, 
               nnet, gbm, randomForest, # Dependencies for SDMtune algos
               presenceabsence, mda, gam, earth, xgboost # Dependencies for other BIOMOD2 algos
)

# --- 3. Source Helper Functions ---
# ... (same as before) ...
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))

# --- 4. Define Group Specifics & Predictor Type ---
# ... (same as before) ...
cat(paste(Sys.time(), "[INFO] --- Starting Script sdmtune_all_algos_biomod2_combined_run.R ---\n"))

group_name <- "anemone" 
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors 
predictor_type_suffix_sdmtune <- ifelse(use_pca, "_pca", "_vif") 

cat(paste(Sys.time(), "[INFO] --- Processing Group:", group_name, "---\n"))
cat(paste(Sys.time(), "[INFO] --- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---\n"))


# --- 5. Load Predictor Information (for the 'current' scenario) ---
# ... (same as before) ...
current_scenario_name <- "current"
env_predictor_stack_current <- NULL 
predictor_paths_or_list_sdmtune_files <- NULL 

if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (!file.exists(pca_paths_rds)) stop("PCA raster paths RDS file not found.")
  predictor_paths_or_list_sdmtune_files <- readRDS(pca_paths_rds) 
  current_pca_path <- predictor_paths_or_list_sdmtune_files[[current_scenario_name]]
  if (is.null(current_pca_path) || !file.exists(current_pca_path)) stop("PCA stack for 'current' scenario not found.")
  env_predictor_stack_current <- tryCatch(terra::rast(current_pca_path), error = function(e) {
    cat(paste(Sys.time(), "[ERROR] Failed to load PCA stack:", e$message, "\n")); NULL
  })
} else {
  core_vif_vars <- config$final_vars_vif_anemone 
  predictor_paths_or_list_sdmtune_files <- core_vif_vars 
  scenario_vif_vars <- generate_scenario_variable_list(core_vif_vars, current_scenario_name, config)
  if(length(scenario_vif_vars) < 1) stop("No VIF variables defined for current scenario.")
  env_predictor_stack_current <- load_selected_env_data(current_scenario_name, scenario_vif_vars, config)
}
if(is.null(env_predictor_stack_current)) stop("Failed to load 'current' environmental predictors.")
cat(paste(Sys.time(), "[INFO] Loaded 'current' env stack for SDMtune & BIOMOD2. Layers:", paste(names(env_predictor_stack_current), collapse=", "), "\n"))


# --- 6. Create Intermediate Output Dirs ---
# ... (same as before for SDMtune and BIOMOD2 main group folder) ...
# SDMtune intermediate dirs for EACH algorithm type will be handled within the loop
sdmtune_base_intermediate_models_dir <- config$models_dir_intermediate
sdmtune_base_intermediate_results_dir <- config$results_dir_intermediate
dir.create(config$species_log_dir, recursive = TRUE, showWarnings = FALSE) 

biomod2_intermediate_base_dir <- file.path(config$sdm_output_dir_intermediate, "biomod2_outputs")
biomod2_folder_suffix <- paste0(predictor_type_suffix_sdmtune, "_biomod2_multi_algo_ensemble") # New suffix for this run type
biomod2_species_group_dir_path <- file.path(biomod2_intermediate_base_dir, paste0(group_name, biomod2_folder_suffix))
dir.create(biomod2_species_group_dir_path, recursive = TRUE, showWarnings = FALSE)
cat(paste(Sys.time(), "[INFO] BIOMOD2 intermediate output base for this run type:", biomod2_species_group_dir_path, "\n"))


# --- 7. Load Species List & Select One Species for Tweaking ---
# ... (same as before) ...
species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
if(nrow(species_df) == 0) stop("Species list is empty.")
species_row_to_process <- species_df[1, ] 
species_name_actual <- species_row_to_process$scientificName
species_name_sanitized_for_files <- gsub(" ", "_", species_name_actual) 
species_name_for_biomod_resp <- gsub(" ", ".", species_name_sanitized_for_files) 
species_aphia_id <- species_row_to_process$AphiaID
config$AphiaID_for_seed <- species_aphia_id 

cat(paste(Sys.time(), "[INFO] --- Processing Species:", species_name_actual, "(AphiaID:", species_aphia_id, ") ---\n"))


# --- 8. Define Combined SDMtune (All Algos) + BIOMOD2 (Multi-Algo Ensemble) Processing Function ---
process_species_all_sdms_combined <- function(
    sp_row, cfg_obj, # Pass the full config object
    env_pred_paths_or_list_sdmtune_files_arg, 
    global_current_env_stack, 
    grp_name, pred_suffix_sdmtune_base, # e.g. _pca or _vif
    use_pca_flag, occ_dir_sp, 
    tuning_scenario_sdmtune = "current"
) {
  
  sp_name <- sp_row$scientificName
  sp_name_sanitized_files <- gsub(" ", "_", sp_name)
  sp_name_biomod <- gsub(" ", ".", sp_name_sanitized_files)
  sp_aphia <- sp_row$AphiaID
  
  log_prefix <- paste(Sys.time(), paste0("[", sp_name, "]"))
  cat(log_prefix, "INFO --- Starting COMBINED SDMtune (All Algos) & BIOMOD2 (Multi-Algo Ensemble) ---\n")
  
  # --- Load & Prepare Occurrences (Common for all) ---
  # ... (same occurrence loading and SAC thinning as before) ...
  cat(log_prefix, "DEBUG Loading and cleaning occurrences...\n")
  cfg_occ_load <- cfg_obj; cfg_occ_load$thinning_method <- NULL 
  occ_data_raw_list <- load_clean_individual_occ_coords(sp_aphia, occ_dir_sp, cfg_occ_load, logger=NULL, species_log_file=NULL)
  if (is.null(occ_data_raw_list) || is.null(occ_data_raw_list$coords) || occ_data_raw_list$count < cfg_obj$min_occurrences_sdm) {
    cat(log_prefix, "ERROR Insufficient occurrences after basic cleaning:", occ_data_raw_list$count %||% 0, "\n"); return(list(status="skipped_occurrences_initial"))
  }
  occs_coords_cleaned <- occ_data_raw_list$coords 
  cat(log_prefix, "INFO Occurrences after basic cleaning:", nrow(occs_coords_cleaned), "\n")
  
  occs_coords_thinned_sac <- occs_coords_cleaned
  occurrence_count_final <- nrow(occs_coords_thinned_sac)
  if (cfg_obj$apply_sac_thinning) {
    cat(log_prefix, "DEBUG Applying SAC thinning...\n")
    sac_result <- thin_occurrences_by_sac(occs_coords_cleaned, global_current_env_stack, cfg_obj, logger=NULL, species_log_file=NULL)
    if (!is.null(sac_result) && !is.null(sac_result$coords_thinned)) {
      occs_coords_thinned_sac <- sac_result$coords_thinned
      occurrence_count_final <- sac_result$n_thinned
      cat(log_prefix, "INFO Occs after SAC thinning:", occurrence_count_final, "(Dist:", sac_result$thinning_distance_km %||% "NA", "km)\n")
      if (occurrence_count_final < cfg_obj$min_occurrences_sdm) { cat(log_prefix, "ERROR Insufficient occs after SAC.\n"); return(list(status="skipped_occurrences_sac"))}
    } else { cat(log_prefix, "WARN SAC thinning failed or no change, using pre-SAC coords.\n") }
  } else { cat(log_prefix, "INFO SAC thinning disabled.\n") }
  
  final_occs_coords <- occs_coords_thinned_sac
  final_occs_sf <- sf::st_as_sf(as.data.frame(final_occs_coords), coords=c("longitude","latitude"), crs=cfg_obj$occurrence_crs)
  final_occs_sf$AphiaID <- sp_aphia
  if(sf::st_crs(final_occs_sf) != terra::crs(global_current_env_stack)){ final_occs_sf <- sf::st_transform(final_occs_sf, crs=terra::crs(global_current_env_stack)) }
  
  # --- Background points and species-specific stack (common for all SDMtune runs) ---
  bg_ret_sdm_common <- generate_sdm_background_obis(final_occs_sf, global_current_env_stack, cfg_obj, NULL, NULL, sp_aphia)
  if (is.null(bg_ret_sdm_common) || is.null(bg_ret_sdm_common$background_points) || is.null(bg_ret_sdm_common$species_specific_stack)) {
    cat(log_prefix, "ERROR [SDMtune_Setup] Common Background/stack generation failed.\n"); return(list(status="error_sdmtune_common_bg"))
  }
  sdmtune_bg_pts_common <- bg_ret_sdm_common$background_points
  sdmtune_spec_stack_common <- bg_ret_sdm_common$species_specific_stack # This is the stack masked to the BG extent
  
  # ===========================================================================
  # --- SDMtune Workflow (Iterate for ANN, BRT, Maxnet, RF) ---
  # ===========================================================================
  cat(log_prefix, "INFO --- Starting SDMtune Workflow for Multiple Algorithms ---\n")
  
  sdmtune_algos_to_run <- c("Maxnet", "RF", "ANN", "BRT")
  all_sdmtune_best_hypers <- list() # Store best hypers for each algo
  
  for (sdm_method_current in sdmtune_algos_to_run) {
    cat(log_prefix, "INFO [SDMtune] --- Processing Algorithm:", sdm_method_current, "---\n")
    
    current_algo_pred_suffix <- paste0(pred_suffix_sdmtune_base, "_sdmtune_", tolower(sdm_method_current))
    
    # Create specific intermediate dirs for this SDMtune algo run
    sdmtune_algo_models_dir <- file.path(cfg_obj$models_dir_intermediate, paste0(grp_name, current_algo_pred_suffix))
    sdmtune_algo_results_dir <- file.path(cfg_obj$results_dir_intermediate, paste0(grp_name, current_algo_pred_suffix))
    dir.create(sdmtune_algo_models_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(sdmtune_algo_results_dir, recursive = TRUE, showWarnings = FALSE)
    
    sdmtune_model_file_path_algo <- file.path(sdmtune_algo_models_dir, paste0("sdm_model_", sp_name_sanitized_files, current_algo_pred_suffix, ".rds"))
    sdmtune_tuning_rds_path_algo <- file.path(sdmtune_algo_results_dir, paste0("sdm_tuning_", sp_name_sanitized_files, current_algo_pred_suffix, "_object.rds"))
    
    # Prepare a config copy for this specific SDMtune algo
    cfg_algo_sdmtune <- cfg_obj 
    cfg_algo_sdmtune$sdm_method <- sdm_method_current
    cfg_algo_sdmtune$sdm_tune_grid <- cfg_obj$sdm_tune_grids[[sdm_method_current]]
    if(is.null(cfg_algo_sdmtune$sdm_tune_grid)) {
      cat(log_prefix, "WARN [SDMtune] No specific tune grid found for", sdm_method_current, "in config$sdm_tune_grids. Skipping tuning for this algo.\n")
      all_sdmtune_best_hypers[[sdm_method_current]] <- NULL # Mark as not tuned
      next # Skip to next SDMtune algo
    }
    
    sdmtune_swd_algo <- tryCatch(SDMtune::prepareSWD(species=sp_name, p=final_occs_coords, a=sdmtune_bg_pts_common, env=sdmtune_spec_stack_common, verbose=FALSE), 
                                 error=function(e){cat(log_prefix, "ERROR [SDMtune:", sdm_method_current, "] SWD prep failed.\n"); NULL})
    if(is.null(sdmtune_swd_algo)) {
      all_sdmtune_best_hypers[[sdm_method_current]] <- NULL; next
    }
    
    # Run tuning and training (using your existing helpers, with cfg_algo_sdmtune)
    # We assume run_standard_sdms flag applies to all these SDMtune runs
    sdmtune_tuning_obj_algo <- NULL
    sdmtune_best_hypers_algo <- NULL
    
    if (!cfg_obj$force_rerun$run_standard_sdms && file.exists(sdmtune_model_file_path_algo)) {
      cat(log_prefix, "INFO [SDMtune:", sdm_method_current, "] Loading existing model and hypers.\n")
      # Try to load existing hypers if model exists
      if (file.exists(sdmtune_tuning_rds_path_algo)) {
        loaded_tuning_obj <- readRDS(sdmtune_tuning_rds_path_algo)
        if (!is.null(loaded_tuning_obj) && inherits(loaded_tuning_obj, "SDMtune") && !is.null(attr(loaded_tuning_obj, "best_hypers"))){
          sdmtune_best_hypers_algo <- attr(loaded_tuning_obj, "best_hypers")
        }
      }
    } else {
      cat(log_prefix, "INFO [SDMtune:", sdm_method_current, "] Running new tuning.\n")
      sdmtune_spatial_folds_algo <- create_spatial_cv_folds_simplified(sdmtune_swd_algo, sdmtune_spec_stack_common, cfg_algo_sdmtune, NULL, NULL)
      if(is.null(sdmtune_spatial_folds_algo)) { 
        cat(log_prefix, "ERROR [SDMtune:", sdm_method_current, "] Spatial folds failed.\n"); 
        all_sdmtune_best_hypers[[sdm_method_current]] <- NULL; next
      }
      
      sdmtune_tuning_obj_algo <- run_sdm_tuning_scv(final_occs_coords, sdmtune_spec_stack_common, sdmtune_bg_pts_common, cfg_algo_sdmtune, NULL, sp_name, NULL)
      if(is.null(sdmtune_tuning_obj_algo) || !inherits(sdmtune_tuning_obj_algo, "SDMtune") || is.null(attr(sdmtune_tuning_obj_algo, "best_hypers"))) { 
        cat(log_prefix, "ERROR [SDMtune:", sdm_method_current, "] Tuning failed or no best hypers.\n")
        all_sdmtune_best_hypers[[sdm_method_current]] <- NULL; next
      }
      sdmtune_best_hypers_algo <- attr(sdmtune_tuning_obj_algo, "best_hypers")
      save_tuning_results(sdmtune_tuning_obj_algo, sp_name_sanitized_files, current_algo_pred_suffix, cfg_obj, NULL, NULL) # Save with algo-specific suffix
      
      # Train and save the SDMtune model (optional, but good for direct comparison)
      # temp_sdmtune_model_obj <- train_final_sdm(final_occs_coords, sdmtune_spec_stack_common, sdmtune_bg_pts_common, sdmtune_best_hypers_algo, cfg_algo_sdmtune, NULL, sp_name, NULL)
      # if(!is.null(temp_sdmtune_model_obj)) {
      #    save_final_model(temp_sdmtune_model_obj, sp_name_sanitized_files, current_algo_pred_suffix, grp_name, cfg_obj, NULL, NULL)
      # }
    }
    all_sdmtune_best_hypers[[sdm_method_current]] <- sdmtune_best_hypers_algo
    cat(log_prefix, "INFO [SDMtune:", sdm_method_current, "] Processing finished. Best hypers:", if(!is.null(sdmtune_best_hypers_algo)) paste(names(sdmtune_best_hypers_algo), sdmtune_best_hypers_algo, collapse="; ") else "NONE", "\n")
  } # End SDMtune algo loop
  cat(log_prefix, "INFO --- Finished All SDMtune Algorithm Workflows ---\n")
  
  # ===========================================================================
  # --- BIOMOD2 Workflow (Multi-Algo Ensemble, using SDMtune hypers where available) ---
  # ===========================================================================
  cat(log_prefix, "INFO --- Starting BIOMOD2 Workflow (Multi-Algo Ensemble from SDMtune Hypers) ---\n")
  myBiomodEM <- NULL 
  
  myBiomodData <- format_data_for_biomod2(sp_name_biomod, final_occs_coords, sdmtune_bg_pts_common, sdmtune_spec_stack_common, NULL)
  if (is.null(myBiomodData)) { cat(log_prefix, "ERROR [BIOMOD2] Data formatting failed.\n"); return(list(status="error_biomod_format")) }
  
  cat(log_prefix, "INFO [BIOMOD2] Data formatted.\n")
  myBiomodCVTable <- create_biomod2_block_cv_table(myBiomodData, sdmtune_spec_stack_common, cfg_obj, NULL)
  if (is.null(myBiomodCVTable)) cat(log_prefix, "WARN [BIOMOD2] blockCV table failed. BIOMOD_Modeling will use internal random CV.\n")
  
  # Algorithms for BIOMOD2 (from paper, ensuring names match BIOMOD2)
  # biomod2_algos_to_run <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'RF', 'MAXNET', 'MARS') # XGBOOST often needs separate install of xgboost R package
  biomod2_algos_to_run <- c('GBM', 'ANN', 'RF', 'MAXNET')
  # biomod2_algos_to_run <- c('RF', 'MAXNET') # For quicker testing
  
  # Prepare user.val list for BIOMOD2 based on successful SDMtune runs
  user_val_biomod2 <- list()
  
  # MAXNET
  if (!is.null(all_sdmtune_best_hypers[["Maxnet"]])) {
    prep_maxnet <- prepare_biomod2_maxnet_user_val(all_sdmtune_best_hypers[["Maxnet"]])
    if (!is.null(prep_maxnet)) user_val_biomod2[['MAXNET.binary.maxnet.maxnet']] <- prep_maxnet
  }
  # RF
  if (!is.null(all_sdmtune_best_hypers[["RF"]])) {
    prep_rf <- prepare_biomod2_rf_user_val(all_sdmtune_best_hypers[["RF"]])
    if (!is.null(prep_rf)) user_val_biomod2[['RF.binary.randomForest.randomForest']] <- prep_rf
  }
  # ANN
  if (!is.null(all_sdmtune_best_hypers[["ANN"]])) {
    prep_ann <- prepare_biomod2_ann_user_val(all_sdmtune_best_hypers[["ANN"]])
    if (!is.null(prep_ann)) user_val_biomod2[['ANN.binary.nnet.nnet']] <- prep_ann
  }
  # BRT for GBM
  if (!is.null(all_sdmtune_best_hypers[["BRT"]])) {
    prep_gbm <- prepare_biomod2_gbm_user_val(all_sdmtune_best_hypers[["BRT"]])
    if (!is.null(prep_gbm)) user_val_biomod2[['GBM.binary.gbm.gbm']] <- prep_gbm
  }
  
  cat(log_prefix, "INFO [BIOMOD2] User-defined options prepared for:", paste(names(user_val_biomod2), collapse=", "), "\n")
  cat(log_prefix, "INFO [BIOMOD2] BIOMOD2 Models to attempt:", paste(biomod2_algos_to_run, collapse=", "), "\n")
  
  cfg_obj$group_name_for_biomod_output <- grp_name 
  
  biomod_options_obj <- biomod2::bm_ModelingOptions(
    data.type = 'binary',
    models = biomod2_algos_to_run,
    strategy = 'user.defined', 
    user.val = user_val_biomod2, 
    user.base = 'bigboss', 
    bm.format = myBiomodData 
  )
  
  myBiomodModelOut <- run_biomod2_models_with_blockcv(
    myBiomodData, biomod_options_obj, biomod2_algos_to_run, myBiomodCVTable,
    sp_name_sanitized_files, pred_suffix_sdmtune_base, # Use base suffix for folder structure
    grp_name, cfg_obj, NULL
  )
  
  # --- Ensemble Modeling ---
  if (!is.null(myBiomodModelOut) && length(myBiomodModelOut@models.computed) > 0) {
    cat(log_prefix, "INFO [BIOMOD2] Individual modeling completed. Building ensemble model...\n")
    
    cat(log_prefix, "INFO [BIOMOD2] All Individual Model Evaluations from BIOMOD_Modeling:\n")
    all_individual_evals <- get_evaluations(myBiomodModelOut, as.data.frame = TRUE)
    print(all_individual_evals)
    
    myBiomodEM <- tryCatch({
      biomod2::BIOMOD_EnsembleModeling(
        bm.mod = myBiomodModelOut,
        models.chosen = 'all', 
        em.by = 'all', 
        em.algo = c('EMmean'), 
        metric.select = c('TSS', 'ROC'), # Use both, BIOMOD2 will apply AND logic
        metric.select.thresh = c(0.7, 0.8), 
        metric.eval = c('TSS', 'ROC'), 
        var.import = 0, 
        seed.val = cfg_obj$AphiaID_for_seed %||% 4242 
      )
    }, error = function(e) {
      cat(log_prefix, "ERROR [BIOMOD2] BIOMOD_EnsembleModeling failed:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(myBiomodEM) && length(myBiomodEM@em.models_kept) > 0) {
      cat(log_prefix, "INFO [BIOMOD2] Ensemble model built. Models kept:", paste(myBiomodEM@em.models_kept, collapse=", "), "\n")
      
      cat(log_prefix, "INFO [BIOMOD2] Ensemble Model Evaluations:\n")
      ensemble_evals_df <- biomod2::get_evaluations(myBiomodEM, as.data.frame = TRUE)
      print(ensemble_evals_df)
      
      # --- Ensemble Model Projection ---
      cat(log_prefix, "INFO [BIOMOD2] Projecting ensemble model to current scenario...\n")
      myBiomodEMProj <- tryCatch({
        biomod2::BIOMOD_EnsembleForecasting(
          bm.em = myBiomodEM,
          proj.name = paste0("Ensemble_Current_", format(Sys.time(), "%Y%m%d%H%M%S")),
          new.env = global_current_env_stack, 
          models.chosen = get_built_models(myBiomodEM), 
          metric.binary = NULL, 
          compress = TRUE,
          output.format = ".tif"
        )
      }, error = function(e) {
        cat(log_prefix, "ERROR [BIOMOD2] BIOMOD_EnsembleForecasting failed:", e$message, "\n"); return(NULL)
      })
      
      if (!is.null(myBiomodEMProj)) {
        ens_proj_raster_raw <- biomod2::get_predictions(myBiomodEMProj)
        if (!is.null(ens_proj_raster_raw) && inherits(ens_proj_raster_raw, "SpatRaster") && terra::nlyr(ens_proj_raster_raw) > 0) {
          ens_proj_raster <- ens_proj_raster_raw / 1000 
          names(ens_proj_raster) <- paste0("suitability_EnsembleMean_filtered_current")
          
          ensemble_projection_suffix <- paste0(pred_suffix_sdmtune_base, "_biomod2_EnsembleMean_filtered_AllSDMtune")
          save_biomod2_projection(ens_proj_raster, sp_name_sanitized_files, tuning_scenario_sdmtune, ensemble_projection_suffix, cfg_obj, NULL, NULL)
        } else { cat(log_prefix, "ERROR [BIOMOD2] get_predictions for ensemble returned invalid raster.\n") }
      }
    } else { cat(log_prefix, "WARN [BIOMOD2] Ensemble modeling resulted in no models being kept or failed.\n") }
  } else { cat(log_prefix, "ERROR [BIOMOD2] Individual modeling failed or no models computed. Cannot proceed to ensembling.\n") }
  
  cat(log_prefix, "INFO --- Finished BIOMOD2 Workflow (Multi-Algo Ensemble from SDMtune Hypers) ---\n")
  
  cat(log_prefix, "INFO --- Combined processing finished ---\n")
  return(list(status = "success_sdmtune_all_biomod2_ensemble", 
              species = sp_name, 
              occ_count = occurrence_count_final, 
              message = "SDMtune (All Algos) and BIOMOD2 (Multi-Algo Ensemble) steps completed.",
              sdmtune_all_hypers = all_sdmtune_best_hypers,
              biomod_individual_models_out = myBiomodModelOut,
              biomod_ensemble_model_out = myBiomodEM 
  ))
} # End process_species_all_sdms_combined


# --- 9. Run for the Selected Species ---
config$use_parallel <- FALSE; config$num_cores <- 1
if(exists("future")) future::plan(future::sequential)

cat(paste(Sys.time(), "[INFO] --- Calling COMBINED SDMtune (All Algos) & BIOMOD2 (Multi-Algo Ensemble) processing for species:", species_name_actual, "---\n"))
combined_run_result <- process_species_all_sdms_combined(
  sp_row = species_row_to_process,
  cfg_obj = config, # Pass the full config object
  env_pred_paths_or_list_sdmtune_files_arg = predictor_paths_or_list_sdmtune_files, 
  global_current_env_stack = env_predictor_stack_current, 
  grp_name = group_name,
  pred_suffix_sdmtune_base = predictor_type_suffix_sdmtune, 
  use_pca_flag = use_pca,
  occ_dir_sp = occurrence_dir
)

# --- 10. Process Results ---
cat(paste(Sys.time(), "[INFO]--- Combined SDMtune (All Algos) & BIOMOD2 (Multi-Algo Ensemble) Tweaking Run Result ---\n"))
if (!is.null(combined_run_result)) {
  cat(paste(Sys.time(), "[STATUS:", combined_run_result$status, "]", combined_run_result$species, "-", combined_run_result$message %||% "Completed.", "\n"))
  cat("  SDMtune Hyperparameters obtained:\n")
  print(combined_run_result$sdmtune_all_hypers)
  if (!is.null(combined_run_result$biomod_individual_models_out)) {
    cat("  BIOMOD2 Individual models computed:", paste(combined_run_result$biomod_individual_models_out@models.computed, collapse=", "),"\n")
  }
  if (!is.null(combined_run_result$biomod_ensemble_model_out)) {
    cat("  BIOMOD2 Ensemble models built:", paste(get_built_models(combined_run_result$biomod_ensemble_model_out), collapse=", "),"\n")
    cat("  BIOMOD2 Ensemble evaluations:\n")
    print(get_evaluations(combined_run_result$biomod_ensemble_model_out))
  }
} else {
  cat(paste(Sys.time(), "[ERROR] Combined processing returned NULL overall.\n"))
}

gc(full=TRUE)
cat(paste(Sys.time(), "[INFO]--- Script sdmtune_all_algos_biomod2_combined_run.R finished. ---\n"))
#-------------------------------------------------------------------------------