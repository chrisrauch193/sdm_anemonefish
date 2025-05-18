# scripts/sdm_runs/experiments/biomod2_tweaking_base_run.R
#-------------------------------------------------------------------------------
# Tweaking script for a basic BIOMOD2 run with default/bigboss parameters.
# Focuses on running MAXNET for a single species.
# Placeholder for blockCV.
#-------------------------------------------------------------------------------

cat("--- Running Script biomod2_tweaking_base_run.R ---\n")

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
pacman::p_load(terra, sf, dplyr, readr, tools, stringr,
               biomod2, maxnet, # Core biomod2 and maxnet for MAXNET
               presenceabsence, randomForest, gbm, mda, gam, earth, xgboost, # Common biomod2 algo deps
               blockCV # For when we add it
)

# --- 3. Source Helper Functions ---
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))

# --- 4. Define Group Specifics & Predictor Type ---
cat(paste(Sys.time(), "[INFO] --- Starting Script biomod2_tweaking_base_run.R ---\n"))

group_name <- "anemone" 
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors 
predictor_type_suffix_base <- ifelse(use_pca, "_pca", "_vif") # Base for predictors

cat(paste(Sys.time(), "[INFO] --- Processing Group:", group_name, "---\n"))
cat(paste(Sys.time(), "[INFO] --- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---\n"))

# --- 5. Load Predictor Information (for the 'current' scenario) ---
current_scenario_name <- "current"
env_predictor_stack_current <- NULL

if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (!file.exists(pca_paths_rds)) stop("PCA raster paths RDS file not found.")
  predictor_paths_list <- readRDS(pca_paths_rds)
  current_pca_path <- predictor_paths_list[[current_scenario_name]]
  if (is.null(current_pca_path) || !file.exists(current_pca_path)) stop("PCA stack for 'current' scenario not found.")
  env_predictor_stack_current <- tryCatch(terra::rast(current_pca_path), error = function(e) {
    cat(paste(Sys.time(), "[ERROR] Failed to load PCA stack:", e$message, "\n")); NULL
  })
} else {
  core_vif_vars <- config$final_vars_vif_anemone 
  scenario_vif_vars <- generate_scenario_variable_list(core_vif_vars, current_scenario_name, config)
  if(length(scenario_vif_vars) < 1) stop("No VIF variables defined for current scenario.")
  env_predictor_stack_current <- load_selected_env_data(current_scenario_name, scenario_vif_vars, config)
}
if(is.null(env_predictor_stack_current)) stop("Failed to load 'current' environmental predictors.")
cat(paste(Sys.time(), "[INFO] Loaded 'current' env stack. Layers:", paste(names(env_predictor_stack_current), collapse=", "), "\n"))

# --- 6. Create BIOMOD2 Specific Intermediate Output Dirs ---
biomod2_intermediate_base_dir <- file.path(config$sdm_output_dir_intermediate, "biomod2_outputs")
# Suffix for BIOMOD2 output folders will be like _pca_biomod2
biomod2_folder_suffix <- paste0(predictor_type_suffix_base, "_biomod2") 
biomod2_species_group_dir <- file.path(biomod2_intermediate_base_dir, paste0(group_name, biomod2_folder_suffix))
dir.create(biomod2_species_group_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste(Sys.time(), "[INFO] BIOMOD2 intermediate output base for this run type:", biomod2_species_group_dir, "\n"))

# --- 7. Load Species List & Select One Species for Tweaking ---
species_df <- readr::read_csv(species_list_file, show_col_types = FALSE)
if(nrow(species_df) == 0) stop("Species list is empty.")
species_row_to_process <- species_df[1, ]
species_name_actual <- species_row_to_process$scientificName
species_name_sanitized_for_files <- gsub(" ", "_", species_name_actual) 
species_name_for_biomod_resp <- gsub(" ", ".", species_name_sanitized_for_files) # BIOMOD_FormatingData needs this
species_aphia_id <- species_row_to_process$AphiaID
config$AphiaID_for_seed <- species_aphia_id # Store in config for helper seed access

cat(paste(Sys.time(), "[INFO] --- Processing Species:", species_name_actual, "(AphiaID:", species_aphia_id, ") ---\n"))

# --- 8. Define and Run process_species_biomod2_basic Function ---
process_species_biomod2_basic <- function(sp_name_actual, sp_name_biomod, sp_aphia, 
                                          occ_dir_sp, env_stack_curr, 
                                          pred_type_suffix_base_val, # e.g. _pca
                                          grp_name_val, cfg_obj) {
  
  log_prefix <- paste(Sys.time(), paste0("[", sp_name_actual, "]"))
  cat(log_prefix, "INFO [BIOMOD2_Process] --- Starting basic BIOMOD2 processing ---\n")
  
  # 8.1 Load Occurrences
  cat(log_prefix, "DEBUG [BIOMOD2_Process] Loading occurrences...\n")
  cfg_occ <- cfg_obj; cfg_occ$thinning_method <- NULL 
  occ_data_res <- load_clean_individual_occ_coords(sp_aphia, occ_dir_sp, cfg_occ, logger=NULL, species_log_file=NULL)
  if (is.null(occ_data_res) || is.null(occ_data_res$coords) || occ_data_res$count < cfg_obj$min_occurrences_sdm) {
    cat(log_prefix, "ERROR [BIOMOD2_Process] Insufficient occurrences:", occ_data_res$count %||% 0, "\n"); return(NULL)
  }
  pres_coords_current <- occ_data_res$coords
  cat(log_prefix, "INFO [BIOMOD2_Process] Presences:", nrow(pres_coords_current), "\n")
  
  # 8.2 Background Points
  cat(log_prefix, "DEBUG [BIOMOD2_Process] Generating background points...\n")
  set.seed(sp_aphia %||% 123); n_bg <- cfg_obj$background_points_n %||% 10000
  bg_pts_raw <- terra::spatSample(env_stack_curr, size=n_bg, method="random", na.rm=TRUE, xy=TRUE, warn=FALSE)
  if (nrow(bg_pts_raw) == 0) { cat(log_prefix, "ERROR [BIOMOD2_Process] Failed BG sampling.\n"); return(NULL)}
  bg_coords_current <- as.data.frame(bg_pts_raw[, c("x", "y")]); colnames(bg_coords_current) <- c("longitude", "latitude")
  cat(log_prefix, "INFO [BIOMOD2_Process] Background points:", nrow(bg_coords_current), "\n")
  
  # 8.3 Format Data
  myBiomodDataCurrent <- format_data_for_biomod2(sp_name_biomod, pres_coords_current, bg_coords_current, env_stack_curr, NULL)
  if (is.null(myBiomodDataCurrent)) { cat(log_prefix, "ERROR [BIOMOD2_Process] Data formatting failed.\n"); return(NULL)}
  
  # 8.4 Placeholder for blockCV (currently uses BIOMOD2 internal random CV)
  cat(log_prefix, "INFO [BIOMOD2_Process] Using BIOMOD2 internal random CV (blockCV placeholder).\n")
  myBiomodCVTableCurrent <- NULL # create_biomod2_block_cv_table(myBiomodDataCurrent, env_stack_curr, cfg_obj, NULL)
  
  # 8.5 Run MAXNET with 'bigboss' options
  cfg_obj$group_name_for_biomod_output <- grp_name_val # Ensure helper can make paths
  myBiomodModelOutCurrent <- run_biomod2_models_basic(
    myBiomodDataCurrent, myBiomodCVTableCurrent, 'MAXNET',
    gsub("\\.", "_", sp_name_biomod), # Use file-safe name for dir
    pred_type_suffix_base_val, # e.g., _pca (not the combined _pca_biomod2)
    grp_name_val, cfg_obj, NULL
  )
  if (is.null(myBiomodModelOutCurrent)) { cat(log_prefix, "ERROR [BIOMOD2_Process] Modeling failed.\n"); return(NULL)}
  
  # 8.6 Project Current
  algo_to_project_current <- 'MAXNET'
  cat(log_prefix, "INFO [BIOMOD2_Process] Projecting", algo_to_project_current, "to current...\n")
  projection_current <- project_biomod2_models_current(
    myBiomodModelOutCurrent, env_stack_curr,
    gsub("\\.", "_", sp_name_biomod), # Pass file-safe name for logging
    pred_type_suffix_base_val,
    algo_to_project_current, cfg_obj, NULL
  )
  if (is.null(projection_current)) { cat(log_prefix, "ERROR [BIOMOD2_Process] Projection failed.\n"); return(NULL)}
  
  # 8.7 Save Projection
  # Suffix for this basic biomod2 MAXNET run
  biomod2_model_run_suffix <- paste0(pred_type_suffix_base_val, "_biomod2_MAXNET_default")
  save_success <- save_biomod2_projection(
    projection_current, gsub("\\.", "_", sp_name_biomod), current_scenario_name,
    biomod2_model_run_suffix, cfg_obj, NULL, NULL
  )
  if(save_success) cat(log_prefix, "INFO [BIOMOD2_Process] Current projection saved.\n")
  else cat(log_prefix, "ERROR [BIOMOD2_Process] Failed to save current projection.\n")
  
  cat(log_prefix, "INFO [BIOMOD2_Process] --- Basic BIOMOD2 processing finished ---\n")
  return(list(status = "success_biomod2_basic", species = sp_name_actual, biomodModel=myBiomodModelOutCurrent))
} # End process_species_biomod2_basic

# --- 9. Run for the Selected Species ---
# For tweaking script, ensure parallel is off
config$use_parallel <- FALSE 
config$num_cores <- 1
if(exists("future")) future::plan(future::sequential) # Ensure sequential for this script

cat(paste(Sys.time(), "[INFO] --- Calling BIOMOD2 processing for species:", species_name_actual, "---\n"))
biomod2_run_result <- process_species_biomod2_basic(
  sp_name_actual = species_name_actual,
  sp_name_biomod = species_name_for_biomod_resp, # The dot-separated one for BIOMOD_FormatingData
  sp_aphia = species_aphia_id,
  occ_dir_sp = occurrence_dir,
  env_stack_curr = env_predictor_stack_current,
  pred_type_suffix_base_val = predictor_type_suffix_base, # e.g., "_pca"
  grp_name_val = group_name,
  cfg_obj = config
)

# --- 10. Process Results (Very Simplified for this Tweaking Script) ---
cat(paste(Sys.time(), "[INFO]--- BIOMOD2 Tweaking Run Result ---\n"))
if (!is.null(biomod2_run_result)) {
  cat(paste(Sys.time(), "[STATUS:", biomod2_run_result$status, "]", biomod2_run_result$species, "-", biomod2_run_result$message %||% "Completed.", "\n"))
  if (!is.null(biomod2_run_result$biomodModel)) {
    cat("  BIOMOD_Modeling output object generated.\n")
    # Example: Get evaluations if you want to see them
    # evals <- biomod2::get_evaluations(biomod2_run_result$biomodModel)
    # print(evals)
  }
} else {
  cat(paste(Sys.time(), "[ERROR] BIOMOD2 processing returned NULL overall.\n"))
}

gc(full=TRUE)
cat(paste(Sys.time(), "[INFO]--- Script biomod2_tweaking_base_run.R finished. ---\n"))
#-------------------------------------------------------------------------------