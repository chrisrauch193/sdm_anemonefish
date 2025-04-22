# scripts/final/basic_run.R

source("~/a0236995/sdm_anemonefish/scripts/final/cleanup.R", echo=TRUE)

cat("--- Running basic_run.R ---\n")

# --- 1. Setup: Load Config FIRST ---
# ... (keep as is) ...
if (file.exists("../config.R")) { source("../config.R") } else if (file.exists("scripts/config.R")) { source("scripts/config.R") } else { stop("FATAL: Configuration file 'scripts/config.R' not found.") }
if (!exists("config") || !is.list(config)) { stop("FATAL: 'config' list object not found or invalid after sourcing config.R") }


# --- 2. Load Required Packages ---
# ... (keep as is) ...
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(terra, sf, dplyr, readr, tools, stringr, log4r, future, furrr, progressr, ENMeval, predicts)


# --- 3. Source Helper Functions ---
# ... (keep as is) ...
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers_enmeval.R"))


# --- 4. Setup Logging ---
# ... (keep as is) ...
logger <- setup_logger(log_file = config$log_file_path, log_level = config$log_level, append = config$log_append, log_to_console = config$log_to_console, console_level = config$log_console_level)
log4r::info(logger, "--- Starting basic_run.R ---")


# --- 5. Define Group Specifics & Predictor Type ---
# ... (keep as is) ...
group_name <- "anemone"
config$group_name <- group_name
species_list_file <- config$anemone_species_list_file
occurrence_dir <- config$anemone_occurrence_dir
use_pca <- config$use_pca_predictors
predictor_type_suffix <- ifelse(use_pca, "_enmeval_pca", "_enmeval_vif")
log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors:", ifelse(use_pca, "PCA Components", "VIF-Selected Variables"), "---"))
log4r::info(logger, paste("--- ENMeval Algorithm(s):", paste(config$enmeval_algorithms, collapse=", "), "---"))
log4r::info(logger, paste("--- ENMeval Partition Method:", config$enmeval_partitions, "---"))
log4r::info(logger, paste("--- ENMeval Model Selection Metric:", config$enmeval_selection_metric, "---"))


# --- 6. Load Predictor Information ---
# ... (Keep this section as is) ...
predictor_paths_or_list <- NULL
if(use_pca) {
  pca_paths_rds <- config$pca_raster_paths_rds_path
  if (is.null(pca_paths_rds) || !file.exists(pca_paths_rds)) { log4r::fatal(logger, "PCA paths RDS missing/invalid."); stop("PCA paths RDS invalid.") }
  predictor_paths_or_list <- tryCatch(readRDS(pca_paths_rds), error = function(e) { log4r::fatal(logger, paste("Failed load PCA paths RDS:", e$message)); NULL })
  if (!is.list(predictor_paths_or_list) || length(predictor_paths_or_list) == 0) { log4r::fatal(logger, "PCA paths list empty/invalid."); stop("PCA paths list invalid.") }
  log4r::info(logger, paste("Loaded PCA raster paths for scenarios:", paste(names(predictor_paths_or_list), collapse=", ")))
} else {
  predictor_paths_or_list <- config$final_vars_vif_anemone # Use VIF list defined for anemone
  if(is.null(predictor_paths_or_list) || !is.character(predictor_paths_or_list) || length(predictor_paths_or_list) < 2) { log4r::fatal(logger, "Final VIF var list `final_vars_vif_anemone` invalid."); stop("VIF vars missing.") }
  log4r::info(logger, paste("Using VIF-selected core variables:", paste(predictor_paths_or_list, collapse=", ")))
}

# --- 7. Create Intermediate Output Dirs ---
# ... (Keep this section as is) ...
base_intermediate_model_path <- config$models_dir_intermediate
base_intermediate_results_path <- config$results_dir_intermediate
base_species_log_path <- config$species_log_dir
intermediate_models_dir <- file.path(base_intermediate_model_path, paste0(group_name, predictor_type_suffix))
intermediate_results_dir <- file.path(base_intermediate_results_path, paste0(group_name, predictor_type_suffix))
tryCatch({ dir.create(intermediate_models_dir, recursive = TRUE, showWarnings = FALSE); dir.create(intermediate_results_dir, recursive = TRUE, showWarnings = FALSE); dir.create(base_species_log_path, recursive = TRUE, showWarnings = FALSE); log4r::debug(logger, paste("Intermediate dirs created/checked:", group_name, predictor_type_suffix))
}, error = function(e) { log4r::fatal(logger, paste("Failed create intermediate dirs:", e$message)); stop("Directory creation failed.") })

# --- 8. Load Species List ---
# ... (Keep this section as is) ...
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE) }, error = function(e) { log4r::fatal(logger, paste("Failed load species list:", e$message)); stop("Species list failed.") })
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))



species_row <- species_df[1, ]
tuning_scenario <- "current"

species_row
config
predictor_paths_or_list
group_name
predictor_type_suffix
use_pca
occurrence_dir
tuning_scenario





species_name <- species_row$scientificName
species_name_sanitized <- gsub(" ", "_", species_name)
species_aphia_id <- species_row$AphiaID

species_log_file <- file.path(config$species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), paste0("[", species_name, "]"), paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
slog("INFO", paste0("--- Starting ENMeval processing (", predictor_type_suffix, ") ---"))

# --- Define File Paths ---
# ... (Keep as is) ...
intermediate_models_dir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix))
optimal_model_file <- file.path(intermediate_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
enmeval_rds_file <- file.path(intermediate_results_dir, paste0("enmevaluation_object_", species_name_sanitized, predictor_type_suffix, ".rds"))


# --- Load Predictors for Tuning Scenario ---
# ... (Keep as is) ...
slog("DEBUG", "Loading tuning predictors for scenario:", tuning_scenario)
tuning_predictor_stack <- NULL
if(use_pca){
  tuning_predictor_path <- predictor_paths_or_list[[tuning_scenario]]
  if (is.null(tuning_predictor_path) || !file.exists(tuning_predictor_path)) { msg <- paste0("Skipping: PCA stack path for tuning scenario '", tuning_scenario, "' not found."); slog("ERROR", msg); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
  tuning_predictor_stack <- tryCatch(terra::rast(tuning_predictor_path), error = function(e) {slog("ERROR", "Failed to load tuning PCA stack:", e$message); NULL})
} else {
  current_vif_vars <- generate_scenario_variable_list(predictor_paths_or_list, tuning_scenario, config)
  tuning_predictor_stack <- load_selected_env_data(tuning_scenario, current_vif_vars, config)
}
if(is.null(tuning_predictor_stack)) { msg <- paste0("Skipping: Failed load predictor stack for tuning scenario '", tuning_scenario, "'."); slog("ERROR", msg); return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg)) }
if(terra::crs(tuning_predictor_stack) == "") { msg <- paste0("Skipping: Tuning predictor stack has no CRS assigned."); slog("ERROR", msg); return(list(status = "error_tuning_predictors_crs", species = species_name, occurrence_count = NA, message = msg)) }
slog("DEBUG", "Tuning predictor stack loaded. Names:", paste(names(tuning_predictor_stack), collapse=", "))


# --- Load, Clean, and Thin Occurrences ---
# ... (Keep as is) ...
config_for_occ_load <- config; config_for_occ_load$predictor_stack_for_thinning <- tuning_predictor_stack
slog("DEBUG", "Loading/cleaning/thinning occurrences.")
occ_data_result <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger=NULL, species_log_file=species_log_file)
if (is.null(occ_data_result) || is.null(occ_data_result$coords) || occ_data_result$count < config$min_occurrences_sdm) { msg <- paste0("Skipping: Insufficient valid/thinned occurrences (found ", occ_data_result$count %||% 0, "). Required: ", config$min_occurrences_sdm); slog("WARN", msg); return(list(status = "skipped_occurrences", species = species_name, occurrence_count = occ_data_result$count %||% 0, message = msg)) }
occs_coords_df <- as.data.frame(occ_data_result$coords); colnames(occs_coords_df) <- c("longitude", "latitude")
occurrence_count_after_thinning <- occ_data_result$count
slog("INFO", "Occurrence count after clean/thin:", occurrence_count_after_thinning)


# --- Generate Background Points ---
# ... (Keep as is) ...
slog("DEBUG", "Generating background points.")
background_points_df <- generate_sdm_background(tuning_predictor_stack, config$background_points_n, config, logger=NULL, species_log_file=species_log_file, seed = species_aphia_id)
if (is.null(background_points_df)) { msg <- paste0("Skipping: Failed background point generation."); slog("ERROR", msg); return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg)) }
colnames(background_points_df) <- c("longitude", "latitude")
slog("DEBUG", "Background points generated:", nrow(background_points_df))




# enmeval_args <- list(
#   occs = occs_coords_df,
#   envs = tuning_predictor_stack,
#   bg = background_points_df,
#   tune.args = config$enmeval_tuning_settings,
#   algorithm = config$enmeval_algorithms, # Vector of algorithms
#   partitions = config$enmeval_partitions, # The method name
#   other.settings = list(pred.type = config$enmeval_pred_type, abs.auc.diff = FALSE, validation.bg = "partition"),
#   #partition.settings = list(), # Initialize empty, specific args added below
#   clamp = config$enmeval_clamp,
#   parallel = config$use_parallel && config$num_cores > 1,
#   numCores = if(config$use_parallel && config$num_cores > 1) config$num_cores else 1,
#   quiet = TRUE
# )




# prep data for custom block sampling
pts <- rbind(occs_coords_df, background_points_df)
pts$occ <- c(rep(1, nrow(occs_coords_df)), rep(0, nrow(background_points_df)))

# investigate spatial autocorrelation in the landscape to choose a suitable size for spatial blocks
# folds.cor <- cv_spatial_autocor(r = tuning_predictor_stack, x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326), column = 'occ', num_sample = 10000, plot = T, progress = T)

# generate folds
# scv <- cv_spatial(x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326), column = 'occ', r = tuning_predictor_stack, k = 5, hexagon = T, flat_top = F, size = folds.cor$range,
#                   selection = 'random', iteration = 50, progress = T, report = T, plot = T, raster_colors = terrain.colors(10, rev = T))
scv <- cv_spatial(x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326), column = 'occ', r = tuning_predictor_stack, k = 5, hexagon = T, flat_top = F, size = 20000,
                  selection = 'random', iteration = 50, progress = T, report = T, plot = T, raster_colors = terrain.colors(10, rev = T))

# plot folds
cv_plot(cv = scv, x = st_as_sf(pts, coords = c('longitude', 'latitude'), crs = 4326))

# separate occs and bg folds
fold_ids <- data.frame(fold_ids = scv$folds_ids)
folds <-  cbind(pts, fold_ids)
head(folds)

occs.folds <- folds %>% dplyr::filter(occ == 1) %>% dplyr::select('fold_ids')
bg.folds <- folds %>% dplyr::filter(occ == 0) %>% dplyr::select('fold_ids')

# export folds
saveRDS(as.vector(occs.folds)$fold_ids, 'outputs/folds/occs_fold.rds')
saveRDS(as.vector(bg.folds)$fold_ids, 'outputs/folds/bg_fold.rds')



enmeval_results <- ENMevaluate(taxon.name = 'Radianthus magnifica',
                               occs = occs_coords_df,
                               envs = tuning_predictor_stack,
                               bg = background_points_df,
                               tune.args = list(rm = seq(0.5, 5, by = 0.5),
                                                fc = c('L', 'LQ', 'H', 'LQH', 'LQHP', 'LQHPT')),
                               partitions = 'user',
                               user.grp = list(occs.grp = as.vector(occs.folds)$fold_ids,
                                               bg.grp = as.vector(bg.folds)$fold_ids),
                               algorithm = 'maxnet',
                               doClamp = T,
                               parallel = config$use_parallel && config$num_cores > 1,
                               numCores = if(config$use_parallel && config$num_cores > 1) config$num_cores else 1,
                               updateProgress = T)


i.mods <- enmeval_results

i.tune.res <- eval.results(i.mods)
i.find.opt <- i.tune.res %>%
  dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
  dplyr::filter(cbi.val.avg == max(cbi.val.avg)) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg)) %>%
  print()







# # ENM Test 2
# # # get checkerboard
# # cb2 <- get.checkerboard2(occs = occs_coords_df, bg = background_points_df, envs = tuning_predictor_stack, aggregation.factor = c(5, 5))
# # evalplot.grps(pts = occs_coords_df, pts.grp = cb2$occs.grp, envs = tuning_predictor_stack)
# # evalplot.grps(pts = background_points_df, pts.grp = cb2$bg.grp, envs = tuning_predictor_stack)
# 
# 
# # --- Define the partition methods to loop through ---
# partition_methods <- c("randomkfold", "jackknife", "block", "checkerboard")
# 
# # --- Initialize a list to store results for each partition method ---
# all_partition_results <- list()
# 
# # --- Start the loop ---
# for (current_partition in partition_methods) {
#   
#   print(paste("--- Running ENMevaluate with partition method:", current_partition, "---"))
#   
#   # Run ENMevaluate for the current partition method
#   # Use tryCatch to handle potential errors for specific partitions gracefully
#   i.mods <- tryCatch({
#     ENMevaluate(
#       # taxon.name = taxon.name, # Optional if you don't need it later
#       occs = occs_coords_df,
#       envs = tuning_predictor_stack,
#       bg = background_points_df,
#       algorithm = 'maxnet',
#       partitions = current_partition, # Use the loop variable here
#       tune.args = list(fc = c("L", "LQ", "LQH", "H"), rm = 1:3), # Corrected: Removed duplicate LQ
#       parallel = config$use_parallel && config$num_cores > 1,
#       numCores = if (config$use_parallel && config$num_cores > 1) config$num_cores else 1,
#       quiet = TRUE # Suppress ENMeval messages inside the loop if desired
#     )
#   }, error = function(e) {
#     print(paste("ERROR running ENMevaluate for partition", current_partition, ":", e$message))
#     return(NULL) # Return NULL if there was an error
#   })
#   
#   # Proceed only if ENMevaluate ran successfully
#   if (!is.null(i.mods)) {
#     # Extract results table
#     i.tune.res <- eval.results(i.mods)
#     
#     # Find the 'optimal' model based on your criteria
#     # Note: This sequential filtering might result in zero rows if no single model
#     # is the absolute best across all three metrics simultaneously.
#     # Consider alternative selection strategies if needed.
#     i.find.opt <- i.tune.res %>%
#       dplyr::filter(or.10p.avg == min(or.10p.avg, na.rm = TRUE)) %>%
#       dplyr::filter(cbi.val.avg == max(cbi.val.avg, na.rm = TRUE)) %>%
#       dplyr::filter(auc.val.avg == max(auc.val.avg, na.rm = TRUE))
#     
#     print(paste("Optimal tuning result(s) for", current_partition, ":"))
#     print(i.find.opt)
#     
#     # Store the full ENMevaluate object and the optimal results in the list
#     all_partition_results[[current_partition]] <- list(
#       enmeval_object = i.mods,
#       results_table = i.tune.res,
#       optimal_tune = i.find.opt
#     )
#     
#   } else {
#     # Store NULL if ENMevaluate failed for this partition
#     all_partition_results[[current_partition]] <- NULL
#     print(paste("Skipping result processing for partition", current_partition, "due to error."))
#   }
#   
#   # Clean up intermediate objects (optional, helps manage memory in long loops)
#   rm(i.mods, i.tune.res, i.find.opt)
#   gc() # Garbage collection
#   
#   print(paste("--- Finished processing partition:", current_partition, "---"))
#   cat("\n") # Add a newline for better readability
#   
# } # --- End of loop ---
# 
# # --- Accessing results ---
# # You can now access the results for each partition method like this:
# # print(all_partition_results$randomkfold$optimal_tune)
# # print(all_partition_results$block$results_table)
# # plot(all_partition_results$checkerboard$enmeval_object@predictions[[1]]) # Example plot
# 
# print("--- Loop finished. All specified partitions processed. ---")





# 
# ### SDMTune!!!
# # get checkerboard
# ckb2 <- get.checkerboard(occs = occs_coords_df, bg = background_points_df, envs = tuning_predictor_stack, aggregation.factor = c(10,10))
# evalplot.grps(pts = occs_coords_df, pts.grp = ckb2$occs.grp, envs = tuning_predictor_stack)
# evalplot.grps(pts = background_points_df, pts.grp = ckb2$bg.grp, envs = tuning_predictor_stack)
# 
# 
# #####  part 3 ::: run models at 20km spatial resolution ----------
# # format data for SDMtune
# glob.swd <- prepareSWD(species = 'Radianthus magnifica', env = tuning_predictor_stack, p = occs_coords_df, a = background_points_df, verbose = T)
# print(glob.swd)
# 
# # train a base model
# glob.mod <- SDMtune::train(method = 'Maxnet', data = glob.swd, folds = ckb2, iter = 5000, progress = T)
# 
# # tune model parameters // save model = F to prevent crashing
# glob.tune <- gridSearch(model = glob.mod,
#                         hypers = list(reg = seq(1,5, by = 0.5),
#                                       fc = c('l', 'lq', 'h', 'lqh', 'lqhp', 'lqhpt')),
#                         metric = 'tss',
#                         interactive = F,
#                         progress = T,
#                         save_models = F)
# 
# print(glob.tune@results)
# 
# print(glob.tune)
# print(glob.tune@models)
# 
# # get the model with optimal parameters == lqhp 1.0 == default model
# glob.opt <- glob.tune@results %>% dplyr::filter(test_TSS == max(test_TSS))
# print(glob.opt)
# 
# 
