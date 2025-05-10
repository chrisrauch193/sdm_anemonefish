# scripts/data_retrieval/04_download_env_data.R
#-------------------------------------------------------------------------------
# Download Environmental Data using logic from the original script.
# Uses obissdm::get_env_data (ensure package is installed).
# Uses variable lists and paths defined in config.R.
# Prepares terrain layers and copies missing future vars from current.
# Assumes obissdm downloads to data/env/current, data/env/future/sspXXX, data/env/terrain
#-------------------------------------------------------------------------------
cat("--- Running Script 04: Download Environmental Data (User's Original Logic) ---\n")

# Ensure config is loaded if running standalone
if (!exists("config")) { source("scripts/config.R"); if (!exists("config")) stop("Failed load config.") }

# Ensure required packages are loaded
pacman::p_load(terra, fs)
if (!requireNamespace("obissdm", quietly = TRUE)) {
  warning("Package 'obissdm' not found. Trying install from GitHub..."); tryCatch({
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes"); remotes::install_github("iobis/obissdm"); library(obissdm)
  }, error = function(e){ stop("Failed install/load 'obissdm'. Error: ", e$message) })
} else { library(obissdm) }

# --- Parameters Directly from Config ---
datasets_mean <- config$original_datasets_mean
range_datasets <- config$original_range_datasets
future_scenarios_download <- config$original_future_scenarios # SSP names for download
time_steps <- config$original_time_steps
terrain_vars_for_download <- config$original_terrain_vars_for_download

# --- Download Logic (Matching Original Script Calls) ---
cat("\nDownloads starting... Check messages from obissdm::get_env_data for details.\n")
# 1. Download Mean Vars + Terrain
tryCatch({ obissdm::get_env_data(datasets = datasets_mean, future_scenarios = future_scenarios_download, time_steps = time_steps, variables = c("mean"), terrain_vars = terrain_vars_for_download, average_time = TRUE)
  cat("Initial download call completed.\n") }, error = function(e) { warning("Initial download Error: ", e$message) })
# 2. Download Range, Min, Max
tryCatch({ obissdm::get_env_data(datasets = range_datasets, future_scenarios = future_scenarios_download, time_steps = time_steps, variables = c("range", "ltmin", "ltmax"), terrain_vars = NULL, average_time = TRUE)
  cat("Range/min/max download call completed.\n") }, error = function(e) { warning("Range download Error: ", e$message) })
cat("\nDownloads attempted.\n")

# --- Post-processing Steps ---

# 3. Rename terrain_ruggedness
cat("\nProcessing Terrain: Renaming ruggedness to rugosity...\n")
# (Keep the existing renaming logic - it should work if files are in config$terrain_folder)
terrain_target_folder <- config$terrain_folder
ruggedness_download_name <- "terrain_ruggedness_index"
final_rugosity_name <- "rugosity"
if (dir.exists(terrain_target_folder)) {
  pattern_to_find <- paste0(ruggedness_download_name, ".*\\.tif$")
  to_rename <- list.files(terrain_target_folder, pattern = pattern_to_find, recursive = TRUE, full.names = TRUE)
  to_rename <- to_rename[!grepl("\\.aux\\.xml$", to_rename, ignore.case = TRUE)]
  if (length(to_rename) > 0) {
    renamed_count <- 0; for(f in to_rename){ new_name <- gsub(ruggedness_download_name, final_rugosity_name, f, ignore.case = TRUE); tryCatch({ edit_r <- terra::rast(f); names(edit_r) <- final_rugosity_name; terra::writeRaster(edit_r, new_name, overwrite = TRUE, gdal=c("COMPRESS=LZW")); if(file.exists(new_name) && file.size(new_name) > 0){ fs::file_delete(f); aux_file_old <- paste0(tools::file_path_sans_ext(f), ".aux.xml"); if(fs::file_exists(aux_file_old)) fs::file_delete(aux_file_old); renamed_count <- renamed_count + 1 } else { warning("Failed write ", new_name) } }, error = function(e){ warning("Error processing ", f, ": ", e$message) }) }; cat("Processed", renamed_count, "ruggedness file(s).\n")
  } else { cat("No '", ruggedness_download_name, ".tif' files found in ", terrain_target_folder, "\n", sep="") }
} else { warning("Target terrain folder not found: ", terrain_target_folder) }


# 4. Distance to coast layer
cat("\nProcessing Terrain: Calculating distance to coast...\n")
# (Keep the existing distance calculation logic - relies on files being in correct place)
distcoast_file <- file.path(config$terrain_folder, "distcoast.tif")
if (!file.exists(distcoast_file)) {
  base_layer_filename <- "thetao_baseline_2000_2019_depthsurf_mean.tif" # Check name
  base_layer_path <- file.path(config$scenario_folder_map$current, base_layer_filename)
  if (file.exists(base_layer_path)) {
    tryCatch({ base <- terra::rast(base_layer_path); cat("Using base raster:", base_layer_path, "\n"); land_mask <- terra::classify(base[[1]], cbind(NA, 1), others = NA)
    if(is.na(terra::global(land_mask, "max", na.rm=TRUE)$max) || terra::global(land_mask, "max", na.rm=TRUE)$max != 1){ warning("Cannot create land mask.")
    } else { coast_agg <- terra::aggregate(land_mask, fact = 4, fun = "max", na.rm = TRUE); coast_dist_agg <- terra::distance(coast_agg); coast_dist <- terra::disagg(coast_dist_agg, fact = 4, method = "bilinear"); coast_dist <- terra::resample(coast_dist, base, method="bilinear"); coast_dist <- terra::mask(coast_dist, base[[1]]); coast_dist_km <- coast_dist / 1000; names(coast_dist_km) <- "distcoast"; dir.create(dirname(distcoast_file), showWarnings = FALSE, recursive = TRUE); terra::writeRaster(coast_dist_km, distcoast_file, overwrite = TRUE, gdal=c("COMPRESS=LZW")); cat("Distance to coast saved:", distcoast_file, "\n") }
    }, error = function(e){ warning("Error calculating distance: ", e$message) })
  } else { warning("Base layer not found: ", base_layer_path) }
} else { cat("Distance to coast file exists. Skipping calculation.\n") }


# 5. --- REVISED v2: Copy specified current vars & RENAME for future directories ---
cat("\nCopying & renaming current-only variables (CHL, PAR) to future directories...\n")
vars_to_copy <- config$vars_to_copy_to_future # Basenames without extension
current_dir <- config$scenario_folder_map$current
time_tags_map <- list(dec50 = "_dec50_", dec100 = "_dec100_") # Tags corresponding to time_steps keys

copied_count_total <- 0
skipped_count_total <- 0

# Loop through each unique future *SSP* scenario name (e.g., "ssp119", "ssp585")
for (ssp_scenario in config$original_future_scenarios) {
  
  # Find the corresponding target directory (e.g., data/env/future/ssp119)
  # Assumes the directory name matches the ssp_scenario name
  future_ssp_dir <- file.path(config$env_data_dir, "future", ssp_scenario)
  
  if (!dir.exists(future_ssp_dir)) {
    warning("Target future directory does not exist for SSP ", ssp_scenario, ", cannot copy files: ", future_ssp_dir)
    next # Skip to next SSP scenario
  }
  cat("  Checking/copying to SSP directory:", future_ssp_dir, "\n")
  copied_count_dir <- 0
  skipped_count_dir <- 0
  
  # Loop through each variable to copy (e.g., chl_..., par_...)
  for (var_basename in vars_to_copy) {
    source_file <- file.path(current_dir, paste0(var_basename, ".tif"))
    
    if (file.exists(source_file)) {
      # Loop through the time tags (_dec50_, _dec100_)
      for (time_tag_name in names(time_tags_map)) {
        time_tag_value <- time_tags_map[[time_tag_name]]
        time_tag_clean <- gsub("^_|_$", "", time_tag_value) # e.g., "dec50"
        
        # Construct the NEW filename for the future directory
        # Replace 'baseline_YYYY_YYYY' with 'sspXXX_decYY'
        # 1. Remove baseline year part
        new_basename_step1 <- stringr::str_replace(var_basename, "_baseline_\\d{4}_\\d{4}", "")
        # 2. Insert SSP and time tag (e.g., before _depthsurf/_depthmax or before _mean)
        #    We need a robust way based on expected structure. Assume structure is VAR_DEPTH_STAT
        parts <- stringr::str_split(new_basename_step1, "_")[[1]]
        if (length(parts) >= 3) { # e.g., chl, depthmax, mean
          var_part <- parts[1]
          depth_part <- parts[grepl("depth", parts)] # Find depth part
          stat_part <- parts[length(parts)] # Assume last part is stat
          # Construct: var_ssp_depth_tag_stat
          new_basename <- paste(var_part, ssp_scenario, depth_part, time_tag_clean, stat_part, sep = "_")
        } else {
          warning("Unexpected filename structure for '", var_basename, "', using simpler renaming.")
          new_basename <- paste(new_basename_step1, ssp_scenario, time_tag_clean, sep = "_")
        }
        
        dest_file <- file.path(future_ssp_dir, paste0(new_basename, ".tif"))
        
        # Copy only if destination doesn't exist
        if (!file.exists(dest_file)) {
          tryCatch({
            fs::file_copy(path = source_file, new_path = dest_file, overwrite = TRUE)
            aux_source <- paste0(source_file, ".aux.xml")
            if(fs::file_exists(aux_source)) fs::file_copy(aux_source, paste0(dest_file, ".aux.xml"), overwrite = TRUE)
            copied_count_dir <- copied_count_dir + 1
            cat("    Copied", basename(source_file), "as", basename(dest_file), "\n")
          }, error = function(e){
            warning("Failed to copy/rename ", basename(source_file), " to ", basename(dest_file), ": ", e$message)
          })
        } else {
          skipped_count_dir <- skipped_count_dir + 1
        }
      } # End loop time_tags
    } else {
      warning("Source file to copy not found in 'current' directory: ", basename(source_file))
    }
  } # End loop vars_to_copy
  cat("    --> Copied:", copied_count_dir, "file versions. Skipped:", skipped_count_dir, "(already exist).\n")
  copied_count_total <- copied_count_total + copied_count_dir
  skipped_count_total <- skipped_count_total + skipped_count_dir
}
cat("Total file versions copied & renamed to future directories:", copied_count_total, "\n")
# 6. Verification Step (Optional)
cat("\nVerifying expected directory structure...\n")
# (Keep verification logic as is)
if(dir.exists(config$scenario_folder_map$current)){ cat("  Found:", config$scenario_folder_map$current, "\n") } else { cat("  MISSING:", config$scenario_folder_map$current, "\n") }
future_base_dir = dirname(config$scenario_folder_map$ssp119_2050); if(dir.exists(future_base_dir)){ cat("  Found:", future_base_dir, "\n"); if(length(config$original_future_scenarios) > 0){ ssp119_dir = file.path(future_base_dir, "ssp119"); if(dir.exists(ssp119_dir)) cat("  Found:", ssp119_dir, "\n") else cat("  MISSING:", ssp119_dir, "\n"); ssp585_dir = file.path(future_base_dir, "ssp585"); if(dir.exists(ssp585_dir)) cat("  Found:", ssp585_dir, "\n") else cat("  MISSING:", ssp585_dir, "\n") } } else if (length(config$original_future_scenarios) > 0) { cat("  MISSING base future directory:", future_base_dir, "\n") }
if(dir.exists(config$terrain_folder)){ cat("  Found:", config$terrain_folder, "\n") } else { cat("  MISSING:", config$terrain_folder, "\n") }

cat("\n--- Script 04 finished. ---\n")
#-------------------------------------------------------------------------------