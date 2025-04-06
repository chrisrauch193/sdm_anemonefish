# scripts/04_download_env_data.R
# Download all environmental data form Bio-Oracle using obissdm library
# https://github.com/iobis/mpaeu_sdm/blob/main/codes/get_env_data.R


# List datasets to download
datasets <- c(
  "par_mean_baseline_2000_2020_depthsurf",
  "sws_baseline_2000_2019_depthsurf",
  "thetao_baseline_2000_2019_depthsurf",
  "thetao_baseline_2000_2019_depthmax",
  "so_baseline_2000_2019_depthmax",
  "no3_baseline_2000_2018_depthmax",
  "chl_baseline_2000_2018_depthmax",
  "phyc_baseline_2000_2020_depthmax",
  "o2_baseline_2000_2018_depthmax"
)

# List scenarios to download
# Most optimistic and least optimistic scenarios
future_scenarios <- c("ssp119", "ssp585")


# Define time steps
time_steps <- list(
  current = c("2000-01-01T00:00:00Z", "2010-01-01T00:00:00Z"), #2000-2010/2010-2020
  dec50 = c("2030-01-01", "2040-01-01"), #2030-2040/2040-2050
  dec100 = c("2080-01-01", "2090-01-01") #2080-2090/2090-2100
)

# Define variables to be downloaded
variables <- c("mean")

obissdm::get_env_data(datasets = datasets, future_scenarios = future_scenarios,
                      time_steps = time_steps, variables = variables,
                      terrain_vars = c(
                        "bathymetry_mean",
                        "slope",
                        "terrain_ruggedness_index"
                      ), average_time = T)

# List datasets to download
range_datasets <- c(
  "thetao_baseline_2000_2019_depthmax",
  "o2_baseline_2000_2018_depthmax",
  "no3_baseline_2000_2018_depthmax"
)

# Download also the range, ltmin and ltmax
obissdm::get_env_data(datasets = range_datasets,
                      future_scenarios = future_scenarios,
                      time_steps = time_steps, variables = c("range", "ltmin", "ltmax"),
                      average_time = T)

# Rename terrain_ruggedness
to_rename <- list.files("data/env/terrain/", recursive = T, full.names = T)
to_rename <- to_rename[grepl("rugg", to_rename)]
to_rename <- to_rename[!grepl("aux", to_rename)]
new_names <- gsub("terrain_ruggedness_index", "rugosity", to_rename)
edit_r <- terra::rast(to_rename)
names(edit_r) <- "rugosity"
terra::writeRaster(edit_r, new_names, overwrite = T)
fs::file_delete(to_rename)

# Distance to coast layer
# Load a base layer
base <- rast("data/env/current/thetao_baseline_depthsurf_mean.tif")

coast <- base
coast[] <- NA

coast_mask <- mask(coast, base, updatevalue = 1, inverse = F)
plot(coast_mask)

coast_agg <- aggregate(coast_mask, 4, na.rm = T)

coast_dist <- terra::distance(coast_agg)
plot(coast_dist)

coast_dist <- disagg(coast_dist, 4)
coast_dist <- mask(coast_dist, base)

coast_dist <- coast_dist/1000 # to km

names(coast_dist) <- "coastdist"

writeRaster(coast_dist, "data/env/terrain/distcoast.tif", overwrite = T)






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


#-----------------------------------------------------------------------
# Organize Future Scenario Files into Decade Subdirectories
#-----------------------------------------------------------------------

message("Organizing future scenario files into dec50 and dec100 subdirectories...")

# Define the base directory where future scenario folders (sspXXX) are located
base_future_dir <- "data/env/future"

# List all scenario directories (e.g., ssp119, ssp585) within the base future directory
scenario_dirs <- list.dirs(base_future_dir, full.names = TRUE, recursive = FALSE)

# Check if any scenario directories were found
if (length(scenario_dirs) == 0) {
  warning("No scenario directories found in ", base_future_dir, ". Skipping file organization.")
} else {
  # Loop through each scenario directory found
  for (scenario_dir in scenario_dirs) {
    message(paste("Processing directory:", scenario_dir))
    
    # Define the target subdirectories for dec50 and dec100
    dec50_target_dir <- file.path(scenario_dir, "dec50")
    dec100_target_dir <- file.path(scenario_dir, "dec100")
    
    # Create the target subdirectories if they don't exist
    # recursive = TRUE ensures parent dirs are created if needed (though unlikely here)
    # showWarnings = FALSE prevents warnings if the directories already exist
    dir.create(dec50_target_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(dec100_target_dir, showWarnings = FALSE, recursive = TRUE)
    
    # List all files directly within the current scenario directory
    # Using pattern = "\\.tif$" to only list .tif files, adjust if other extensions are needed
    # recursive = FALSE ensures we only get files in the main scenario dir, not already moved ones
    files_in_scenario <- list.files(scenario_dir, full.names = TRUE, recursive = FALSE, pattern = "\\.tif$")
    
    # Filter out any potential directories listed (safer practice)
    files_in_scenario <- files_in_scenario[!file.info(files_in_scenario)$isdir]
    
    # Check if there are any files to move
    if (length(files_in_scenario) == 0) {
      message("  No .tif files found directly in ", scenario_dir, " to move.")
      next # Skip to the next scenario directory
    }
    
    message(paste("  Found", length(files_in_scenario), "files to check."))
    
    # Loop through each file found in the scenario directory
    for (current_file_path in files_in_scenario) {
      file_name <- basename(current_file_path) # Extract just the filename
      
      # Check if the filename contains "_dec50_"
      if (grepl("_dec50_", file_name, fixed = TRUE)) {
        # Construct the new path inside the dec50 directory
        new_file_path <- file.path(dec50_target_dir, file_name)
        # Move the file
        message(paste("    Moving", file_name, "to", dec50_target_dir))
        file.rename(from = current_file_path, to = new_file_path)
        
        # Check if the filename contains "_dec100_"
      } else if (grepl("_dec100_", file_name, fixed = TRUE)) {
        # Construct the new path inside the dec100 directory
        new_file_path <- file.path(dec100_target_dir, file_name)
        # Move the file
        message(paste("    Moving", file_name, "to", dec100_target_dir))
        file.rename(from = current_file_path, to = new_file_path)
        
      } else {
        # Optional: Notify about files that don't match either pattern
        # message(paste("    Skipping file (no dec50/dec100 pattern):", file_name))
      }
    } # End loop through files in one scenario
  } # End loop through scenario directories
  
  message("File organization complete.")
} # End else block (scenario directories found)

#-----------------------------------------------------------------------
# End of Script
#-----------------------------------------------------------------------




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