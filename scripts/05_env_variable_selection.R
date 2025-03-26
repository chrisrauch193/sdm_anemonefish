# scripts/05_env_variable_selection.R

# This script selects the environmental variables to use for the SDM.

library(raster)

# Define the environmental data folder
env_folder <- "data/env/current"

# Load environmental data (current conditions) - Load this ONCE
# List all raster files in the env_folder
raster_files <- list.files(env_folder, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# Filter for "mean" variables
raster_files <- raster_files[grepl("_mean\\.tif$", raster_files)]  # Only keep files ending in _mean.tif

# Stack the raster files into a single RasterStack
env_stack <- raster::stack(raster_files)

# Print the names of the layers in the RasterStack
cat("Available environmental variables:\n")
print(names(env_stack))

#-------------------------------------------------------------------------------
# USER-DEFINED SELECTION OF ENVIRONMENTAL VARIABLES
#-------------------------------------------------------------------------------

# Specify the names of the environmental variables you want to use.
# Comment out the variables you don't want to use.

selected_env_vars <- c(
  "BO_chlomean_mean",
  "BO_damean_mean",
  "BO_nitrate_mean",
  "BO_phosmean_mean",
  "BO_salinity_mean",
  "BO_sfeo_mean",
  "BO_silicate_mean",
  "BO_tempmean_mean"
  # Add or remove variables as needed
)

#-------------------------------------------------------------------------------

# Subset the env_stack to only include the selected variables
env_stack <- env_stack[[selected_env_vars]]

cat("\nSelected environmental variables:\n")
print(names(env_stack))

# Get CRS (Coordinate Reference System) after subsetting the env_stack
crs <- raster::crs(env_stack)

# You can optionally save the selected env_stack to a file
# saveRDS(env_stack, file = "data/env/selected_env_stack.rds")