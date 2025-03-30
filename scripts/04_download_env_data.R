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