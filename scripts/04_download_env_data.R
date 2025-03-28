# scripts/04_download_env_data.R
# Download all environmental data form Bio-Oracle using obissdm library
# https://github.com/iobis/mpaeu_sdm/blob/main/codes/get_env_data.R

# List datasets to download
datasets <- c(
  "thetao_baseline_2000_2019_depthsurf",
  "so_baseline_2000_2019_depthsurf",
  "PAR_mean_baseline_2000_2020_depthsurf",
  "phyc_baseline_2000_2020_depthsurf",
  "ph_baseline_2000_2018_depthsurf",
  "sws_baseline_2000_2019_depthsurf",
  "o2_baseline_2000_2018_depthsurf",
  "KDPAR_mean_baseline_2000_2020_depthsurf",
  "no3_baseline_2000_2018_depthsurf",
  "chl_baseline_2000_2018_depthsurf",
  "tas_baseline_2000_2020_depthsurf",
  "si_baseline_2000_2018_depthsurf",
  "mlotst_baseline_2000_2019_depthsurf"
)

datasets <- c(datasets,
              gsub("depthsurf", "depthmean", datasets))


# List scenarios to download
# Most optimistic and least optimistic scenarios
future_scenarios <- c("ssp245", "ssp370")


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

# For just temperature, download also the range, ltmin and ltmax
obissdm::get_env_data(datasets = "thetao_baseline_2000_2019_depthsurf",
                      future_scenarios = future_scenarios,
                      time_steps = time_steps, variables = c("range", "ltmin", "ltmax"),
                      average_time = T)

# For Chlorophyll-a we remove the depthmean and depthmax, as for the future is
# not available
to_remove <- list.files("data/env/current", full.names = T)
to_remove <- to_remove[grepl("chl", to_remove)]
to_remove <- to_remove[grepl("depthmean|depthmax", to_remove)]
fs::file_delete(to_remove)

# Rename KDPAR for kd
to_rename <- list.files("data/env", recursive = T, full.names = T)
to_rename <- to_rename[grepl("kdpar", to_rename)]
new_names <- gsub("kdpar", "kd", to_rename)
file.rename(to_rename, new_names)

# Rename terrain_ruggedness
to_rename <- list.files("data/env/terrain/", recursive = T, full.names = T)
to_rename <- to_rename[grepl("rugg", to_rename)]
to_rename <- to_rename[!grepl("aux", to_rename)]
new_names <- gsub("terrain_ruggedness_index", "rugosity", to_rename)
edit_r <- terra::rast(to_rename)
names(edit_r) <- "rugosity"
terra::writeRaster(edit_r, new_names, overwrite = T)
fs::file_delete(to_rename)