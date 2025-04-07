# scripts/01_install_requirements.R
#-------------------------------------------------------------------------------
# Install and load required packages using pacman
#-------------------------------------------------------------------------------
cat("--- Loading/Installing Required Packages ---\n")

if (!require("pacman")) install.packages("pacman")

# List all required packages (add any others identified during development)
pkg_list <- c(
  "here", "dplyr", "terra", "sf", "stringr", "ggplot2", "readr", "readxl",
  "tidyr", "vegan", "ggvegan", "corrplot", "car", "usdm",
  "robis", "rgbif", "worrms", "obistools", # Added worrms/obistools for script 02
  "ENMeval", "dismo", "raster", "parallel", # Added parallel
  "devtools", "tools", "ggtext",
  "sdmtune"
)

# Check for biooracler (specific install from GitHub if needed)
# May require Rtools (Windows) or equivalent build tools (Mac/Linux)
# You might need to run this manually once if devtools fails within pacman
# if (!require("biooracler")) {
#   cat("Attempting to install biooracler from GitHub...\n")
#   tryCatch({
#     devtools::install_github("bio-oracle/biooracler", build_vignettes = FALSE)
#   }, error = function(e) {
#     warning("Failed to install biooracler automatically. Please try installing manually:\n",
#             "devtools::install_github('bio-oracle/biooracler')", call. = FALSE)
#   })
# }
# Add biooracler to the list if installation is confirmed/handled
# pkg_list <- c(pkg_list, "biooracler") # Uncomment if biooracler needed and installed


# Use pacman to load packages, installing missing ones
pacman::p_load(char = pkg_list, install = TRUE, update = FALSE)

# Check if all packages were loaded successfully
loaded_pkgs <- names(sessionInfo()$otherPkgs)
missing_pkgs <- setdiff(pkg_list, loaded_pkgs)

# Attempt to load biooracler separately if it was installed manually or previously
if (!"biooracler" %in% loaded_pkgs && "biooracler" %in% installed.packages()[,"Package"]) {
  library(biooracler)
  loaded_pkgs <- c(loaded_pkgs, "biooracler")
  missing_pkgs <- setdiff(missing_pkgs, "biooracler") # Remove from missing if loaded
}


if (length(missing_pkgs) > 0) {
  warning("The following packages could not be loaded: ",
          paste(missing_pkgs, collapse = ", "),
          ". Please install them manually.", call. = FALSE)
} else {
  cat("All required packages loaded successfully.\n")
}

# Set terra options for memory management if needed
# terraOptions(memfrac = 0.7, tempdir = file.path(base_dir, "temp")) # Example: Use 70% of RAM

#-------------------------------------------------------------------------------