# scripts/05_env_variable_selection.R

# This script selects the environmental variables to use for the SDM,
# and performs collinearity analysis using VIF.

library(raster)
library(usdm)
library(terra)
library(dplyr)
library(ggplot2)

# Define the environmental data folder
env_folder <- "data/env/current"
save_location <- "data/log"

# Create log directory if it doesn't exist
if (!dir.exists(save_location)) {
  dir.create(save_location, recursive = TRUE)
}

#-------------------------------------------------------------------------------
# Helper function for visualizing VIF results
#-------------------------------------------------------------------------------

plot_vif_results <- function(vif_result, save_path = NULL) {
  df <- data.frame(Variable = vif_result@results[, "Variables"], VIF = vif_result@results[, "VIF"])
  
  p <- ggplot(df, aes(x = Variable, y = VIF)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = vif_result@threshold, color = "red", linetype = "dashed") +
    annotate("text", x = Inf, y = vif_result@threshold, label = paste("VIF Threshold =", vif_result@threshold), vjust = -0.5, hjust = 1, color = "red") +
    coord_flip() +
    labs(title = "VIF Analysis Results", x = "Environmental Variable", y = "VIF Value") +
    theme_minimal()
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 8, height = 6)
    cat("VIF plot saved to", save_path, "\n")
  }
  
  return(p)
}

#-------------------------------------------------------------------------------
# Function to load and prepare environmental data
#-------------------------------------------------------------------------------

load_and_prepare_env_data <- function(env_folder) {
  tryCatch({
    # List all raster files in the env_folder
    raster_files <- list.files(env_folder, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
    
    # Stack the raster files into a single RasterStack
    env_stack <- raster::stack(raster_files)
    
    # Print the names of the layers in the RasterStack
    cat("Available environmental variables:\n")
    print(names(env_stack))
    
    # Convert raster stack to terra rast
    env_rast <- terra::rast(env_stack)
    return(env_rast)
    
  }, error = function(e) {
    cat("\nError loading environmental data:", e$message, "\n")
    return(NULL) # or handle the error as appropriate
  })
}

# Load the environmental data
env_rast <- load_and_prepare_env_data(env_folder)

if(is.null(env_rast)){
  stop("Failed to load environmental data.")
}


#-------------------------------------------------------------------------------
# USER-DEFINED SELECTION OF ENVIRONMENTAL VARIABLES (BROAD)
#-------------------------------------------------------------------------------

# Specify the names of the environmental variables you initially want to consider.
# Comment out the variables you don't want to consider at all.

initial_env_vars <- c(
  "chl_mean.1",
  "no3_mean.1",
  "o2_mean.1",
  "ph_mean.1",
  "phyc_mean.1",
  "so_mean.1",
  "sws_mean.1",
  "thetao_mean.1"
  # "bathymetry_mean", # Make sure bathymetry data is available
  # "slope",           # Make sure slope data is available
  # "rugosity"          # Make sure rugosity data is available
  # Add or remove variables as needed
)

# Subset the env_rast to only include the initially selected variables
env_rast <- env_rast[[initial_env_vars[initial_env_vars %in% names(env_rast)]]]

cat("\nInitially selected environmental variables:\n")
print(names(env_rast))


#-------------------------------------------------------------------------------
# COLLINEARITY ANALYSIS USING VIF
#-------------------------------------------------------------------------------

# Perform VIF analysis to identify and remove collinear variables
vif_threshold <- 5 # Adjust this threshold as needed


#convert the SpatRaster to a data frame so usdm::vif can use it
env_df <- terra::as.data.frame(env_rast, xy = FALSE, na.rm = TRUE)

# Check for missing values in the data frame
if (any(is.na(env_df))) {
  cat("\nMissing values found in environmental data. Imputing with column means.\n")
  for (i in 1:ncol(env_df)) {
    if (any(is.na(env_df[, i]))) {
      env_df[is.na(env_df[, i]), i] <- mean(env_df[, i], na.rm = TRUE)
    }
  }
}

#Check if any columns have zero variance
zero_variance_cols <- apply(env_df, 2, var) == 0
if (any(zero_variance_cols)) {
  cat("\nVariables with zero variance found:\n")
  print(names(env_df)[zero_variance_cols])
  cat("\nRemoving these variables before VIF analysis.\n")
  env_df <- env_df[, !zero_variance_cols]
  env_rast <- env_rast[[names(env_df)]]
}


# Perform VIF analysis
vif_result <- tryCatch({
  usdm::vifstep(env_df, th = vif_threshold)
}, error = function(e) {
  cat("\nError during VIF analysis:", e$message, "\n")
  return(NULL)
})


# Print VIF results
if (!is.null(vif_result)) {
  cat("\nVIF analysis results:\n")
  print(vif_result)
  
  # Generate and save the VIF plot
  vif_plot_path <- file.path(save_location, paste0("vif_plot_", format(Sys.Date(), "%Y%m%d"), ".png"))
  vif_plot <- plot_vif_results(vif_result, save_path = vif_plot_path)
  
  # Get the names of the variables to keep
  variables_to_keep <- vif_result@variables
  
  # Subset the raster stack to keep only the selected variables
  env_rast <- env_rast[[variables_to_keep]]
  
  cat("\nEnvironmental variables retained after VIF analysis:\n")
  print(names(env_rast))
} else {
  cat("\nVIF analysis failed. Keeping all initially selected variables.\n")
  variables_to_keep <- names(env_rast) # Keep all initial variables if VIF fails
}

# Save the VIF results to a file
if (!is.null(vif_result)) {
  saveRDS(vif_result, file = file.path(save_location, paste0("vif_result_", format(Sys.Date(), "%Y%m%d"), ".rds")))
  cat("\nVIF results saved to:", file.path(save_location, paste0("vif_result_", format(Sys.Date(), "%Y%m%d"), ".rds")), "\n")
}

# Save the final environmental variables (raster stack)
saveRDS(env_rast, file = "data/env/final_env_rast.rds")

cat("\nFinal environmental raster saved to: data/env/final_env_rast.rds\n")

# The env_stack is now the raster stack with collinear variables removed.