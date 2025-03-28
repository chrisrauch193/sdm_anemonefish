# scripts/06_env_variable_pca.R
# This script performs PCA on the selected environmental variables.

library(raster)
library(stats)
library(vegan) # For decostand and rda
library(terra)

# Load the selected environmental variables
if (!exists("env_rast")) {
  env_rast <- readRDS("data/env/final_env_rast.rds")
  if (!inherits(env_rast, "SpatRaster")) {
    stop("env_rast is not a SpatRaster object. Ensure it was created correctly in 05_1_env_variable_selection.R.")
  }
}

#-------------------------------------------------------------------------------
# USER-DEFINED PCA OPTION
#-------------------------------------------------------------------------------
perform_pca <- TRUE  # Set to TRUE to perform PCA, FALSE to skip

#-------------------------------------------------------------------------------

if (perform_pca) {
  # Convert raster to data frame
  env_df <- terra::as.data.frame(env_rast, xy = FALSE, na.rm = TRUE)
  
  # Handle missing values (replace with mean) - IMPORTANT for PCA
  for (i in 1:ncol(env_df)) {
    if (any(is.na(env_df[, i]))) {
      env_df[is.na(env_df[, i]), i] <- mean(env_df[, i], na.rm = TRUE)
    }
  }
  
  # Standardize the data (mean 0, sd 1) using vegan::decostand
  env_st <- vegan::decostand(env_df, method = "standardize")
  
  # Perform PCA using vegan::rda (more suitable for ecological data)
  env_pca <- vegan::rda(env_st)
  
  # Get summary of PCA results
  pca_summary <- summary(env_pca)
  print(pca_summary)
  
  # Calculate explained variance for each PC axis
  explained_variance <- pca_summary$cont$importance[2, ] * 100
  cat("\nExplained variance per PC axis:\n")
  print(explained_variance)
  
  # Determine the number of PCs to retain (e.g., based on explained variance)
  cumulative_variance <- cumsum(explained_variance)
  n_pcs_to_retain <- which(cumulative_variance >= 95)[1]
  
  cat("\nNumber of PCs retained (explaining >= 95% variance):", n_pcs_to_retain, "\n")
  
  # Extract the scores for the retained PCs
  pca_scores <- scores(env_pca, choices = 1:n_pcs_to_retain, display = "sites")
  print(head(pca_scores))
  
  # Create a new SpatRaster from the selected PCs
  pca_rasters <- terra::rast(nrow = nrow(env_rast), ncol = ncol(env_rast),
                             xmin = xmin(env_rast), xmax = xmax(env_rast),
                             ymin = ymin(env_rast), ymax = ymax(env_rast),
                             crs = crs(env_rast), nlyr = n_pcs_to_retain)
  names(pca_rasters) <- paste0("PC", 1:n_pcs_to_retain)
  
  # Fill the SpatRaster layers with PCA scores using a loop
  for (i in 1:n_pcs_to_retain) {
    pca_values <- predict(env_pca, newdata = env_st)[, i]  # Predict PCA values
    terra::values(pca_rasters[[i]]) <- pca_values             # Assign values to raster layer
  }
  
  cat("\nPCA Raster:\n")
  print(pca_rasters)
  
  # Optionally save the PCA results and the PCA raster stack
  saveRDS(env_pca, file = "data/env/env_pca.rds")
  saveRDS(pca_rasters, file = "data/env/pca_rasters.rds")
  
  # Update env_rast to use the PCA rasters
  env_rast <- pca_rasters
  
} else {
  cat("\nPCA skipped. Using selected environmental variables.\n")
}

#The final list of covariates
cat("\nThe final list of covariates:\n")
print(names(env_rast))

# The env_rast is now either the selected environmental variables or the PCA rasters.