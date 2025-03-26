# scripts/06_env_variable_pca.R
# This script performs PCA on the selected environmental variables.

library(raster)
library(stats)

# Load the selected environmental variables (if not already in memory)
# env_stack <- readRDS("data/env/selected_env_stack.rds") #Uncomment this line if you chose to save the output from the previous script
if (!exists("env_stack")) {
  stop("env_stack not found. Please run 05_1_env_variable_selection.R first.")
}

#-------------------------------------------------------------------------------
# USER-DEFINED PCA OPTION
#-------------------------------------------------------------------------------
perform_pca <- TRUE  # Set to TRUE to perform PCA, FALSE to skip

#-------------------------------------------------------------------------------

if (perform_pca) {
  # Convert raster stack to a data frame
  env_df <- raster::rasterToPoints(env_stack)
  env_df <- as.data.frame(env_df[, -(1:2)]) # Remove x and y coordinates
  
  # Handle missing values (replace with mean) - IMPORTANT for PCA
  for (i in 1:ncol(env_df)) {
    if (any(is.na(env_df[, i]))) {
      env_df[is.na(env_df[, i]), i] <- mean(env_df[, i], na.rm = TRUE)
    }
  }
  
  # Perform PCA
  pca_result <- prcomp(env_df, scale. = TRUE) # scale. = TRUE for standardized PCA
  
  # Print PCA summary
  print(summary(pca_result))
  
  # Determine the number of PCs to retain (e.g., based on explained variance)
  # For example, keep PCs that explain 95% of the variance:
  cumulative_variance <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
  n_pcs_to_retain <- which(cumulative_variance >= 0.95)[1]
  
  cat("\nNumber of PCs retained:", n_pcs_to_retain, "\n")
  
  # Create a new raster stack from the selected PCs
  pca_rasters <- raster::predict(env_stack, pca_result, na.rm = TRUE, index = 1:n_pcs_to_retain)
  
  # Rename the layers
  names(pca_rasters) <- paste0("PC", 1:n_pcs_to_retain)
  
  cat("\nPCA Raster Stack:\n")
  print(pca_rasters)
  
  # Optionally save the PCA results and the PCA raster stack
  # saveRDS(pca_result, file = "data/env/pca_result.rds")
  # saveRDS(pca_rasters, file = "data/env/pca_rasters.rds")
  
  # Update env_stack to use the PCA rasters
  env_stack <- pca_rasters
  
} else {
  cat("\nPCA skipped. Using selected environmental variables.\n")
}

#The final list of covariates
cat("\nThe final list of covariates:\n")
print(names(env_stack))

# The env_stack is now either the selected environmental variables or the PCA rasters.