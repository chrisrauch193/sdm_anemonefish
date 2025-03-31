# scripts/06_pca_analysis.R
# Performs PCA analysis on environmental data using vegan package,
# broken out by anemones and anemonefish, then combined

library(vegan)
library(ggvegan)
library(dplyr)
library(ggplot2)
library(terra) # For raster operations
library(tools)  # file_path_sans_ext
library(tidyr) #for pivot_longer


#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------

# Define global variables (adjust as needed)
env_folder      <- "data/env/current"
save_location   <- "data/log"
coral_shapefile <- "data/shapefiles/WCMC008_CoralReef2018_Py_v4_1.shp"
occurrence_crs  <- "EPSG:4326"
env_crs         <- "EPSG:4326"

# Create log directory if it doesn't exist
if (!dir.exists(save_location)) {
  dir.create(save_location, recursive = TRUE)
}

# Load the extract_env.R functions
source("helpers/extract_env.R")

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------

#' @title Perform PCA and generate plots using vegan
#' @description This function performs PCA on environmental data extracted at
#'   species occurrence points using the vegan package. It generates and saves
#'   a PCA biplot.
#' @param env_values_df Data frame containing the extracted environmental values.
#' @param output_prefix Prefix for the output file names (e.g., "anemone").
#' @param save_location Directory to save the PCA plot.
#' @return None (saves the PCA plot to a file).
perform_vegan_pca <- function(env_values_df, output_prefix, save_location) {
  
  # Check if the dataframe is empty
  if (nrow(env_values_df) == 0) {
    cat(paste0("Skipping PCA for ", output_prefix, " because no occurrence data was extracted.\n"))
    return()
  }
  
  # --- 1. Data Standardization (Decostand) ---
  env.d <- decostand(env_values_df, method = "standardize")
  
  # --- 2. Perform PCA (rda) ---
  env.pca <- rda(env.d)
  
  # --- 3. Extract Variance Explained ---
  a <- summary(env.pca)
  b <- a$cont$importance
  prop_explained <- b[2, ] # Get the proportion of variance explained for all axes
  num_components <- length(prop_explained)
  
  # --- 4. Extract Scores and Scale Loadings ---
  
  # Site Scores (Occurrence Points)
  site_scores <- scores(env.pca, display = "sites", choices = c(1, 2))
  site_df <- data.frame(site_scores)
  
  # Species Scores (Environmental Variables - Loadings)
  species_scores <- scores(env.pca, display = "species", choices = c(1, 2))
  species_df <- data.frame(species_scores)
  
  # Scaling factor for species scores (loadings)
  
  # Calculate scaling factor: ratio of site SD to species SD
  scaling_factor <- min(
    (max(site_df$PC1) - min(site_df$PC1)) / (max(species_df$PC1) - min(species_df$PC1)),
    (max(site_df$PC2) - min(site_df$PC2)) / (max(species_df$PC2) - min(species_df$PC2))
  ) * 0.9 #scale down by 10%
  # Apply scaling to species scores
  species_df$PC1 <- species_df$PC1 * scaling_factor
  species_df$PC2 <- species_df$PC2 * scaling_factor
  
  # --- 5. Create ggplot2 Plot ---
  pca_plot <- ggplot() +
    geom_point(data = site_df, aes(x = PC1, y = PC2), color = "salmon", size = 1) + # Site scores
    
    geom_segment(data = species_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.25, "cm")), color = "black") + # Species scores as arrows
    geom_text(data = species_df, aes(x = PC1, y = PC2, label = rownames(species_df)),
              hjust = 0.5, vjust = 0.5, size = 3, color = "black") + # Variable labels
    
    labs(x = paste("PCA Axis 1 (", round(prop_explained[1] * 100, 2), "%)", sep = ""),
         y = paste("PCA Axis 2 (", round(prop_explained[2] * 100, 2), "%)", sep = ""),
         title = paste0(output_prefix, " PCA Biplot")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          text = element_text(size=13))
  
  
  # --- 6. Stacked Bar Plot (Contribution per Variable per Component) ---
  
  # Extract variable loadings for the first four PCs
  var_loadings <- scores(env.pca, display = "species", choices = 1:4)
  var_loadings_df <- as.data.frame(var_loadings)
  var_loadings_df$Variable <- rownames(var_loadings_df)
  
  # Reshape the data for easier plotting
  var_loadings_long <- var_loadings_df %>%
    pivot_longer(cols = starts_with("PC"),
                 names_to = "Principal Component",
                 values_to = "Contribution")
  
  # Factorize Principal Component for consistent order
  var_loadings_long$`Principal Component` <- factor(var_loadings_long$`Principal Component`, levels = paste0("PC", 1:4))
  
  # Define color palette
  #Color scheme
  num_vars <- length(unique(var_loadings_long$Variable))
  colors <- colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars)
  
  contribution_plot <- ggplot(var_loadings_long, aes(x = `Principal Component`, y = Contribution, fill = Variable)) +
    geom_bar(stat = "identity", position = "stack") + # Stack bars
    scale_fill_manual(values = colors) +
    labs(title = paste0(output_prefix, " - Variable Contribution to PCs"),
         x = "Principal Component",
         y = "Contribution") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(),
          panel.grid.minor = element_blank()) #Remove grid lines
  
  # --- 7. Save the Plots ---
  pca_save_path <- file.path(save_location, paste0(output_prefix, "_vegan_pca.png"))
  ggsave(pca_save_path, plot = pca_plot, width = 8, height = 6)
  
  contribution_save_path <- file.path(save_location, paste0(output_prefix, "_contribution_plot.png"))
  ggsave(contribution_save_path, plot = contribution_plot, width = 10, height = 7)
  
  
  cat(paste("Vegan PCA plot saved to", pca_save_path, "\n"))
  cat(paste("Variable contribution plot saved to", contribution_save_path, "\n"))
}
#-------------------------------------------------------------------------------
# Main Analysis
#-------------------------------------------------------------------------------

# 1. Load and prepare environmental data
env_rast <- load_and_prepare_env_data(env_folder, coral_shapefile)
if (is.null(env_rast)) {
  stop("Failed to load and prepare environmental data.")
}

#-------------------------------------------------------------------------------
# Anemone Analysis
#-------------------------------------------------------------------------------

cat("--- Anemone Analysis ---\n")

# Define species-specific settings
occurrence_folder_anemone <- "data/occurrence/anemone"
output_prefix_anemone   <- "anemone"

# 2a. Load occurrence data and extract environmental values
env_values_df_anemone <- load_occurrence_and_extract_env(occurrence_folder_anemone, env_rast, occurrence_crs)

# 3a. Perform PCA and generate plots
perform_vegan_pca(env_values_df_anemone, output_prefix_anemone, save_location)

#-------------------------------------------------------------------------------
# Anemonefish Analysis
#-------------------------------------------------------------------------------

cat("--- Anemonefish Analysis ---\n")

# Define species-specific settings
occurrence_folder_anemonefish <- "data/occurrence/anemonefish"
output_prefix_anemonefish   <- "anemonefish"

# 2b. Load occurrence data and extract environmental values
env_values_df_anemonefish <- load_occurrence_and_extract_env(occurrence_folder_anemonefish, env_rast, occurrence_crs)

# 3b. Perform PCA and generate plots
perform_vegan_pca(env_values_df_anemonefish, output_prefix_anemonefish, save_location)

#-------------------------------------------------------------------------------
# Combined Analysis (Anemone + Anemonefish)
#-------------------------------------------------------------------------------

cat("--- Combined Analysis (Anemone + Anemonefish) ---\n")

# 1. Combine the dataframes
# (Assuming both dataframes have the same column structure)
if (nrow(env_values_df_anemone) > 0 && nrow(env_values_df_anemonefish) > 0) {
  env_values_combined <- rbind(env_values_df_anemone, env_values_df_anemonefish)
  output_prefix_combined <- "anemone_anemonefish_combined"
  
  # 2. Perform PCA and generate plots
  perform_vegan_pca(env_values_combined, output_prefix_combined, save_location)
} else {
  cat("Skipping combined analysis because one or both species have no occurrence data.\n")
}

cat("Code complete.\n")