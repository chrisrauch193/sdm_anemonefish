# scripts/05_env_variable_selection.R

# This script selects the environmental variables to use for the SDM,
# and performs collinearity analysis for anemones and anemonefish, manually broken out.

library(raster)
library(usdm)
library(terra)
library(dplyr)
library(ggplot2)
library(tools)
library(ggcorrplot)
library(corrplot)
library(car)

# Load the extract_env.R functions
source("helpers/extract_env.R")

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------

# Define global variables (adjust as needed)
env_folder      <- "data/env/current"
save_location   <- "data/log"
vif_threshold   <- 5
coral_shapefile <- "data/shapefiles/WCMC008_CoralReef2018_Py_v4_1.shp"
occurrence_crs  <- "EPSG:4326"
env_crs         <- "EPSG:4326"

# Create log directory if it doesn't exist
if (!dir.exists(save_location)) {
  dir.create(save_location, recursive = TRUE)
}

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------

plot_vif_results <- function(vif_result, model, save_path = NULL) {
  df <- data.frame(Variable = names(vif_result), VIF = vif_result)
  
  # Define the color palette
  num_vars <- nrow(df)
  colors <- colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars) # Dark blue to dark red
  
  # Ensure VIF values are non-negative
  df$VIF <- pmax(0, df$VIF)
  
  p <- ggplot(df, aes(x = Variable, y = VIF)) +
    geom_bar(stat = "identity", aes(fill = Variable), show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    labs(title = "VIF Analysis Results",
         x     = "Environment Variables",
         y     = "Vif Values") + # Corrected axis labels to match image
    scale_y_continuous(limits = c(0, ceiling(max(df$VIF))),
                       breaks = seq(0, ceiling(max(df$VIF)), by = 1)) + # Match y-axis
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 8, height = 6)
    cat("VIF plot saved to", save_path, "\n")
  }
  
  return(p)
}

plot_correlation_results <- function(env_extract, save_path = NULL) {
  
  # Pearson correlation analysis
  env.cor <- round(cor(env_extract, method = "pearson"), 3)
  env.p   <- round(cor_pmat(env_extract, method = "pearson"), 3) # Get p-values
  #Setting this parameter
  par(mfrow = c(1, 1), mar = c(2, 2, 1, 1))
  # GGcorplot implementation
  cor.plot <- corrplot(corr = env.cor, type = "upper", tl.pos = "tp",
                       tl.col = "black", tl.cex = 1.2, cl.cex = 1, p.mat = env.p, insig =
                         "label_sig", sig.level = c(.01, .05), pch.cex = 1.5, pch.col = "black",
                       order = "original")
  cor.plot <- corrplot(
    corr = env.cor, type = "lower", add = TRUE, method = "number",
    tl.pos = "n", tl.col = "black",
    col = "black", tl.cex = 1.3, diag = FALSE, cl.pos = "n", pch.col =
      "black",
    number.cex = 1, number.font = 1, order = "original")
  
  # GGcorplot implementation
  if (!is.null(save_path)) {
    png(filename = save_path, width = 8, height = 6, units = "in", res = 300)
    corrplot(corr = env.cor, type = "upper", tl.pos = "tp",
             tl.col = "black", tl.cex = 1.2, cl.cex = 1, p.mat = env.p, insig =
               "label_sig", sig.level = c(.01, .05), pch.cex = 1.5, pch.col = "black",
             order = "original")
    corrplot(
      corr = env.cor, type = "lower", add = TRUE, method = "number",
      tl.pos = "n", tl.col = "black",
      col = "black", tl.cex = 1.3, diag = FALSE, cl.pos = "n", pch.col =
        "black", # cl.cex图例字体大小
      number.cex = 1, number.font = 1, order = "original"
    )
    dev.off()
    cat("Correlation plot saved to", save_path, "\n")
  }
  
  return(recordPlot()) #Returns the plot object
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

# 2b. Load occurrence data and extract environmental values
env_values_df_anemone <- load_occurrence_and_extract_env(occurrence_folder_anemone, env_rast, occurrence_crs)

# 3b. Collinearity Analysis (VIF)
env_values_df_anemone <- env_values_df_anemone[, !apply(env_values_df_anemone, MARGIN = 2, FUN = anyNA), drop = FALSE]

#Visualize scatter matrix
#Adding this new step to do scatterplot matrix now
scatterplotMatrix(env_values_df_anemone, main = "anemone correlation")

# Define a numerical code with matrix
present_anemone <- as.matrix(seq(from = 1.0, to = 1.0, length.out = nrow(env_values_df_anemone)))
env_extract1_anemone <- cbind(present_anemone, env_values_df_anemone)

# fit linear model for VIF calculation using car::vif
result_anemone <- lm(present_anemone ~ ., data = env_extract1_anemone)

# VIF calculation with car
vif_result_anemone <- car::vif(result_anemone)

# Visualize Pearson corerlation
save_path_cor <- file.path(save_location, paste0(output_prefix_anemone, "_pearson_cor.png"))
plot_correlation_results(env_values_df_anemone, save_path_cor)

# Plot VIF with new function
save_path_vif <- file.path(save_location, paste0(output_prefix_anemone, "_vif_analysis.png"))
plot_vif_results(vif_result_anemone, model = result_anemone, save_path_vif)

#Code to step the variables - this is similar to VIF step but it's for lm not VIF directly
myStep_anemone <- step(result_anemone, direction = "both")
summary(myStep_anemone)
selected_vars_anemone <- names(car::vif(myStep_anemone))

# 4. Save Results
variables_file <- file.path(save_location, paste0(output_prefix_anemone, "_selected_variables.txt"))
writeLines(selected_vars_anemone, variables_file)
cat("\nAnemone selected variables saved to", variables_file, "\n")

cat("Code complete, now running next script")

#-------------------------------------------------------------------------------
# Anemonefish Analysis
#-------------------------------------------------------------------------------

cat("--- Anemonefish Analysis ---\n")

# Define species-specific settings
occurrence_folder_anemonefish <- "data/occurrence/anemonefish"
output_prefix_anemonefish   <- "anemonefish"

# 2a. Load occurrence data and extract environmental values
env_values_df_anemonefish <- load_occurrence_and_extract_env(occurrence_folder_anemonefish, env_rast, occurrence_crs)

# 3a. Collinearity Analysis (VIF)
env_values_df_anemonefish <- env_values_df_anemonefish[, !apply(env_values_df_anemonefish, MARGIN = 2, FUN = anyNA), drop = FALSE]

#Visualize scatter matrix
#Adding this new step to do scatterplot matrix now
scatterplotMatrix(env_values_df_anemonefish, main = "anemonefish correlation")

# Define a numerical code with matrix
present_anemonefish <- as.matrix(seq(from = 1.0, to = 1.0, length.out = nrow(env_values_df_anemonefish)))
env_extract1_anemonefish <- cbind(present_anemonefish, env_values_df_anemonefish)

# fit linear model for VIF calculation using car::vif
result_anemonefish <- lm(present_anemonefish ~ ., data = env_extract1_anemonefish)

# VIF calculation with car
vif_result_anemonefish <- car::vif(result_anemonefish)

# Visualize Pearson corerlation
save_path_cor <- file.path(save_location, paste0(output_prefix_anemonefish, "_pearson_cor.png"))
plot_correlation_results(env_values_df_anemonefish, save_path_cor)

# Plot VIF with new function
save_path_vif <- file.path(save_location, paste0(output_prefix_anemonefish, "_vif_analysis.png"))
plot_vif_results(vif_result_anemonefish, model = result_anemonefish, save_path_vif)

#Code to step the variables - this is similar to VIF step but it's for lm not VIF directly
myStep_anemonefish <- step(result_anemonefish, direction = "both")
summary(myStep_anemonefish)
selected_vars_anemonefish <- names(car::vif(myStep_anemonefish))

# 4. Save Results
variables_file <- file.path(save_location, paste0(output_prefix_anemonefish, "_selected_variables.txt"))
writeLines(selected_vars_anemonefish, variables_file)
cat("\nAnemonefish selected variables saved to", variables_file, "\n")
