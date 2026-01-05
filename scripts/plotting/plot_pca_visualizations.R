# scripts/plotting/plot_pca_visualizations.R
#-------------------------------------------------------------------------------
# Generate PCA Visualization Plots
# - Scree Plot (Variance Explained)
# - Biplot (PC1 vs PC2)
# - Individual PC Loadings Bar Plots
# - Eigenvalue check (Kaiser-Guttman)
#-------------------------------------------------------------------------------
cat("--- Running Script: PCA Visualization Plotting ---\n")

# --- 1. Setup ---
# ... (your existing setup code remains the same) ...
if (!exists("config")) {
  if (file.exists("scripts/config.R")) {
    source("scripts/config.R")
  } else if (file.exists("../scripts/config.R")) { # If in scripts/plotting
    source("../scripts/config.R")
  } else if (file.exists("../../scripts/config.R")) { # If in scripts/plotting and base is one level up
    source("../../scripts/config.R")
  } else {
    stop("config.R not found. Please ensure it's in the 'scripts' directory relative to your project root or adjust path.")
  }
  if (!exists("config")) stop("Failed to load config object.")
}


pacman::p_load(ggplot2, dplyr, tidyr, tools, ggrepel, scales)

# --- 2. Define Output Directory ---
# ... (your existing code remains the same) ...
pca_plots_output_dir <- file.path(config$log_dir_base, "pca_visualizations")
if (!dir.exists(pca_plots_output_dir)) {
  dir.create(pca_plots_output_dir, recursive = TRUE)
  cat("Created PCA plot output directory:", pca_plots_output_dir, "\n")
}

# --- 3. Load PCA Model ---
# ... (your existing code remains the same) ...
pca_model_path <- config$pca_models_rds_path # From config.R
if (!file.exists(pca_model_path)) {
  stop("PCA model RDS file not found at:", pca_model_path,
       "\nPlease run '05_preprocess_env_pca_only.R' first.")
}
pca_model <- readRDS(pca_model_path)
cat("Loaded PCA model from:", pca_model_path, "\n")

if (!inherits(pca_model, "prcomp")) {
  stop("Loaded object is not a 'prcomp' PCA model.")
}

# Print standard PCA summary
cat("\n--- Standard PCA Summary (summary(pca_model)) ---\n")
print(summary(pca_model))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# --- ADD THIS SECTION: Kaiser-Guttman Criterion Check ---
cat("\n--- Kaiser-Guttman Criterion (Eigenvalues > 1) ---\n")
if (!is.null(pca_model$sdev)) {
  eigenvalues <- pca_model$sdev^2
  cat("Eigenvalues for each Principal Component:\n")
  # Print eigenvalues with their corresponding PC number
  eigen_df <- data.frame(PC = paste0("PC", 1:length(eigenvalues)), Eigenvalue = eigenvalues)
  print(eigen_df, row.names = FALSE)
  
  num_pcs_eigen_gt_1 <- sum(eigenvalues > 1)
  cat("\nNumber of Principal Components with Eigenvalue > 1 (Kaiser-Guttman Criterion):", num_pcs_eigen_gt_1, "\n")
  
  # You can also explicitly list which PCs meet the criterion
  pcs_meeting_criterion <- paste0("PC", which(eigenvalues > 1))
  cat("PCs meeting the criterion:", paste(pcs_meeting_criterion, collapse=", "), "\n")
} else {
  cat("Could not find '$sdev' component in pca_model. Unable to calculate eigenvalues.\n")
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

n_components_to_plot <- config$n_pca_components %||% 4 # Default to 4 if not in config
# --- 4. Plot 1: Variance Explained (Scree Plot) ---
cat("Generating Scree Plot...\n")
# It's often easier to work with the summary$importance as a matrix for this
pca_importance_matrix <- summary(pca_model)$importance

variance_explained_df <- data.frame(
  PC = colnames(pca_importance_matrix), # PC names from columns of the importance matrix
  Proportion_of_Variance = as.numeric(pca_importance_matrix["Proportion of Variance", ]), # Extract row as numeric vector
  Cumulative_Proportion = as.numeric(pca_importance_matrix["Cumulative Proportion", ])  # Extract row as numeric vector
)
# Ensure PC is treated as a factor to maintain order in ggplot
variance_explained_df$PC_factor <- factor(variance_explained_df$PC, levels = variance_explained_df$PC)

p_scree <- ggplot(variance_explained_df, aes(x = PC_factor)) +
  geom_col(aes(y = Proportion_of_Variance), fill = "steelblue", alpha = 0.8) +
  geom_line(aes(y = Cumulative_Proportion, group = 1), color = "red3", linewidth = 1) +
  geom_point(aes(y = Cumulative_Proportion), color = "red3", size = 2.5) +
  geom_text(aes(y = Cumulative_Proportion, label = scales::percent(Cumulative_Proportion, accuracy = 0.1)),
            vjust = -0.8, color = "red3", size = 3.5) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  annotate("text", x = Inf, y = 0.92, label = "90% Variance Threshold", hjust = 1.05, vjust = 0, color = "darkgreen", size = 4) +
  scale_y_continuous(name = "Proportion of Variance Explained",
                     labels = scales::percent_format(accuracy = 1),
                     sec.axis = sec_axis(~., name = "Cumulative Proportion", labels = scales::percent_format(accuracy = 1)),
                     expand = expansion(mult = c(0, 0.05))) + # Ensure y-axis starts at 0
  labs(title = "PCA Variance Explained by Component", x = "Principal Component") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y.left = element_text(color = "steelblue"),
    axis.text.y.left = element_text(color = "steelblue"),
    axis.title.y.right = element_text(color = "red3"),
    axis.text.y.right = element_text(color = "red3"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

scree_plot_filename <- file.path(pca_plots_output_dir, "pca_variance_explained.png")
ggsave(scree_plot_filename, p_scree, width = 10, height = 7, dpi = 300)
cat("Saved Scree Plot to:", scree_plot_filename, "\n")


# --- 5. Plot 2: Biplot (PC1 vs PC2) ---
cat("Generating Biplot for PC1 vs PC2...\n")
loadings_data <- as.data.frame(pca_model$rotation[,1:2]) # Get PC1 and PC2 loadings
loadings_data$OriginalVariable <- rownames(loadings_data)
# Use your config's display name function
loadings_data$Variable <- sapply(loadings_data$OriginalVariable, config$get_display_name, USE.NAMES = FALSE)

# Calculate scaling factor for arrows to make them visible
# This is a simple heuristic. More sophisticated scaling might be needed depending on data.
arrow_scale <- 6 * mean(sqrt(pca_model$sdev[1:2]^2)) # Scale by sqrt of eigenvalues

p_biplot <- ggplot() +
  geom_segment(data = loadings_data,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.25, "cm")), color = "red4", linewidth=0.8) +
  geom_text_repel(data = loadings_data,
                  aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = Variable),
                  size = 3.5, color = "black",
                  segment.color = "grey50", # Color of line from point to text
                  segment.alpha = 0.7,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf) + # Allow overlaps if necessary
  labs(title = "PCA Biplot (Loadings of PC1 vs PC2)",
       x = paste0("PC1 (", scales::percent(summary(pca_model)$importance[2,1], accuracy = 0.1), ")"),
       y = paste0("PC2 (", scales::percent(summary(pca_model)$importance[2,2], accuracy = 0.1), ")")) +
  coord_fixed() + # Maintain 1:1 aspect ratio for PCs
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey70", linewidth=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey70", linewidth=0.5)


biplot_filename <- file.path(pca_plots_output_dir, "pca_biplot_pc1_pc2.png")
ggsave(biplot_filename, p_biplot, width = 10, height = 8, dpi = 300)
cat("Saved Biplot to:", biplot_filename, "\n")



# --- 6. Plot 2.5: Biplot (PC3 vs PC4) ---
cat("Generating Biplot for PC3 vs PC4...\n")
# Get PC3 and PC4 loadings by selecting columns 3 and 4
loadings_data_34 <- as.data.frame(pca_model$rotation[,3:4]) 
loadings_data_34$OriginalVariable <- rownames(loadings_data_34)
# Use your config's display name function
loadings_data_34$Variable <- sapply(loadings_data_34$OriginalVariable, config$get_display_name, USE.NAMES = FALSE)

# Calculate scaling factor for arrows using the standard deviations of PC3 and PC4
arrow_scale_34 <- 6 * mean(sqrt(pca_model$sdev[3:4]^2))

p_biplot_34 <- ggplot() +
  # Use PC3 for x-axis and PC4 for y-axis
  geom_segment(data = loadings_data_34,
               aes(x = 0, y = 0, xend = PC3 * arrow_scale_34, yend = PC4 * arrow_scale_34),
               arrow = arrow(length = unit(0.25, "cm")), color = "blue4", linewidth=0.8) +
  geom_text_repel(data = loadings_data_34,
                  aes(x = PC3 * arrow_scale_34, y = PC4 * arrow_scale_34, label = Variable),
                  size = 3.5, color = "black",
                  segment.color = "grey50", 
                  segment.alpha = 0.7,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf) +
  # Update title and axis labels for PC3 and PC4
  labs(title = "PCA Biplot (Loadings of PC3 vs PC4)",
       x = paste0("PC3 (", scales::percent(summary(pca_model)$importance[2,3], accuracy = 0.1), ")"),
       y = paste0("PC4 (", scales::percent(summary(pca_model)$importance[2,4], accuracy = 0.1), ")")) +
  coord_fixed() + # Maintain 1:1 aspect ratio for PCs
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey70", linewidth=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey70", linewidth=0.5)

# Save the new plot with a new filename
biplot_filename_34 <- file.path(pca_plots_output_dir, "pca_biplot_pc3_pc4.png")
ggsave(biplot_filename_34, p_biplot_34, width = 10, height = 8, dpi = 300)
cat("Saved Biplot to:", biplot_filename_34, "\n")





# --- 6. Plot 3: Individual PC Loadings Bar Plots ---
cat("Generating individual PC loading plots...\n")
loadings_matrix_all <- pca_model$rotation
loadings_df_all <- as.data.frame(loadings_matrix_all)
loadings_df_all$OriginalVariable <- rownames(loadings_df_all)
# Use your config's display name function
loadings_df_all$Variable <- sapply(loadings_df_all$OriginalVariable, config$get_display_name, USE.NAMES = FALSE)

# Reshape to long format for ggplot
loadings_long_df_all <- tidyr::pivot_longer(loadings_df_all,
                                            cols = starts_with("PC"),
                                            names_to = "Principal_Component",
                                            values_to = "Loading")

# Filter for the first n_pca_components
loadings_to_plot_df_all <- loadings_long_df_all %>%
  filter(Principal_Component %in% paste0("PC", 1:n_components_to_plot))

# Determine overall min and max loading for consistent x-axis limits
max_abs_loading <- max(abs(loadings_to_plot_df_all$Loading), na.rm = TRUE)
x_limits <- c(-max_abs_loading * 1.1, max_abs_loading * 1.1)


for (pc_choice in paste0("PC", 1:n_components_to_plot)) {
  pc_data <- loadings_to_plot_df_all %>%
    filter(Principal_Component == pc_choice) %>%
    mutate(Loading_Direction = ifelse(Loading > 0, "Positive", "Negative"))
  
  # Order variables by absolute loading for better visual
  pc_data$Variable <- factor(pc_data$Variable, levels = pc_data$Variable[order(abs(pc_data$Loading))])
  
  p_loadings_individual <- ggplot(pc_data, aes(x = Loading, y = Variable, fill = Loading_Direction)) +
    geom_col(alpha = 0.9) +
    scale_fill_manual(values = c("Positive" = "darkgreen", "Negative" = "firebrick3"), name = "Loading Direction") +
    scale_x_continuous(limits = x_limits, breaks = scales::pretty_breaks(n=5)) +
    labs(title = paste(pc_choice, "Loadings Plot"),
         x = "Loading Value",
         y = "Original Environmental Variable") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth=0.5)
    ) +
    geom_vline(xintercept = 0, linetype="solid", color="grey50", linewidth=0.6)
  
  
  # Dynamically adjust height based on number of variables
  num_vars_current_pc <- nrow(pc_data)
  plot_height <- max(5, 2 + num_vars_current_pc * 0.25) # Base height + per variable
  
  plot_filename_individual <- file.path(pca_plots_output_dir, paste0("pca_loadings_", tolower(pc_choice), ".png"))
  ggsave(plot_filename_individual, p_loadings_individual, width = 8, height = plot_height, units = "in", dpi = 300, limitsize = FALSE)
  cat("Saved plot:", plot_filename_individual, "\n")
}

cat("--- PCA Visualization Plotting Script Finished ---\n")


# --- 5. Plot 1: Variance Explained (Scree Plot - Individual and Cumulative) ---
cat("\n--- Generating Scree Plot (Individual & Cumulative Variance) ---\n")
pca_importance_matrix <- summary(pca_model)$importance

variance_explained_df <- data.frame(
  PC_Num = 1:ncol(pca_importance_matrix), # Numeric PC number for plotting line
  PC_Factor = factor(colnames(pca_importance_matrix), levels = colnames(pca_importance_matrix)),
  Proportion_of_Variance = as.numeric(pca_importance_matrix["Proportion of Variance", ]),
  Cumulative_Proportion = as.numeric(pca_importance_matrix["Cumulative Proportion", ])
)


# Plot for "Elbow" detection using individual variances
p_scree_individual <- ggplot(variance_explained_df, aes(x = PC_Num, y = Proportion_of_Variance)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 3) +
  geom_text(aes(label = scales::percent(Proportion_of_Variance, accuracy = 0.1)),
            vjust = -0.8, color = "darkblue", size = 3.5) +
  scale_x_continuous(breaks = variance_explained_df$PC_Num, labels = variance_explained_df$PC_Factor) + # Show all PC labels
  scale_y_continuous(name = "Proportion of Variance Explained (Individual PC)",
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Scree Plot: Individual Variance per PC (for Elbow Method)",
       x = "Principal Component") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

scree_individual_plot_filename <- file.path(pca_plots_output_dir, "pca_scree_plot_individual_variance.png")
ggsave(scree_individual_plot_filename, p_scree_individual, width = 10, height = 6, dpi = 300)
cat("Saved Scree Plot (Individual Variance for Elbow) to:", scree_individual_plot_filename, "\n")
cat("  -> Visually inspect this plot for an 'elbow' to help determine number of PCs to retain.\n")









cat("\n--- Kaiser-Guttman Criterion (Eigenvalues > 1) ---\n")
if (!is.null(pca_model$sdev)) {
  eigenvalues <- pca_model$sdev^2
  cat("Eigenvalues for each Principal Component:\n")
  # Print eigenvalues with their corresponding PC number
  eigen_df <- data.frame(PC = paste0("PC", 1:length(eigenvalues)), Eigenvalue = eigenvalues)
  print(eigen_df, row.names = FALSE)
  
  num_pcs_eigen_gt_1 <- sum(eigenvalues > 1)
  cat("\nNumber of Principal Components with Eigenvalue > 1 (Kaiser-Guttman Criterion):", num_pcs_eigen_gt_1, "\n")
  
  # You can also explicitly list which PCs meet the criterion
  pcs_meeting_criterion <- paste0("PC", which(eigenvalues > 1))
  cat("PCs meeting the criterion:", paste(pcs_meeting_criterion, collapse=", "), "\n")
} else {
  cat("Could not find '$sdev' component in pca_model. Unable to calculate eigenvalues.\n")
}