# scripts/99_pca_plots.R
# ------------------------------------------------------------------------------
# PCA VISUALIZATIONS (Scree Plot, Biplot, Loadings)
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, tidyr, tools, ggrepel, scales, factoextra)

# --- CONFIG ---
BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
OUTPUT_DIR  <- file.path(BASE_DIR, "outputs", "figures", "pca_plots")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

PCA_MODEL_PATH <- file.path(DATA_DIR, "env", "pca_model.rds")

if(!file.exists(PCA_MODEL_PATH)) stop("PCA Model not found. Run 0c_env_var_selection.R first.")

cat("Loading PCA Model...\n")
pca_model <- readRDS(PCA_MODEL_PATH)

# --- 1. SCREE PLOT (Variance Explained) ---
cat("Generating Scree Plot...\n")

# Use factoextra for a cleaner, publication-ready scree plot
p_scree <- fviz_eig(pca_model, addlabels = TRUE, ylim = c(0, 50)) +
  labs(title = "Variance Explained by PC Axes", 
       subtitle = "Justification for selecting top 5 PCs (>95% Variance)",
       y = "Percentage of Variance Explained", x = "Principal Component") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank())

ggsave(file.path(OUTPUT_DIR, "01_Scree_Plot.png"), p_scree, width = 8, height = 6, bg = "white")

# --- 2. BIPLOT (PC1 vs PC2) ---
cat("Generating Biplot (PC1 vs PC2)...\n")

# Variable Loadings (Arrows)
loadings <- as.data.frame(pca_model$rotation[, 1:2])
loadings$var <- rownames(loadings)

# Scaling factor for arrows (heuristic)
scale_factor <- 8 

p_biplot <- ggplot() +
  # Zero lines
  geom_hline(yintercept = 0, linetype="dashed", color="grey70") +
  geom_vline(xintercept = 0, linetype="dashed", color="grey70") +
  
  # Arrows
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * scale_factor, yend = PC2 * scale_factor),
               arrow = arrow(length = unit(0.2, "cm")), color = "firebrick", size = 1) +
  
  # Labels
  geom_text_repel(data = loadings,
                  aes(x = PC1 * scale_factor, y = PC2 * scale_factor, label = var),
                  color = "black", size = 4, fontface = "bold", box.padding = 0.5) +
  
  # Styling
  labs(title = "PCA Biplot (PC1 vs PC2)",
       subtitle = "Contribution of original variables to principal components",
       x = paste0("PC1 (", round(summary(pca_model)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_model)$importance[2,2]*100, 1), "%)")) +
  coord_fixed() +
  theme_minimal(base_size = 14)

ggsave(file.path(OUTPUT_DIR, "02_Biplot_PC1_PC2.png"), p_biplot, width = 10, height = 8, bg = "white")

# --- 3. LOADINGS BAR PLOTS (PC1 to PC5) ---
cat("Generating Loadings Bar Plots...\n")

# Extract Rotation Matrix
rot <- as.data.frame(pca_model$rotation[, 1:5])
rot$Variable <- rownames(rot)

# Reshape
rot_long <- pivot_longer(rot, cols = starts_with("PC"), names_to = "PC", values_to = "Loading")

# Plot
p_loadings <- ggplot(rot_long, aes(x = Variable, y = Loading, fill = Loading > 0)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = c("firebrick", "steelblue")) +
  facet_wrap(~PC, ncol = 1) +
  coord_flip() +
  labs(title = "Variable Loadings per Principal Component", x = NULL, y = "Loading Strength") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(OUTPUT_DIR, "03_Loadings_BarPlot.png"), p_loadings, width = 8, height = 12, bg = "white")

# --- 4. KAISER-GUTTMAN CHECK ---
cat("\n--- Kaiser-Guttman Criterion (Eigenvalues > 1) ---\n")
eigenvalues <- pca_model$sdev^2
pcs_kept <- sum(eigenvalues > 1)
cat("Number of PCs with Eigenvalue > 1:", pcs_kept, "\n")
print(data.frame(PC = paste0("PC", 1:length(eigenvalues)), Eigenvalue = round(eigenvalues, 3)))

cat("\nDone. Plots saved to:", OUTPUT_DIR, "\n")