# scripts/99_pca_plots.R
# ------------------------------------------------------------------------------
# FIGURE S1-S3: PCA VISUALIZATIONS (Thesis Ready)
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

# --- CALCULATE VARIANCE METRICS ---
eig_val <- get_eigenvalue(pca_model)
cum_var_4 <- eig_val[4, "cumulative.variance.percent"]

# --- 1. SCREE PLOT (Variance Explained) ---
cat("Generating Scree Plot...\n")

p_scree <- fviz_eig(pca_model, addlabels = TRUE, ylim = c(0, 50), 
                    barfill = ifelse(1:nrow(eig_val) <= 4, "#2E9FDF", "#999999"), 
                    barcolor = ifelse(1:nrow(eig_val) <= 4, "#2E9FDF", "#999999")) +
  labs(title = "Variance Explained by Principal Components", 
       subtitle = paste0("Selected: PC1–PC4 (Cumulative Variance: ", round(cum_var_4, 1), "%)"),
       y = "Percentage of Variance (%)", x = "Principal Component") +
  geom_hline(yintercept = 100/nrow(eig_val), linetype="dashed", color="red") + 
  annotate("text", x=6, y=12, label="Kaiser Criterion cutoff", color="red", size=3) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(), plot.title = element_text(face="bold"))

ggsave(file.path(OUTPUT_DIR, "01_Scree_Plot.png"), p_scree, width = 8, height = 6, bg = "white")

# --- 2. BIPLOT (PC1 vs PC2) ---
cat("Generating Biplot (PC1 vs PC2)...\n")

loadings <- as.data.frame(pca_model$rotation[, 1:2])
loadings$var <- rownames(loadings)
scale_factor <- 8 

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype="dashed", color="grey70") +
  geom_vline(xintercept = 0, linetype="dashed", color="grey70") +
  # Arrows
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * scale_factor, yend = PC2 * scale_factor),
               arrow = arrow(length = unit(0.3, "cm")), color = "firebrick", size = 1) +
  # Labels
  geom_text_repel(data = loadings,
                  aes(x = PC1 * scale_factor, y = PC2 * scale_factor, label = var),
                  color = "black", size = 4, fontface = "bold", box.padding = 0.5) +
  labs(title = "Environmental Gradients (PC1 vs PC2)",
       subtitle = "Variable contributions to the primary axes",
       x = paste0("PC1 (", round(eig_val[1,"variance.percent"], 1), "%)"),
       y = paste0("PC2 (", round(eig_val[2,"variance.percent"], 1), "%)")) +
  coord_fixed() +
  theme_minimal(base_size = 14)

ggsave(file.path(OUTPUT_DIR, "02_Biplot_PC1_PC2.png"), p_biplot, width = 10, height = 8, bg = "white")

# --- 3. LOADINGS BAR PLOTS (FIXED LAYOUT) ---
cat("Generating Loadings Bar Plots (PC1-4)...\n")

rot <- as.data.frame(pca_model$rotation[, 1:4]) 
rot$Variable <- rownames(rot)
rot_long <- pivot_longer(rot, cols = starts_with("PC"), names_to = "PC", values_to = "Loading")

p_loadings <- ggplot(rot_long, aes(x = reorder(Variable, abs(Loading)), y = Loading, fill = Loading > 0)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black", linewidth=0.5) +
  scale_fill_manual(values = c("firebrick", "steelblue")) +
  
  # --- FIX: Use 2 columns instead of 4 to give bars room to breathe ---
  facet_wrap(~PC, ncol = 2) + 
  
  coord_flip() +
  labs(title = "Variable Loadings for Selected Components", 
       subtitle = "Interpretation of environmental drivers for PC1–PC4",
       x = NULL, y = "Loading Strength") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12),
        panel.border = element_rect(color="grey80", fill=NA),
        axis.text.y = element_text(size=10)) # Ensure labels are readable

# Increased height to accommodate the grid layout
ggsave(file.path(OUTPUT_DIR, "03_Loadings_BarPlot.png"), p_loadings, width = 12, height = 10, bg = "white")

cat("\nPlots saved to:", OUTPUT_DIR, "\n")