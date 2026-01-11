# scripts/99_env_justification.R
# ------------------------------------------------------------------------------
# FIGURE S2: ENVIRONMENTAL VARIABLE JUSTIFICATION
# - Scree Plot (PCA Variance)
# - Variable Contribution
# - VIF Check
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, ggplot2, factoextra, usdm, corrplot)

# --- CONFIG ---
BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
RAW_ENV_DIR <- file.path(DATA_DIR, "env", "current")
OUTPUT_PLOT <- file.path(BASE_DIR, "outputs", "figures", "env_justification")
dir.create(OUTPUT_PLOT, recursive = TRUE, showWarnings = FALSE)

CLIMATE_VARS <- c("sws_baseline_depthsurf_mean", "so_baseline_depthmax_mean", 
                  "thetao_baseline_depthmax_mean", "no3_baseline_depthmax_mean", 
                  "no3_baseline_depthmax_range", "chl_baseline_depthmax_mean", 
                  "o2_baseline_depthmax_range", "phyc_baseline_depthmax_mean")

# --- 1. LOAD DATA ---
cat("Loading Data...\n")
clim_list <- list()
for(v in CLIMATE_VARS) {
  r <- terra::rast(file.path(RAW_ENV_DIR, paste0(v, ".tif")))
  names(r) <- v
  clim_list[[v]] <- r
}
clim_stack <- terra::rast(clim_list)

# Sample for PCA (50k points)
set.seed(42)
samp <- terra::spatSample(clim_stack, size=50000, method="random", na.rm=TRUE, xy=FALSE)
samp <- samp[complete.cases(samp), ]

# --- 2. RUN PCA ---
cat("Running PCA...\n")
pca_res <- prcomp(samp, center=TRUE, scale.=TRUE)

# --- 3. JUSTIFICATION PLOTS ---

# A. Scree Plot (Variance Explained)
p1 <- fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50)) +
  labs(title = "Variance Explained by PC Axes", 
       subtitle = "Justification for selecting top 5 PCs (>95% Variance)") +
  theme_minimal()
ggsave(file.path(OUTPUT_PLOT, "Scree_Plot.png"), p1, width=8, height=6, bg="white")

# B. Variable Contributions (Circle Plot)
p2 <- fviz_pca_var(pca_res, col.var = "contrib", 
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                   repel = TRUE) +
  labs(title = "Variable Contributions (PC1 & PC2)") +
  theme_minimal()
ggsave(file.path(OUTPUT_PLOT, "Variable_Contributions.png"), p2, width=8, height=8, bg="white")

# --- 4. VIF ANALYSIS (Multi-collinearity Check) ---
cat("Calculating VIF for Final Predictors (PC1-5 + Rugosity)...\n")

# Load Rugosity (Static Variable)
rug_r <- terra::rast(file.path(DATA_DIR, "env", "terrain", "rugosity.tif"))

# Resample Rugosity to match Climate resolution if needed
if(!compareGeom(rug_r, clim_stack, stopOnError=FALSE)) rug_r <- resample(rug_r, clim_stack)

# Extract new sample with both Climate and Rugosity
full_stack <- c(clim_stack, rug_r)
vif_samp <- terra::spatSample(full_stack, size=10000, method="random", na.rm=TRUE, xy=FALSE)
vif_samp <- vif_samp[complete.cases(vif_samp), ]

# 1. Project Climate columns to PC1-5
clim_only <- vif_samp[, CLIMATE_VARS]
pca_scores <- predict(pca_res, clim_only)[, 1:5]

# 2. Combine with Rugosity
final_predictors <- data.frame(pca_scores, rugosity = vif_samp[, "rugosity"])

# 3. Run VIF
vif_res <- usdm::vif(final_predictors)
print(vif_res)

write.csv(vif_res, file.path(OUTPUT_PLOT, "Final_VIF_Scores.csv"), row.names=FALSE)

cat("--- FINISHED ---\n")
cat("Plots saved to:", OUTPUT_PLOT, "\n")