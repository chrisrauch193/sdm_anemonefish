# scripts/99_env_justification.R
# ------------------------------------------------------------------------------
# STATISTICAL JUSTIFICATION: Variance & VIF Analysis
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, usdm, factoextra)

# --- CONFIG ---
BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
RAW_ENV_DIR <- file.path(DATA_DIR, "env", "current")
OUTPUT_DIR  <- file.path(BASE_DIR, "outputs", "figures", "env_justification")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PCA_MODEL_PATH <- file.path(DATA_DIR, "env", "pca_model.rds")
RUG_PATH       <- file.path(DATA_DIR, "env", "terrain", "rugosity.tif")

# Variables (Must match what PCA was trained on)
CLIMATE_VARS <- c("sws_baseline_depthsurf_mean", "so_baseline_depthmax_mean", 
                  "thetao_baseline_depthmax_mean", "no3_baseline_depthmax_mean", 
                  "no3_baseline_depthmax_range", "chl_baseline_depthmax_mean", 
                  "o2_baseline_depthmax_range", "phyc_baseline_depthmax_mean")

# --- EXECUTE ---
cat("Loading Data...\n")
if(!file.exists(PCA_MODEL_PATH)) stop("PCA Model not found.")
pca_model <- readRDS(PCA_MODEL_PATH)

# --- 1. KAISER-GUTTMAN CHECK (The Proof) ---
cat("\n======================================================\n")
cat("          PCA SELECTION JUSTIFICATION \n")
cat("======================================================\n")

eig_val <- get_eigenvalue(pca_model)
print(round(eig_val, 3))

cat("\n[DECISION LOGIC]\n")
cat("Kaiser Criterion (Eigenvalue > 1):\n")
pcs_keep <- rownames(eig_val)[eig_val$eigenvalue > 1]
cat(" -> Retain:", paste(pcs_keep, collapse=", "), "\n")

# Check PC4/PC5 specifically
pc4_cum <- eig_val["Dim.4", "cumulative.variance.percent"]
pc5_val <- eig_val["Dim.5", "eigenvalue"]

cat(paste0("\nPC4 Cumulative Variance: ", round(pc4_cum, 2), "%\n"))
cat(paste0("PC5 Eigenvalue: ", round(pc5_val, 3), " (Should be < 1 to drop)\n"))

if(pc5_val < 1) {
  cat("\nVERDICT: STOP at PC4. PC5 is noise.\n")
} else {
  cat("\nVERDICT: Consider keeping PC5 (Eigenvalue > 1).\n")
}
cat("======================================================\n")

# --- 2. PREPARE DATA FOR VIF ---
cat("\nPreparing data for VIF check on selected variables (PC1-4 + Rugosity)...\n")

# Load Raw Climate (to project) + Rugosity
clim_list <- list()
for(v in CLIMATE_VARS) clim_list[[v]] <- terra::rast(file.path(RAW_ENV_DIR, paste0(v, ".tif")))
clim_stack <- terra::rast(clim_list)

rug_r <- terra::rast(RUG_PATH)
if(!terra::compareGeom(rug_r, clim_stack, stopOnError=FALSE)) rug_r <- terra::resample(rug_r, clim_stack)

# Sample points for VIF calculation
set.seed(42)
full_stack <- c(clim_stack, rug_r)
samp <- terra::spatSample(full_stack, size=10000, method="random", na.rm=TRUE, xy=FALSE)
samp <- samp[complete.cases(samp), ]

# --- 3. CALCULATE VIF ---
cat("Projecting PC1-4...\n")
# 1. Predict PC1-4
clim_only <- samp[, CLIMATE_VARS]
pca_scores <- predict(pca_model, clim_only)[, 1:4] # EXPLICITLY 1:4

# 2. Combine with Rugosity
final_df <- data.frame(pca_scores, rugosity = samp[, "rugosity"])

# 3. Calculate VIF
cat("Calculating VIF...\n")
vif_res <- usdm::vif(final_df)

print(vif_res)
write.csv(vif_res, file.path(OUTPUT_DIR, "Final_VIF_Table.csv"), row.names=FALSE)

cat("\n--- CHECK COMPLETE ---\n")
if(max(vif_res$VIF) < 5) {
  cat("SUCCESS: All VIF values < 5. Predictors are independent.\n")
} else {
  cat("WARNING: High VIF detected.\n")
}