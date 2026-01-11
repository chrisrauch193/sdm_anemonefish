# scripts/98_plot_meow_regions.R
# ------------------------------------------------------------------------------
# FIGURE 1: THESIS STUDY AREA (MEOW PROVINCES)
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(sf, ggplot2, dplyr, readr, ggrepel, viridis, patchwork, stringr)

# --- CONFIG ---
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "outputs", "figures")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

# INPUT: Your finalized CSV from 0a_marine_regions.R
REGION_FILE <- file.path(DATA_DIR, "meow_ecos_df.csv")

if(!file.exists(REGION_FILE)) stop("Missing meow_ecos_df.csv. Run 0a_marine_regions.R first.")

# --- 1. DATA PREP ---
cat("Loading Region Data...\n")
df_raw <- read_csv(REGION_FILE, show_col_types = FALSE)

# Standardize Columns
colnames(df_raw) <- tolower(colnames(df_raw))
df <- df_raw %>% 
  rename(long = long, lat = lat, PROVINCE = province, REALM = realm, group = id)

# Shift 0-360 if needed
if(min(df$long, na.rm=TRUE) < 0) df <- df %>% mutate(long = ifelse(long < 0, long + 360, long))

# Calculate Centroids for Labels
centroids <- df %>% 
  group_by(PROVINCE, REALM) %>% 
  summarise(X_cent = mean(long), Y_cent = mean(lat), .groups="drop") %>% 
  arrange(X_cent) %>% 
  mutate(
    Map_Index = row_number(),
    Clean_Name = str_trunc(PROVINCE, 25, "right")
  )

# Join Index back to Polygons
df_final <- df %>% left_join(centroids %>% dplyr::select(PROVINCE, Map_Index), by="PROVINCE")

# --- 2. PLOTTING ---
cat("Generating Thesis Study Area Map...\n")

# Color Palette (Realms)
realm_cols <- c(
  "Central Indo-Pacific" = "#E6AB02", 
  "Eastern Indo-Pacific" = "#A6761D",
  "Western Indo-Pacific" = "#66A61E", 
  "Temperate Northern Pacific" = "#7570B3",
  "Temperate Australasia" = "#E7298A", 
  "Temperate Southern Africa" = "#D95F02"
)

# A. MAP
p_map <- ggplot() +
  # Background World (0-360)
  borders("world2", colour="gray90", fill="gray95", size=0.1) +
  
  # The Regions
  geom_polygon(data = df_final, aes(x = long, y = lat, group = group, fill = REALM),
               color = "white", size = 0.1, alpha = 0.9) +
  
  # Labels
  geom_text_repel(data = centroids, aes(x = X_cent, y = Y_cent, label = Map_Index),
                  size = 3, fontface = "bold", min.segment.length = 0,
                  box.padding = 0.4, segment.size = 0.3, segment.color = "gray40") +
  
  scale_fill_manual(values = realm_cols) +
  
  # Zoom: South Africa to Easter Island Buffer, Japan to NZ
  coord_fixed(xlim = c(20, 260), ylim = c(-50, 50), expand = FALSE) +
  
  labs(title = "Study Area: Marine Ecoregions of the World (MEOW)", 
       subtitle = paste("Selected Provinces:", nrow(centroids))) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color="gray40"),
    panel.background = element_rect(fill = "#E0F3F8", color = NA),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

# B. KEY
n_rows <- ceiling(nrow(centroids) / 4)
centroids$col_group <- rep(1:4, each = n_rows, length.out = nrow(centroids))

p_key <- ggplot(centroids, aes(x = 0, y = -Map_Index, label = paste0(Map_Index, ". ", Clean_Name))) +
  geom_text(aes(color = REALM), hjust = 0, size = 3, fontface = "bold") +
  scale_color_manual(values = realm_cols) +
  facet_wrap(~col_group, scales = "free") +
  theme_void() +
  theme(legend.position = "none", strip.text = element_blank()) +
  scale_x_continuous(expand = c(0, 0.1))

# C. COMBINE
final_plot <- p_map / p_key + plot_layout(heights = c(4, 1))

# SAVE
out_file <- file.path(OUTPUT_DIR, "Fig1_Study_Area_Map.png")
ggsave(out_file, final_plot, width = 12, height = 10, dpi = 300, bg="white")
cat("Saved:", out_file, "\n")