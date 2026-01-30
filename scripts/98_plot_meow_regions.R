# scripts/98_plot_meow_regions.R
# ------------------------------------------------------------------------------
# FIGURE 1: THESIS STUDY AREA (MEOW PROVINCES) - FINAL POLISH
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(sf, ggplot2, dplyr, readr, ggrepel, viridis, patchwork, stringr, grid, maps)

# --- CONFIG ---
BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
OUTPUT_DIR  <- file.path(BASE_DIR, "outputs", "figures")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

REGION_FILE <- file.path(DATA_DIR, "meow_ecos_df.csv")
if(!file.exists(REGION_FILE)) stop("Missing meow_ecos_df.csv")

# --- 1. DATA PREP ---
cat("Loading Region Data...\n")
df_raw <- read_csv(REGION_FILE, show_col_types = FALSE)
colnames(df_raw) <- tolower(colnames(df_raw))

# Define Order (West to East)
realm_levels <- c(
  "Western Indo-Pacific", 
  "Central Indo-Pacific", 
  "Eastern Indo-Pacific", 
  "Temperate Northern Pacific",
  "Temperate Australasia", 
  "Temperate Southern Africa"
)

# Standardize
df <- df_raw %>% 
  rename(long = long, lat = lat, PROVINCE = province, REALM = realm, group = id) %>%
  mutate(long = ifelse(long < 0, long + 360, long)) %>%
  mutate(REALM = factor(REALM, levels = realm_levels))

# Calculate Centroids & Sort
centroids <- df %>% 
  group_by(PROVINCE, REALM) %>% 
  summarise(X_cent = mean(long), Y_cent = mean(lat), .groups="drop") %>% 
  arrange(REALM, X_cent) %>% 
  mutate(
    Map_Index = row_number(),
    Clean_Name = str_replace(PROVINCE, " Shelf", ""), 
    # Increased truncation length to fit long names
    Clean_Name = str_trunc(Clean_Name, 35, "right") 
  )

df_final <- df %>% left_join(centroids %>% dplyr::select(PROVINCE, Map_Index), by="PROVINCE")

# Custom Thesis Palette
realm_cols <- c(
  "Western Indo-Pacific"        = "#E6AB02", # Gold
  "Central Indo-Pacific"        = "#66A61E", # Green
  "Eastern Indo-Pacific"        = "#D95F02", # Orange
  "Temperate Northern Pacific" = "#7570B3", # Purple
  "Temperate Australasia"       = "#E7298A", # Pink
  "Temperate Southern Africa"   = "#A6761D"  # Brown
)

# --- 2. COMPONENT 1: THE MAP (Top) ---
# NOTE: Plot Order Matters! Regions first (bottom), then Land (top).

p_map <- ggplot() +
  # LAYER 1: Marine Regions (Bottom)
  geom_polygon(data = df_final, aes(x = long, y = lat, group = group, fill = REALM),
               color = "white", linewidth = 0.05, alpha = 0.85) +
  
  # LAYER 2: Land Mass (Top) - "world2" wraps 0-360 correctly
  borders("world2", colour = "gray60", fill = "gray90", size = 0.1) +
  
  # LAYER 3: Labels
  geom_label_repel(data = centroids, aes(x = X_cent, y = Y_cent, label = Map_Index),
                   size = 3, fontface = "bold", 
                   label.padding = unit(0.1, "lines"), label.r = unit(0.1, "lines"),
                   fill = "white", alpha = 0.9, segment.size = 0.2,
                   min.segment.length = 0, box.padding = 0.2) +
  
  scale_fill_manual(values = realm_cols) +
  
  # Use quickmap for better aspect ratio preservation on lat/long data
  coord_quickmap(xlim = c(25, 260), ylim = c(-45, 45), expand = FALSE) +
  
  theme_void() + 
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(t=5, b=2)),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size=10, margin = margin(b=2)),
    panel.background = element_rect(fill = "#E0F3F8", color = "black", linewidth = 0.5), # Ocean color
    plot.margin = margin(5, 5, 2, 5) 
  )

# --- 3. COMPONENT 2: REALM LEGEND (Bottom Left) ---
realm_legend_df <- data.frame(REALM = factor(realm_levels, levels = realm_levels))

p_realm_legend <- ggplot(realm_legend_df, aes(x=1, y=REALM, color=REALM)) +
  geom_point(size = 4) +
  geom_text(aes(label = REALM), hjust = 0, nudge_x = 0.2, size = 3, fontface = "bold", color="black") +
  scale_color_manual(values = realm_cols) +
  scale_y_discrete(limits = rev(realm_levels)) + 
  xlim(1, 4) + 
  labs(title = "Realms") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold", size=11, hjust=0, margin=margin(b=5)),
    plot.margin = margin(10, 0, 10, 10) 
  )

# --- 4. COMPONENT 3: PROVINCE LIST (Bottom Right) ---
n_cols <- 4
centroids$col_id <- ceiling(seq_len(nrow(centroids)) / (nrow(centroids)/n_cols))

p_prov_list <- ggplot(centroids, aes(x = 0, y = -Map_Index)) +
  geom_text(aes(label = paste0(Map_Index, ". ", Clean_Name), color = REALM), 
            hjust = 0, size = 2.5, fontface="plain") + 
  scale_color_manual(values = realm_cols) +
  facet_wrap(~col_id, scales = "free_y", ncol = n_cols) +
  labs(title = "Provinces") +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    plot.title = element_text(face="bold", size=11, hjust=0, margin=margin(b=5)),
    plot.margin = margin(10, 10, 10, 0) 
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1))

# --- 5. ASSEMBLE ---
bottom_row <- p_realm_legend + p_prov_list + plot_layout(widths = c(1, 2.8))

final_plot <- p_map / bottom_row + plot_layout(heights = c(3, 1.2))

# SAVE
out_file <- file.path(OUTPUT_DIR, "Fig1_Study_Area_Map_Final.png")
ggsave(out_file, final_plot, width = 12, height = 9, dpi = 300, bg="white")

cat("Success! Saved final composite map to:", out_file, "\n")