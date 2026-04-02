# ============================================================
# Climate profiles and map for three parental genotypes
# Repository root: clover_drought
#
# External files not included in GitHub:
# - CLOVER_DROUGHT_ENV_DIR/                  -> Please download from DRYAD (10.5061/dryad.6djh9w1hc)

#
# Before running:
# Set the environment variable CLOVER_DROUGHT_ENV_DIR to the folder
# containing the external environment files.
# Example:
Sys.setenv(CLOVER_DROUGHT_ENV_DIR = "/Users/kuowenhsi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Drought_F3_paper/Env")
setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Drought_F3_paper/clover_drought")
# ============================================================

library(tidyverse)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(scatterpie)
library(ggpubr)
library(raster)

# ----------------------------
# Paths
# ----------------------------

repo_data_file <- file.path("data", "cyano_buf_climate_data_419.csv")
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

external_env_dir <- Sys.getenv("CLOVER_DROUGHT_ENV_DIR")

if (external_env_dir == "") {
  stop(
    paste(
      "Please set CLOVER_DROUGHT_ENV_DIR to the folder containing external",
      "geospatial environment data. These files are not included in GitHub.",
      "Please download from Internet."
    )
  )
}

tavg_dir <- file.path(external_env_dir, "wc2.1_2.5m_tavg")
prec_dir <- file.path(external_env_dir, "wc2.1_2.5m_prec")
def_file <- file.path(external_env_dir, "Genotype_monthly_DEF.csv")
gdd_raster_file <- file.path(external_env_dir, "current_30arcsec_growingDegDays5.tif")

# ----------------------------
# Input checks
# ----------------------------

if (!file.exists(repo_data_file)) {
  stop("Missing file in repository: data/cyano_buf_climate_data_419.csv")
}

if (!dir.exists(tavg_dir)) {
  stop("Missing external folder: wc2.1_2.5m_tavg. Please download from Internet.")
}

if (!dir.exists(prec_dir)) {
  stop("Missing external folder: wc2.1_2.5m_prec. Please download from Internet.")
}

if (!file.exists(def_file)) {
  stop("Missing external file: Genotype_monthly_DEF.csv. Please download from Internet.")
}

if (!file.exists(gdd_raster_file)) {
  stop("Missing external file: current_30arcsec_growingDegDays5.tif. Please download from Internet.")
}

# ----------------------------
# Parent locations
# ----------------------------

parent_loc <- read_csv(repo_data_file) %>%
  dplyr::select(genotype, Latitude, Longitude) %>%
  filter(genotype %in% c("DMN_010", "STL_0701", "GFL_007")) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), agr = "constant", crs = 4326)

usa_state <- ne_states(
  country = "United States of America",
  returnclass = "sf"
) %>%
  st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60) %>%
  st_geometry()

canada_state <- ne_states(
  country = "Canada",
  returnclass = "sf"
) %>%
  st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60) %>%
  st_geometry()

mexico_state <- ne_states(
  country = "Mexico",
  returnclass = "sf"
) %>%
  st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60) %>%
  st_geometry()

p <- ggplot(data = parent_loc) +
  geom_sf(data = usa_state, fill = NA, color = "gray75") +
  geom_sf(data = canada_state, fill = NA, color = "gray75") +
  geom_sf(data = mexico_state, fill = NA, color = "gray75") +
  geom_sf(aes(color = genotype, shape = genotype), size = 2) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid = element_blank()
  ) +
  coord_sf(ylim = c(2e6, 7e6), expand = FALSE, crs = 3857)

p

ggsave(
  filename = file.path(output_dir, "just_parents_location.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 600
)

# ----------------------------
# Monthly temperature
# ----------------------------

crop_extent_large <- extent(-135, -55, 15, 60)

tmp_files <- list.files(tavg_dir, full.names = TRUE)
tmp_layer <- stack(tmp_files) %>%
  raster::crop(crop_extent_large)

plot(tmp_layer)

tmp_extracted <- raster::extract(tmp_layer, parent_loc, df = TRUE) %>%
  mutate(ID = parent_loc$genotype) %>%
  pivot_longer(cols = 2:13, names_to = "month", values_to = "avg_temp") %>%
  mutate(month = as.integer(str_remove(month, "wc2.1_2.5m_tavg_"))) %>%
  group_by(ID) %>%
  mutate(GDD5 = avg_temp - 5) %>%
  mutate(GDD5 = ifelse(GDD5 > 0, GDD5 * 30, 0)) %>%
  mutate(GDD5 = cumsum(GDD5)) %>%
  ungroup()

# ----------------------------
# Monthly precipitation
# ----------------------------

prec_files <- list.files(prec_dir, full.names = TRUE)
prep_layer <- stack(prec_files) %>%
  raster::crop(crop_extent_large)

plot(prep_layer)

prep_extracted <- raster::extract(prep_layer, parent_loc, df = TRUE) %>%
  mutate(ID = parent_loc$genotype) %>%
  pivot_longer(cols = 2:13, names_to = "month", values_to = "avg_prep") %>%
  mutate(month = as.integer(str_remove(month, "wc2.1_2.5m_prec_")))

# ----------------------------
# DEF data
# ----------------------------

def_extracted <- read_csv(def_file) %>%
  rename(ID = genotype)

# ----------------------------
# Combined environment table
# ----------------------------

comb_extracted <- left_join(tmp_extracted, prep_extracted, by = c("ID", "month")) %>%
  left_join(def_extracted, by = c("ID", "month")) %>%
  mutate(
    avg_prep = avg_prep / 10,
    def_mean = def_mean / 10
  ) %>%
  pivot_longer(cols = 3:6, names_to = "env_varibles", values_to = "values") %>%
  filter(env_varibles != "GDD5")

# ----------------------------
# GDD plot
# ----------------------------

pGDD <- ggplot(tmp_extracted, aes(x = month, y = GDD5 * 100)) +
  geom_hline(yintercept = 1000, linetype = "dashed") +
  geom_point(aes(color = ID), size = 2) +
  geom_line(aes(color = ID), alpha = 1) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(
    expression("Cumulative Growing Degree Days (> 5" * degree * "C)"),
    n.breaks = 10
  ) +
  scale_color_manual(
    "Locations",
    values = c("blue", "red2", "green4"),
    labels = c(
      DMN_010 = "Duluth, MN",
      GFL_007 = "Gainsville, FL",
      STL_0701 = "St. Louis, MO"
    )
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pGDD

ggsave(
  filename = file.path(output_dir, "GDD5_monthly_20240530.png"),
  plot = pGDD,
  width = 8,
  height = 4,
  dpi = 600
)

# ----------------------------
# Legend plot
# ----------------------------

p_legend <- ggplot(
  data = filter(comb_extracted, ID == "DMN_010"),
  aes(x = month, y = values)
) +
  geom_point(aes(color = env_varibles)) +
  geom_line(aes(color = env_varibles), alpha = 0.7) +
  scale_y_continuous(
    "Temperature (°C)",
    sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"),
    limits = c(-15, 58)
  ) +
  scale_x_continuous(breaks = c(1:6) * 2) +
  scale_color_manual(
    name = "",
    values = c("avg_prep" = "#0095fd", "avg_temp" = "#fd9600", "def_mean" = "#3f4f18"),
    labels = c(
      "avg_prep" = "avg. prep.",
      "avg_temp" = "avg. temp.",
      "def_mean" = "water deficit"
    )
  ) +
  labs(title = "DMN_010") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(color = "#F8766D", hjust = 0.5, face = "bold", size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    plot.background = element_rect(fill = "#FFFFFF90", color = NA),
    legend.background = element_rect(fill = NA, color = NA)
  )

p_legend

# ----------------------------
# Inset plots
# ----------------------------

p_DMN010 <- ggplot(
  data = filter(comb_extracted, ID == "DMN_010", env_varibles != "def_mean"),
  aes(x = month, y = values)
) +
  geom_col(
    data = filter(comb_extracted, ID == "DMN_010", env_varibles == "def_mean"),
    fill = "red3"
  ) +
  geom_point(aes(color = env_varibles)) +
  geom_line(aes(color = env_varibles), alpha = 0.7) +
  scale_y_continuous(
    "Temperature (°C)",
    sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"),
    limits = c(-15, 20)
  ) +
  scale_x_continuous(breaks = c(1:6) * 2) +
  scale_color_manual(
    name = "",
    values = c("avg_prep" = "#0095fd", "avg_temp" = "#fd9600"),
    labels = c("avg_prep" = "avg. prep.", "avg_temp" = "avg. temp.")
  ) +
  labs(title = "DMN_010") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#EE000030"),
    plot.title = element_text(color = "#EE000090", hjust = 0.5, face = "bold", size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    plot.background = element_rect(fill = NA, color = NA)
  )

p_DMN010

ggsave(
  filename = file.path(output_dir, "DMN_monthly.png"),
  plot = p_DMN010,
  width = 3.5,
  height = 2,
  dpi = 600
)

p_STL0701 <- ggplot(
  data = filter(comb_extracted, ID == "STL_0701", env_varibles != "def_mean"),
  aes(x = month, y = values)
) +
  geom_col(
    data = filter(comb_extracted, ID == "STL_0701", env_varibles == "def_mean"),
    fill = "red3"
  ) +
  geom_point(aes(color = env_varibles)) +
  geom_line(aes(color = env_varibles), alpha = 0.7) +
  scale_y_continuous(
    "Temperature (°C)",
    sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"),
    limits = c(-5, 30)
  ) +
  scale_x_continuous(breaks = c(1:6) * 2) +
  scale_color_manual(
    name = "",
    values = c("avg_prep" = "#0095fd", "avg_temp" = "#fd9600"),
    labels = c("avg_prep" = "avg. prep.", "avg_temp" = "avg. temp.")
  ) +
  labs(title = "STL_0701") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#0000FF30"),
    plot.title = element_text(color = "#0000FF90", hjust = 0.5, face = "bold", size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    plot.background = element_rect(fill = NA, color = NA)
  )

p_STL0701

ggsave(
  filename = file.path(output_dir, "STL_monthly.png"),
  plot = p_STL0701,
  width = 3.5,
  height = 2,
  dpi = 600
)

p_GFL007 <- ggplot(
  data = filter(comb_extracted, ID == "GFL_007", env_varibles != "def_mean"),
  aes(x = month, y = values)
) +
  geom_col(
    data = filter(comb_extracted, ID == "GFL_007", env_varibles == "def_mean"),
    fill = "red3"
  ) +
  geom_point(aes(color = env_varibles)) +
  geom_line(aes(color = env_varibles), alpha = 0.7) +
  scale_y_continuous(
    "Temperature (°C)",
    sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"),
    limits = c(0, 35)
  ) +
  scale_x_continuous(breaks = c(1:6) * 2) +
  scale_color_manual(
    name = "",
    values = c("avg_prep" = "#0095fd", "avg_temp" = "#fd9600"),
    labels = c("avg_prep" = "avg. prep.", "avg_temp" = "avg. temp.")
  ) +
  labs(title = "GFL_007") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#008B0030"),
    plot.title = element_text(color = "#008B0090", hjust = 0.5, face = "bold", size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    plot.background = element_rect(fill = NA, color = NA)
  )

p_GFL007

ggsave(
  filename = file.path(output_dir, "GFL_monthly.png"),
  plot = p_GFL007,
  width = 3.5,
  height = 2,
  dpi = 600
)

# ----------------------------
# Optional map with inset climate plots
# ----------------------------

map_with_insets <- p +
  annotation_custom(ggplotGrob(p_DMN010), xmin = -10e6, xmax = -7.7e6, ymin = 5e6, ymax = 7e6) +
  annotation_custom(ggplotGrob(p_STL0701), xmin = -12.8e6, xmax = -10.5e6, ymin = 3.8e6, ymax = 5.8e6) +
  annotation_custom(ggplotGrob(p_GFL007), xmin = -11.5e6, xmax = -9.2e6, ymin = 2e6, ymax = 4e6)

map_with_insets

# ----------------------------
# Growing degree day raster map
# ----------------------------

env_layer <- raster(gdd_raster_file)

env_layer_croped <- crop(env_layer, extent(-135, -65, 25, 50)) / 10

env_layer_df <- as.data.frame(env_layer_croped, xy = TRUE) %>%
  drop_na()

head(env_layer_df)
max(env_layer_croped)

p_env <- ggplot(data = parent_loc) +
  geom_raster(
    data = env_layer_df,
    aes(x = x, y = y, fill = current_30arcsec_growingDegDays5)
  ) +
  geom_sf(data = usa_state, fill = NA) +
  geom_sf(data = canada_state, fill = NA) +
  geom_sf(data = mexico_state, fill = NA) +
  geom_sf(aes(color = genotype, shape = genotype), size = 2, show.legend = FALSE) +
  scale_fill_viridis_c(option = "inferno", name = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.key.height = grid::unit(0.1, "in"),
    legend.key.width = grid::unit(0.4, "in"),
    legend.margin = margin(t = -0.02, b = -0.02, unit = "in"),
    legend.position = "top",
    panel.background = element_rect(fill = "white"),
    axis.title = element_text(size = 8),
    legend.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  ) +
  xlab(expression("Growing Degree Days 5 " * degree * "C")) +
  ylab("") +
  coord_sf(
    xlim = c(-130, -65),
    ylim = c(25, 50),
    expand = FALSE,
    label_axes = ""
  )

p_env

# ----------------------------
# Final composite figure
# ----------------------------

out_p <- p +
  annotation_custom(ggplotGrob(p_DMN010), xmin = -10e6, xmax = -7.7e6, ymin = 5e6, ymax = 7e6) +
  annotation_custom(ggplotGrob(p_STL0701), xmin = -12.8e6, xmax = -10.5e6, ymin = 3.8e6, ymax = 5.8e6) +
  annotation_custom(ggplotGrob(p_GFL007), xmin = -11.5e6, xmax = -9.2e6, ymin = 2e6, ymax = 4e6) +
  annotation_custom(ggplotGrob(p_env), xmin = -16e6, ymin = 2e6, xmax = -11e6, ymax = 4e6)

out_p

ggsave(
  filename = file.path(output_dir, "Figure_1_20260401.png"),
  plot = out_p,
  width = 8,
  height = 6,
  dpi = 600
)
