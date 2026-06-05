# =============================================================================
# 01_eda_maps.R
# Purpose : Exploratory maps of child poverty and county classification
# Inputs  : data/processed/county_merged.rds
# Outputs : figures/ (PNG maps)
# =============================================================================

library(tidyverse)
library(sf)
library(ggplot2)

dat <- readRDS("data/processed/county_merged.rds")
dir.create("figures", showWarnings = FALSE)

# ── Map 1: Continuous child poverty rate ──────────────────────────────────────
ggplot(dat) +
  geom_sf(aes(fill = childpov10, color = childpov10)) +
  scale_fill_viridis_c(name = "Child Poverty %") +
  scale_color_viridis_c(guide = "none") +
  labs(title = "Spatial Variation of County-Level Child Poverty Rate, 2019") +
  theme_minimal(base_size = 11)

ggsave("figures/map_childpov_continuous.png", width = 9, height = 5.5, dpi = 150)

# ── Map 2: Poverty rate categories ────────────────────────────────────────────
dat$povcut <- cut(
  as.numeric(dat$childpov10),
  breaks         = c(0, 15, 30, 100),
  include.lowest = TRUE,
  labels         = c("< 15%", "15–30%", "> 30%")
)

ggplot(dat) +
  geom_sf(aes(fill = factor(povcut)), lwd = 0.1) +
  scale_fill_manual(values = c("white", "grey70", "#1f78b4"), name = "Child Poverty") +
  labs(title = "County-Level Child Poverty Rate, 2019 (Categorized)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("figures/map_childpov_categories.png", width = 9, height = 5.5, dpi = 150)

# ── Map 3: Metro vs non-metro ─────────────────────────────────────────────────
ggplot(dat) +
  geom_sf(aes(fill = factor(metrocut)), lwd = 0.1) +
  scale_fill_manual(values = c("lightblue", "white"), name = "County Type") +
  labs(title = "Metropolitan vs Non-Metropolitan Counties") +
  theme_minimal(base_size = 11)

ggsave("figures/map_metro_nonmetro.png", width = 9, height = 5.5, dpi = 150)
