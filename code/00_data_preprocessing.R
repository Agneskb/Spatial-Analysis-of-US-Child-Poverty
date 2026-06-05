# =============================================================================
# 00_data_preprocessing.R
# Purpose : Load AHRF and ACS data, clean, merge at county level
# Outputs : data/processed/county_merged.rds
# =============================================================================

library(haven)
library(tidycensus)
library(tidyverse)
library(sf)
library(tigris)

# ── Census API key ─────────────────────────────────────────────────────────────
# Run ONCE interactively: census_api_key("YOUR_KEY", install = TRUE)

# ── Load AHRF ──────────────────────────────────────────────────────────────────
ahrf_path <- Sys.getenv("AHRF_PATH")   # set in .Renviron
arf2020   <- read_sas(ahrf_path)

popdat <- arf2020 |>
  mutate(
    cofips       = f00004,
    coname       = f00010,
    state        = f00011,
    childpov10   = f1332219,
    rucc         = as.factor(f0002013),
    singlehh     = 100 * f1160310 / f0453010,
    blkpc        = f0453810,
    hispc        = f0454210,
    unemploy     = f1451215 / f1451015 * 100,
    unemploymale = (f1451215 - f1451515) / (f1451015 - f1451315) * 100,
    workin       = f1460615,
    extr         = f1462115,
    mnf          = f1458715,
    medianincome = f1434615
  ) |>
  select(state, cofips, coname, childpov10, rucc,
         singlehh, blkpc, hispc, unemploy, unemploymale,
         workin, extr, mnf, medianincome) |>
  filter(complete.cases(.)) |>
  filter(!state %in% c("02","15","60","66","69","72","78")) |>
  as.data.frame()

# ── Metro / non-metro classification ──────────────────────────────────────────
popdat$metrocut <- cut(
  as.numeric(popdat$rucc),
  breaks         = c(0, 6, 10),
  include.lowest = TRUE,
  labels         = c("Metropolitan", "Non-Metropolitan")
)

# ── ACS educational attainment (table B06009, 2015) ───────────────────────────
edutable <- get_acs(
  geography   = "county",
  year        = 2015,
  geometry    = TRUE,
  output      = "wide",
  table       = "B06009",
  cache_table = TRUE
)

edudat <- edutable |>
  mutate(
    lesshigh = (B06009_002E + B06009_003E) / B06009_001E * 100,
    cofips   = substr(GEOID, 1, 5)
  ) |>
  select(lesshigh, cofips)

# ── Spatial join ───────────────────────────────────────────────────────────────
dat <- geo_join(
  spatial_data = edudat,
  data_frame   = popdat,
  by_sp        = "cofips",
  by_df        = "cofips",
  how          = "inner"
) |>
  mutate(logmedianincome = log(medianincome))

# ── Save ───────────────────────────────────────────────────────────────────────
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(dat, "data/processed/county_merged.rds")


