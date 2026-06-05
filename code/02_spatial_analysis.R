# =============================================================================
# 02_spatial_analysis.R
# Purpose : Spatial autocorrelation tests and LISA cluster map
# Inputs  : data/processed/county_merged.rds
# Outputs : data/processed/dat_spatial.rds, figures/moran_lisa_map.png
# =============================================================================

library(spdep)
library(sf)
library(sp)
library(tidyverse)

dat <- readRDS("data/processed/county_merged.rds")

# ── Global Moran's I ──────────────────────────────────────────────────────────
knn      <- knearneigh(coordinates(as(dat, "Spatial")), k = 4)
knn_nb   <- knn2nb(knn)
data.wts <- nb2listw(knn_nb, style = "W")

moran_global <- moran.test(
  x     = dat$childpov10[!is.na(dat$childpov10)],
  listw = data.wts
)
print(moran_global)

# ── Local Moran's I / LISA ────────────────────────────────────────────────────
local_moran <- localmoran(
  x           = dat$childpov10[!is.na(dat$childpov10)],
  listw       = data.wts,
  alternative = "two.sided"
)

data_sp          <- as(dat, "Spatial")
data_sp$local_I  <- local_moran[, 1]
data_sp$local_p  <- local_moran[, 5]

# Quadrant classification
sinc <- scale(data_sp$childpov10)
inc  <- lag.listw(var = sinc, x = data.wts)

data_sp$quadsig <- NA_integer_
data_sp$quadsig[sinc >  0 & inc >  0 & data_sp$local_p < 0.05] <- 1L  # High-High
data_sp$quadsig[sinc <  0 & inc <  0 & data_sp$local_p < 0.05] <- 2L  # Low-Low
data_sp$quadsig[sinc >  0 & inc <  0 & data_sp$local_p < 0.05] <- 3L  # High-Low
data_sp$quadsig[sinc <  0 & inc >  0 & data_sp$local_p < 0.05] <- 4L  # Low-High
data_sp$quadsig[data_sp$local_p >= 0.05]                        <- 5L  # Not significant

labels <- c("High-High","Low-Low","High-Low","Low-High","Not Significant")
colors <- c("red","blue","lightpink","lightblue","white")
np     <- findInterval(data_sp$quadsig, seq(1, 5, 1))

png("figures/moran_lisa_map.png", width = 1200, height = 700, res = 130)
plot(data_sp, col = colors[np], main = "Moran LISA Cluster Map — Child Poverty Rate")
legend("bottomleft", legend = labels, fill = colors, bty = "n", cex = 0.8)
dev.off()

# ── Save spatial object ───────────────────────────────────────────────────────
saveRDS(list(dat = dat, data.wts = data.wts), "data/processed/dat_spatial.rds")

