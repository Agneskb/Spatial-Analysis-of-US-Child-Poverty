# =============================================================================
# 03_modeling.R
# Purpose : OLS, spatial error, spatial lag, and spatial regime models
# Inputs  : data/processed/dat_spatial.rds
# Outputs : outputs/tables/ (stargazer HTML/tex), data/processed/model_results.rds
# =============================================================================

library(spatialreg)
library(spdep)
library(stargazer)
library(tidyverse)

spatial_obj <- readRDS("data/processed/dat_spatial.rds")
dat      <- spatial_obj$dat
data.wts <- spatial_obj$data.wts

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

formula_full <- childpov10 ~ singlehh + blkpc + hispc + unemploy +
  unemploymale + workin + extr + mnf + lesshigh + metrocut + logmedianincome

# ── OLS ───────────────────────────────────────────────────────────────────────
fit     <- lm(formula_full, data = dat)
lm.morantest(fit, listw = data.wts)   # test for residual autocorrelation

# ── Spatial error model ────────────────────────────────────────────────────────
fit.err <- errorsarlm(formula_full, data = dat, listw = data.wts)

# ── Spatial lag model ──────────────────────────────────────────────────────────
fit.lag <- lagsarlm(formula_full, data = dat, listw = data.wts, type = "lag")

# ── Spatial regime model (metro vs non-metro) ─────────────────────────────────
formula_regime <- childpov10 ~ singlehh + blkpc + hispc + unemploy +
  unemploymale + workin + extr + mnf + lesshigh + logmedianincome

# Metro subset
knn_m  <- knn2nb(knearneigh(coordinates(as(dat[dat$metrocut == "Metropolitan",],     "Spatial")), k = 4), sym = TRUE)
wts_m  <- nb2listw(knn_m)
efit.1 <- errorsarlm(formula_regime, data = dat[dat$metrocut == "Metropolitan",],     listw = wts_m, method = "MC")

# Non-metro subset
knn_n  <- knn2nb(knearneigh(coordinates(as(dat[dat$metrocut == "Non-Metropolitan",], "Spatial")), k = 4), sym = TRUE)
wts_n  <- nb2listw(knn_n)
efit.2 <- errorsarlm(formula_regime, data = dat[dat$metrocut == "Non-Metropolitan",], listw = wts_n, method = "MC")

# ── Export tables ─────────────────────────────────────────────────────────────
stargazer(fit, fit.err, fit.lag,
          type    = "html",
          out     = "outputs/tables/table2_model_comparison.html",
          title   = "Table 2: OLS vs Spatial Error vs Spatial Lag",
          column.labels = c("OLS","Spatial Error","Spatial Lag"))

stargazer(efit.1, efit.2,
          type  = "html",
          out   = "outputs/tables/table3_regime.html",
          title = "Table 3: Spatial Regime — Metro vs Non-Metro")

# ── Save results ───────────────────────────────────────────────────────────────
saveRDS(list(ols = fit, err = fit.err, lag = fit.lag,
             metro = efit.1, nonmetro = efit.2),
        "data/processed/model_results.rds")

