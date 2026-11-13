# Checks that the beta-binomial is appropriate for the index models fit in
# 01-index-generation.R. Reads the saved fits; run 01 first.
#
# 1. DHARMa QQ residuals from simulations (uniform residuals = adequate).
# 2. DHARMa dispersion test.
# 3. Whether dispersion on the probability scale matters for the index: the
#    index sums exp(eta) = -log(1 - p_bar), but the mean per-hook rate is
#    E(-log(1 - p)) = digamma(phi) - digamma((1 - p_bar) * phi), which isn't
#    proportional to exp(eta) in general. If the ratio of the two, summed over
#    the grid, is nearly constant across years, the relative index is nearly
#    the same either way.
library(ggplot2)
library(tidyverse)
library(sdmTMB)

models <- c("hblldog", "hbll-n-s", "dog-j-circle")

# prediction grids, as in 01-index-generation.R
grid_dog <- gfplot::dogfish_grid$grid %>%
  mutate(UTM.lon = X / 1e3, UTM.lat = Y / 1e3) |>
  distinct()
grid_dog$survey_lumped <- "dog"

hbll_ins_s <- gfplot::hbll_inside_s_grid
hbll_ins_s$grid$area <- "hbll_s"
hbll_ins_n <- gfplot::hbll_inside_n_grid
hbll_ins_n$grid$area <- "hbll_n"
grid_hbll <- rbind(hbll_ins_n$grid, hbll_ins_s$grid)
grid_hbll <- add_utm_columns(
  grid_hbll,
  utm_crs = 32609,
  ll_names = c("X", "Y"),
  utm_names = c("UTM.lon", "UTM.lat")
)
grid_hbll$survey_lumped <- "hbll"

plot_data <- list()

for (model in models) {
message("Checking model: ", model)
fitbb <- readRDS(paste0("data-generated/fit-intonly", model, ".rds"))
d <- fitbb$data

# DHARMa residuals-----------
set.seed(1)
s_bb <- simulate(fitbb, nsim = 300, type = "mle-mvn")
max(s_bb) # bounded by the hooks (or 1 if simulated as proportions)

qq <- dharma_residuals(s_bb, object = fitbb, plot = FALSE)
res_bb <- dharma_residuals(s_bb, object = fitbb, return_DHARMa = TRUE)

# dispersion: variance of the observed residuals about the fitted values vs.
# the same statistic for each simulated dataset (as in DHARMa::testDispersion())
fitted_bb <- res_bb$fittedPredictedResponse
expected_var <- sd(res_bb$simulatedResponse)^2
spread <- function(x) var(x - fitted_bb) / expected_var
disp_sim <- apply(res_bb$simulatedResponse, 2, spread)
disp_obs <- spread(res_bb$observedResponse)
# DHARMa's reported statistic is the observed value over the mean simulated
# value (1 = as dispersed as the model implies); check it matches testDispersion()
stopifnot(isTRUE(all.equal(
  unname(DHARMa::testDispersion(res_bb, plot = FALSE)$statistic),
  disp_obs / mean(disp_sim)
)))
disp_obs <- disp_obs / mean(disp_sim)
disp_sim <- disp_sim / mean(disp_sim)

plot_data[[model]] <- list(qq = qq, disp_sim = disp_sim, disp_obs = disp_obs,
  res = res_bb)

# mean rate vs exp(eta)-----------
grid <- if (model == "dog-j-circle") grid_dog else grid_hbll
years <- seq(min(d$year), max(d$year), 1)
grid <- purrr::map_dfr(years, ~ tibble(grid, year = .x))
grid$julian_c <- 0

phi <- exp(fitbb$model$par[["ln_phi"]])
rates <- predict(fitbb, newdata = grid, offset = rep(0, nrow(grid))) |>
  mutate(
    p_bar = 1 - exp(-exp(est)),
    rate_mean = digamma(phi) - digamma((1 - p_bar) * phi),
    rate_eta = -log(1 - p_bar)
  )
print(quantile(rates$p_bar, c(0.5, 0.99)))
rates |>
  group_by(year) |>
  summarise(ratio = sum(rate_mean) / sum(rate_eta)) |>
  mutate(ratio = round(ratio / mean(ratio), 4)) |>
  as.data.frame() |>
  print()
}

# one figure: QQ (left) and dispersion (right) for each model---------
pdf("figures/dharma-checks.pdf", width = 7, height = 6)
par(mfrow = c(length(models), 3), mar = c(3.4, 3.4, 2, 0.8), mgp = c(2.1, 0.6, 0), cex = 0.7)
panel <- 0
for (model in models) {
  pd <- plot_data[[model]]
  panel <- panel + 1
  plot(pd$qq$expected, pd$qq$observed, pch = 16, cex = 0.5,
    xlab = "Expected", ylab = "Observed", xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, col = "red")
  title(main = paste0("(", letters[panel], ") ", model, ": QQ residuals"),
    adj = 0, cex.main = 0.85)
  panel <- panel + 1
  hist(pd$disp_sim, breaks = 20, col = "grey85", border = "grey60",
    xlim = range(pd$disp_sim, pd$disp_obs), main = "",
    xlab = "Dispersion (observed / mean simulated)")
  abline(v = 1, lty = 2, col = "grey50")
  abline(v = pd$disp_obs, col = "red", lwd = 2)
  title(main = paste0("(", letters[panel], ") ", model, ": dispersion"),
    adj = 0, cex.main = 0.85)
  panel <- panel + 1
  # residuals vs. fitted (rank) with the quantile-regression curves from
  # DHARMa::testQuantiles(); drawn by hand to avoid DHARMa's test text in the title
  x <- rank(pd$res$fittedPredictedResponse) / length(pd$res$fittedPredictedResponse)
  qr <- DHARMa::testQuantiles(pd$res, predictor = x, plot = FALSE)$predictions
  plot(x, pd$res$scaledResiduals, pch = 16, cex = 0.3,
    col = grDevices::adjustcolor("black", 0.4), xlim = c(0, 1), ylim = c(0, 1),
    xlab = "Fitted (rank)", ylab = "Scaled residual")
  abline(h = c(0.25, 0.5, 0.75), lty = 2, col = "grey60")
  for (k in c(1, 3, 5)) lines(qr$pred, qr[[as.character(k)]], col = "red", lwd = 1.5)
  title(main = paste0("(", letters[panel], ") ", model, ": residuals vs. fitted"),
    adj = 0, cex.main = 0.85)
}
dev.off()
