# Calibration factors (survey group effects) from the combined hbll-dogfish
# model (fit-intonlyhblldog.rds), relative to the dogfish J-hook (the original
# gear, 1986-2004). survey_lumped enters the cloglog linear predictor, so
# exp(coef) is a multiplicative effect on catch rate.
#
# Also shown: the HBLL calibration from the paired comparative sets
# (analysis/calibration/01-...) and seasonally paired comparative sets
# (analysis/calibration/02-...), each an intercept-only beta-binomial model with
# rho = HBLL : dogfish circle-hook catch ratio. The paired sets contain no J-hooks, so to put it on the J-hook base it is chained
# with the combined model's dogfish-circle vs. J-hook effect (independent
# estimates; log-scale SEs add in quadrature).
#
# Contrasts are done by hand because emmeans::emmeans() fails on this fit
# (the s(julian_c) smooth breaks emmeans' recover_data()).
library(sdmTMB)
library(tidyverse)

fit <- readRDS("data-generated/fit-intonlyhblldog.rds")

# fixed effects: (Intercept), dog-jhook, hbll (dog circle hook is the reference)
b <- tidy(fit, "fixed")
beta <- c(dog = 0, b$estimate[b$term == "as.factor(survey_lumped)dog-jhook"],
  b$estimate[b$term == "as.factor(survey_lumped)hbll"])
V <- matrix(0, 3, 3)
V[2:3, 2:3] <- vcov(fit)[2:3, 2:3]

L <- diag(3) - matrix(c(0, 1, 0), 3, 3, byrow = TRUE) # beta_i - beta_jhook
est <- as.vector(L %*% beta)
se <- sqrt(diag(L %*% V %*% t(L)))

lev <- c("Dogfish circle hook", "Dogfish J-hook", "HBLL (circle hook)")
src_model <- "Combined model"
src_paired <- "Paired comparative sets"
src_season <- "Seasonally paired comparative sets"

cal <- tibble(
  survey = lev, source = src_model,
  est = est, se = se
)

# paired-set rhos (log scale), chained onto the J-hook base
chain <- function(file, source) {
  rho <- readRDS(file)
  tibble(
    survey = lev[3], source = source,
    est = rho$estimate + est[1], # log(HBLL/dog circle) + log(dog circle/J-hook)
    se = sqrt(rho$std.error^2 + se[1]^2)
  )
}
cal <- bind_rows(
  cal,
  chain("data-generated/calibration_coeffs_paired_intercept.rds", src_paired),
  chain("data-generated/calibration_coeffs_intercept.rds", src_season)
)

cal <- cal |>
  mutate(
    estimate = exp(est), conf.low = exp(est - 1.96 * se), conf.high = exp(est + 1.96 * se),
    survey = factor(survey, levels = rev(lev)),
    source = factor(source, levels = c(src_model, src_paired, src_season))
  )
print(cal)

x <- palette.colors(palette = "Okabe-Ito")

g <- ggplot(cal, aes(estimate, survey,
  xmin = conf.low, xmax = conf.high, colour = survey, shape = source)) +
  geom_vline(xintercept = 1, lty = 2, colour = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.6, linewidth = 0.8) +
  scale_colour_manual(values = unname(x[c(4, 7, 2)]), guide = "none") +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3)) +
  labs(x = "Catchability relative to dogfish J-hook", y = NULL, shape = NULL) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(nrow = 3))
g

ggsave("figures/calibration-factors-hblldog-jhook-base.pdf", g, width = 5.5, height = 3.6)
