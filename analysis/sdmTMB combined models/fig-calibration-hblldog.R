# Calibration factors (survey group effects) from the combined hbll-dogfish
# model (fit-intonlyhblldog.rds). survey_lumped enters the cloglog linear
# predictor, so exp(coef) is a multiplicative effect on catch rate.
# Two panels: (a) relative to the dogfish circle hook, (b) relative to the
# dogfish J-hook (the original gear, 1986-2004).
#
# Also shown for HBLL: the calibration from the paired comparative sets
# (analysis/calibration/01-...) and seasonally paired comparative sets
# (analysis/calibration/02-...), each an intercept-only beta-binomial model with
# rho = HBLL : dogfish circle-hook catch ratio. In (a) rho is shown directly. The
# comparative sets contain no J-hooks, so in (b) it is chained with the combined
# model's dogfish circle vs. J-hook effect (independent estimates; log-scale SEs
# add in quadrature).
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

lev <- c("Dogfish circle hook", "Dogfish J-hook", "HBLL (circle hook)")
src_model <- "Combined model"
src_paired <- "Paired comparative sets"
src_season <- "Seasonally paired comparative sets"

# contrasts beta_i - beta_base (base = index of the reference level in beta)
contrast_tbl <- function(base, panel) {
  L <- diag(3) - matrix(as.numeric(seq_len(3) == base), 3, 3, byrow = TRUE)
  tibble(
    panel = panel, survey = lev, source = src_model,
    est = as.vector(L %*% beta), se = sqrt(diag(L %*% V %*% t(L)))
  )
}

# paired-set rhos are HBLL : dog circle (log scale). On the circle-hook base they
# are used directly; on the J-hook base they are chained with the combined
# model's dog circle : J-hook effect (independent; SEs add in quadrature).
rho_files <- c(
  setNames("data-generated/calibration_coeffs_paired_intercept.rds", src_paired),
  setNames("data-generated/calibration_coeffs_intercept.rds", src_season)
)
paired_tbl <- function(panel, chain_to = NULL) {
  imap_dfr(rho_files, function(file, source) {
    rho <- readRDS(file)
    e <- rho$estimate
    se <- rho$std.error
    if (!is.null(chain_to)) { # chain_to: row of the J-hook-base contrast for dog circle
      e <- e + chain_to$est
      se <- sqrt(se^2 + chain_to$se^2)
    }
    tibble(panel = panel, survey = lev[3], source = source, est = e, se = se)
  })
}

pn_circle <- "(a) Relative to dogfish circle hook"
pn_jhook <- "(b) Relative to dogfish J-hook"
cal_j <- contrast_tbl(2, pn_jhook)
cal <- bind_rows(
  contrast_tbl(1, pn_circle), paired_tbl(pn_circle),
  cal_j, paired_tbl(pn_jhook, chain_to = cal_j[1, ])
) |>
  mutate(
    estimate = exp(est), conf.low = exp(est - 1.96 * se), conf.high = exp(est + 1.96 * se),
    survey = factor(survey, levels = rev(lev)),
    source = factor(source, levels = c(src_model, src_paired, src_season)),
    panel = factor(panel, levels = c(pn_circle, pn_jhook))
  )
print(cal, n = Inf)

x <- palette.colors(palette = "Okabe-Ito")

g <- ggplot(cal, aes(estimate, survey,
  xmin = conf.low, xmax = conf.high, colour = survey, shape = source)) +
  geom_vline(xintercept = 1, lty = 2, colour = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.6, linewidth = 0.8) +
  facet_wrap(~panel, scales = "free_x") +
  scale_colour_manual(values = unname(x[c(4, 7, 2)]), guide = "none") +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_x_log10(breaks = c(0.8, 1, 1.25, 1.5, 2, 3)) +
  labs(x = "Catchability ratio", y = NULL, shape = NULL) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(nrow = 1))
g

ggsave("figures/calibration-factors-hblldog.pdf", g, width = 8, height = 3.2)
