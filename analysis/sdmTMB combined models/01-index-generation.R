library(ggplot2)
library(tidyverse)
library(sdmTMB)
source("analysis/sdmTMB combined models/utils-mesh.R")


# models to run; each is fit, saved, and indexed in turn
models <- c("hblldog", "hbll-n-s", "dog-j-circle")


grid_dog <- gfplot::dogfish_grid$grid %>%
  mutate(UTM.lon = X / 1e3, UTM.lat = Y / 1e3) |>
  distinct()

hbll_ins_s <- gfplot::hbll_inside_s_grid
hbll_ins_s$grid$area <- "hbll_s"
hbll_ins_n <- gfplot::hbll_inside_n_grid
hbll_ins_n$grid$area <- "hbll_n"
grid_hbll <- rbind(hbll_ins_n$grid, hbll_ins_s$grid)
crs = 32609
grid_hbll <- add_utm_columns(
  grid_hbll,
  utm_crs = crs,
  ll_names = c("X", "Y"),
  utm_names = c("UTM.lon", "UTM.lat")
)

df <- readRDS("data-raw/wrangled-hbll-dog-sets.rds") |>
  drop_na(catch_count) |>
  drop_na(offset) |>
  drop_na(julian) |>
  filter(soak < 1000) #not sure what happened but that seems very high so removing

df$log_botdepth2 <- df$log_botdepth * df$log_botdepth
# Likelihood: beta-binomial with hooks as trials, cloglog link and
# offset = log(soak), so exp(eta) is a catch rate per hook-minute and the
# cloglog builds in hook saturation (at most one fish per hook). The index sums
# that rate over the grid (derived_link = "log") at offset 0.
#
# Time structure: year intercept as a random walk with the first year fixed at
# the intercept ("rw0"), IID spatiotemporal fields, and a spatial field.
# survey_lumped enters as a factor alongside the intercept, which "rw0" allows.
df$catch_prop <- df$catch_count / df$lglsp_hook_count

# keeps the RW SD off zero: median 0.3, 95% interval about 0.11 to 0.8
sigma_v_prior <- lognormal_prior(log(0.3), 0.5)


for (model in models) {
message("Running model: ", model)
if (model == "hblldog") {
  d <- df
  sort(unique(d$year))
  years <- seq(min(d$year), max(d$year), 1)

  id_remove <- d %>%
    filter(grepl("COMPARISON", activity_desc) & hooksize_desc == "12/0" & year %in% c(2023, 2024)) |>
    pull(fishing_event_id)

  d <- d |>
    filter(!fishing_event_id %in% id_remove) |>
    drop_na(offset) |>
    drop_na(julian) |>
    drop_na(catch_count)

  meanj <- mean(d$julian)
  d$julian_c <- d$julian - meanj

  ggplot(d, aes(year, catch_prop)) +
    geom_point() +
    facet_wrap(~survey_sep)

  grid <- grid_hbll
  grid <- purrr::map_dfr(years, ~ tibble(grid, year = .x))
  grid$julian_c <- 0
  grid$survey_lumped <- "hbll"

  formulabb <- catch_prop ~ 1 + as.factor(survey_lumped) + s(julian_c)
}

if (model == "hbll-n-s") {
  d <- df |>
    filter(survey_abbrev %in% c("HBLL INS S", "HBLL INS N")) |>
    filter(survey_abbrev != "OTHER") |> # no comparative work, therefore no need for julian but could try
    drop_na(offset) |>
    drop_na(julian) |>
    drop_na(catch_count)

  meanj <- mean(d$julian)
  d$julian_c <- d$julian - meanj

  years <- seq(min(d$year), max(d$year), 1)
  grid <- grid_hbll
  grid <- purrr::map_dfr(years, ~ tibble(grid, year = .x))
  grid$julian_c <- 0

  formulabb <- catch_prop ~ 1
}

if (model == "dog-j-circle") {
  d <- df |>
    filter(survey_sep %in% c("dog", "dog comp", "dog-jhook"))

  id_remove <- d %>%
    filter(grepl("COMPARISON", activity_desc) & hooksize_desc == "12/0" & year %in% c(2023, 2024)) |>
    pull(fishing_event_id)

  d <- d |>
    filter(!fishing_event_id %in% id_remove) |>
    drop_na(offset) |>
    drop_na(julian) |>
    drop_na(catch_count)

  meanj <- mean(d$julian)
  d$julian_c <- d$julian - meanj

  grid_dog <- gfplot::dogfish_grid$grid %>%
    mutate(UTM.lon = X / 1e3, UTM.lat = Y / 1e3) |>
    distinct()
  years <- seq(min(d$year), max(d$year), 1)
  grid_dog <- replicate_df(grid_dog, "year", unique(years))
  grid_dog$julian_c <- 0
  grid_dog$survey_lumped <- "dog"
  grid <- grid_dog

  formulabb <- catch_prop ~ 1 + as.factor(survey_lumped)
}


# index generation--------


mesh <- make_survey_mesh(d, grid)
plot(mesh)
mesh$mesh$n

g <-
  ggplot() +
  geom_point(
    data = d |> arrange(year), aes(UTM.lon, UTM.lat, colour = year),
    size = 1,
    # ), size = 0.25, alpha = 0.25,
    pch = 16
  ) +
  inlabru::gg(mesh$mesh,
    edge.color = "grey60",
    edge.linewidth = 0.15,
    # interior = TRUE,
    # int.color = "blue",
    int.linewidth = 0.25,
    exterior = FALSE,
    # ext.color = "black",
    ext.linewidth = 0.5
  ) +
  scale_colour_viridis_c(direction = -1) +
  labs(colour = "Year") +
  xlab("UTM (km)") +
  ylab("UTM (km)") +
  coord_fixed(expand = FALSE) +
  theme_classic() +
  theme(legend.position = "right")
# theme(legend.position = "inside", legend.position.inside = c(0.2, 0.25))
g
saveRDS(g, paste0("figures/sog-model-mesh-", model, ".rds"))


ggplot() +
  geom_point(
    data = d |> arrange(year), aes(UTM.lon, UTM.lat, colour = year),
    size = 1,
    # ), size = 0.25, alpha = 0.25,
    pch = 16
  ) +
  scale_colour_viridis_c(direction = -1) +
  facet_grid(~survey_lumped)

ggplot() +
  geom_point(
    data = d |> arrange(year), aes(year, catch_prop),
    size = 1,
    # ), size = 0.25, alpha = 0.25,
    pch = 16
  ) +
  scale_colour_viridis_c(direction = -1) +
  facet_grid(~survey_lumped)

ggplot() +
  geom_point(
    data = d |> arrange(year), aes(offset, catch_prop),
    size = 1,
    # ), size = 0.25, alpha = 0.25,
    pch = 16
  ) +
  scale_colour_viridis_c(direction = -1) +
  facet_grid(~survey_sep)

ggplot() +
  geom_point(
    data = d |> arrange(year), aes(month, catch_prop),
    size = 1,
    # ), size = 0.25, alpha = 0.25,
    pch = 16
  ) +
  scale_colour_viridis_c(direction = -1) +
  facet_grid(~survey_sep)


# fit <- sdmTMB(
#     formula = formulanb,
#     data = d,
#     time = "year",
#     mesh = mesh,
#     spatial = "on", #off for hbll
#     spatiotemporal = "rw",
#     family = nbinom2(),
#     offset = d$offset_hksoak,
#     silent = TRUE,
#     share_range = FALSE,
#     extra_time = extratime,
#     control = sdmTMBcontrol(newton_loops = 1L)
#   )
#
# fit$sd_report
# sanity(fit)

fitbb <- sdmTMB(
  formula = formulabb,
  data = d,
  time = "year",
  mesh = mesh,
  spatial = "on",
  spatiotemporal = "iid",
  time_varying = ~1,
  time_varying_type = "rw0",
  weights = d$lglsp_hook_count, # trials = hooks
  offset = log(d$soak), # exposure per hook
  family = betabinomial(link = "cloglog"),
  silent = FALSE,
  share_range = TRUE, # weren't that different before; simpler model
  priors = sdmTMBpriors(
    # PC priors on the Matern fields: P(range < 10 km) = 0.05 keeps the range
    # off scales the mesh can't represent; P(SD > 2) = 0.05 is weakly informative
    matern_s = pc_matern(range_gt = 10, sigma_lt = 2),
    matern_st = pc_matern(range_gt = 10, sigma_lt = 2),
    sigma_V = sigma_v_prior
  ),
  extra_time = seq(min(d$year), max(d$year)), # sdmTMB skips years in the data
  control = sdmTMBcontrol(multiphase = FALSE, collapse_spatial_variance = TRUE)
)

AIC(fitbb)
sanity(fitbb)
# AIC(fit)

ok <- all(unlist(sanity(fitbb)))

if (!ok) {
  warning("sanity() checks failed for model ", model, "; skipping save and index")
  next
}
fitbb |> strip_sdmTMB() |>
  saveRDS(file = paste0("data-generated/fit-intonly", model, ".rds"))

# get index
index <- get_index(fitbb, newdata = grid, offset = rep(0, nrow(grid)),
  derived_link = "log", bias_correct = TRUE)
index <- index |> mutate(model = ifelse(year %in% unique(d$year), "yrs_surved", "yrs_interp"))
index$modelloc <- model
index$type <- "int"
saveRDS(index, file = paste0("data-generated/ind-sog-intonly", model, ".rds"))
}
