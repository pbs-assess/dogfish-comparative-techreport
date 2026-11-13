library(sdmTMB)
library(ggplot2)
library(tidyverse)
source("analysis/sdmTMB combined models/utils-mesh.R")

df <- readRDS("data-raw/wrangled-hbll-dog-sets.rds") |>
  drop_na(catch_count) |>
  drop_na(offset) |>
  drop_na(julian)

df$log_botdepth2 <- df$log_botdepth * df$log_botdepth

depth <- df |>
  dplyr::select(depth_m, grouping_depth_id) |>
  distinct()

df <- df |>
  mutate(depth_bin = case_when(
    depth_m <= 70 ~ 1,
    depth_m > 70 & depth_m <= 110 ~ 2,
    depth_m > 110 & depth_m <= 165 ~ 3,
    depth_m > 165 & depth_m <= 220  ~ 4,
    depth_m > 220 ~ 5 ))

# hbll n and s only -------------------------------------------------------------
d <- df |>
  filter(survey_abbrev %in% c("HBLL INS S", "HBLL INS N")) |>
  #filter(survey_lumped == "hbll") |>
  drop_na(offset) |>
  drop_na(depth_m) |>
  drop_na(julian) |>
  drop_na(catch_count)
range(d$depth_m)
d <- d |> #4 points are deeper than 110, instead of dropping ill change the categorizations
  mutate(depth_bin = case_when(
    depth_m <= 70 ~ 1,
    depth_m > 70 & depth_m <= 150 ~ 2))

# Same model as 01-index-generation.R: beta-binomial with hooks as trials,
# cloglog link, offset = log(soak); rw0 year intercept, IID spatiotemporal
# fields and a spatial field, with PC priors on the fields.
family <- betabinomial(link = "cloglog")
sigma_v_prior <- lognormal_prior(log(0.3), 0.5)






# shallow deep hbll block model ------------------------------------------------------
d$year_factor <- as.factor(d$year)

#mean <- mean(d$log_botdepth)
#d$log_botdepth_c <- d$log_botdepth - mean
#d$log_botdepth_c2 <- d$log_botdepth_c * d$log_botdepth_c

d$catch_prop <- d$catch_count/d$lglsp_hook_count

dsh <- filter(d, depth_bin == 1)

source("analysis/999-load-hbll-n-s-grid.R") #this is the hbll n, s grid from gfdata
years <- seq(min(d$year), max(d$year), 1)
grid <- purrr::map_dfr(years, ~ tibble(grid, year = .x))
gridshallow <- filter(grid, depth_bin == 1)
griddeep <- filter(grid, depth_bin == 2)

meshsh <- make_survey_mesh(dsh, gridshallow)
plot(meshsh)

ddp <- filter(d, depth_bin == 2)
griddeep <- filter(grid, depth_bin == 2)

meshdp <- make_survey_mesh(ddp, griddeep)
plot(meshdp)

ms <- sdmTMB(
  formula = catch_prop ~ 1,
  time = "year",
  spatiotemporal = "iid",
  time_varying = ~1,
  time_varying_type = "rw0",
  weights = dsh$lglsp_hook_count,
  offset = log(dsh$soak),
  silent = FALSE,
  spatial = "on",
  share_range = TRUE,
  family = family,
  priors = sdmTMBpriors(
    # PC priors on the Matern fields: P(range < 10 km) = 0.05 keeps the range
    # off scales the mesh can't represent; P(SD > 2) = 0.05 is weakly informative
    matern_s = pc_matern(range_gt = 10, sigma_lt = 2),
    matern_st = pc_matern(range_gt = 10, sigma_lt = 2),
    sigma_V = sigma_v_prior
  ),
  control = sdmTMBcontrol(collapse_spatial_variance = TRUE, multiphase = FALSE),
  mesh = meshsh,
  data = dsh,
  do_index = FALSE,
  extra_time = seq(min(dsh$year), max(dsh$year))
)

ms |> strip_sdmTMB() |>
  saveRDS("output/fit-depth-bins-shallow.rds")
ms <- readRDS("output/fit-depth-bins-shallow.rds")

md <- sdmTMB(
  formula = catch_prop ~ 1,
  time = "year",
  spatiotemporal = "iid",
  time_varying = ~1,
  time_varying_type = "rw0",
  weights = ddp$lglsp_hook_count,
  offset = log(ddp$soak),
  silent = FALSE,
  spatial = "on",
  share_range = TRUE,
  family = family,
  priors = sdmTMBpriors(
    # PC priors on the Matern fields: P(range < 10 km) = 0.05 keeps the range
    # off scales the mesh can't represent; P(SD > 2) = 0.05 is weakly informative
    matern_s = pc_matern(range_gt = 10, sigma_lt = 2),
    matern_st = pc_matern(range_gt = 10, sigma_lt = 2),
    sigma_V = sigma_v_prior
  ),
  control = sdmTMBcontrol(collapse_spatial_variance = TRUE, multiphase = FALSE),
  mesh = meshdp,
  data = ddp,
  do_index = FALSE,
  extra_time = seq(min(ddp$year), max(ddp$year))
)
md |> strip_sdmTMB() |>
  saveRDS("output/fit-depth-bins-deep.rds")
md <- readRDS("output/fit-depth-bins-deep.rds")

sanity(ms)
sanity(md)
ms$sd_report
md$sd_report


index <- get_index(ms, newdata = gridshallow, offset = rep(0, nrow(gridshallow)),
  derived_link = "log", bias_correct = TRUE)
index$loc <- "shallow"

saveRDS(index, "output/ind-depth-bin-shallow.rds")

index_d <- get_index(md, newdata = griddeep, offset = rep(0, nrow(griddeep)),
  derived_link = "log", bias_correct = TRUE)
index_d$loc <- "deep"

saveRDS(index_d, "output/ind-depth-bin-deep.rds")

# plots: see 01d-plot-indices.R
