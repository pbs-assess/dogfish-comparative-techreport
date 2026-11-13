library(sdmTMB)
library(ggplot2)
library(tidyverse)
source("analysis/sdmTMB combined models/utils-mesh.R")

# Shallow vs. deep indices for the dogfish survey on the dogfish grid.
# Dogfish-survey version of 01b-index-depths.R (which does HBLL INS N/S).
#
# Depth cutoff: 165 m. This is a dogfish survey design-stratum boundary
# (strata 1-55, 56-110, 111-165, 166-220, >220 m), it's close to the median
# set depth (147 m) and grid-cell depth (~166 m), and it leaves at least 12
# sets per bin in every survey year. A 110 m cutoff leaves only 6-13 shallow
# sets in most years since 2005.

df <- readRDS("data-raw/wrangled-hbll-dog-sets.rds") |>
  drop_na(catch_count) |>
  drop_na(offset) |>
  drop_na(julian) |>
  filter(soak < 1000) # as in 01-index-generation.R

depth_cut <- 165

# dogfish survey only (same data as the "dog-j-circle" model) ------------------
d <- df |>
  filter(survey_sep %in% c("dog", "dog comp", "dog-jhook"))

# drop the 12/0 circle-hook comparison sets
id_remove <- d |>
  filter(grepl("COMPARISON", activity_desc) & hooksize_desc == "12/0" & year %in% c(2023, 2024)) |>
  pull(fishing_event_id)

d <- d |>
  filter(!fishing_event_id %in% id_remove) |>
  drop_na(depth_m) |>
  mutate(depth_bin = ifelse(depth_m <= depth_cut, 1, 2))

table(d$depth_bin, d$year)

d$catch_prop <- d$catch_count / d$lglsp_hook_count

# dogfish grid with depth ------------------------------------------------------
# gfplot::dogfish_grid has no depth, so take it from ETOPO 2022 15 arc-second
# bathymetry (nearest raster cell). Against the recorded set depths this gives
# r = 0.8 and a median difference of ~5 m. Cached after the first download.
bathy_file <- "data-raw/etopo-dogfish-bathy.rds"
if (!file.exists(bathy_file)) {
  b <- marmap::getNOAA.bathy(-125.4, -123.0, 48.7, 50.2, resolution = 0.25, keep = FALSE)
  bathy <- marmap::as.xyz(b) |> setNames(c("lon", "lat", "z"))
  saveRDS(bathy, bathy_file)
}
bathy <- readRDS(bathy_file)
bathy_xy <- sf::st_as_sf(bathy, coords = c("lon", "lat"), crs = 4326) |>
  sf::st_transform(32609) |>
  sf::st_coordinates() / 1e3

grid <- gfplot::dogfish_grid$grid |>
  mutate(UTM.lon = X / 1e3, UTM.lat = Y / 1e3) |>
  distinct()
nn <- RANN::nn2(bathy_xy, grid[, c("UTM.lon", "UTM.lat")], k = 1)
grid$depth_m <- -bathy$z[nn$nn.idx[, 1]]
summary(grid$depth_m)
grid <- filter(grid, depth_m > 0) # cells landing on land in the raster
grid <- mutate(grid, depth_bin = ifelse(depth_m <= depth_cut, 1, 2))
table(grid$depth_bin)

ggplot(grid, aes(UTM.lon, UTM.lat, fill = factor(depth_bin))) +
  geom_tile(width = 0.5, height = 0.5) +
  geom_point(data = d, aes(fill = factor(depth_bin)), pch = 21, size = 1.5) +
  coord_fixed() +
  labs(fill = "Depth bin")

years <- seq(min(d$year), max(d$year), 1)
grid <- replicate_df(grid, "year", years)
grid$survey_lumped <- "dog" # predict for the circle-hook survey

# Same model as 01-index-generation.R ("dog-j-circle"): beta-binomial with
# hooks as trials, cloglog link, offset = log(soak); rw0 year intercept, IID
# spatiotemporal fields and a spatial field, with PC priors on the fields.
# survey_lumped (dog vs. dog-jhook) calibrates the J-hook years.
family <- betabinomial(link = "cloglog")
sigma_v_prior <- lognormal_prior(log(0.3), 0.5)
formula <- catch_prop ~ 1 + as.factor(survey_lumped)

fit_bin <- function(dat, mesh) {
  sdmTMB(
    formula = formula,
    time = "year",
    spatiotemporal = "iid",
    time_varying = ~1,
    time_varying_type = "rw0",
    weights = dat$lglsp_hook_count,
    offset = log(dat$soak),
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
    mesh = mesh,
    data = dat,
    extra_time = seq(min(dat$year), max(dat$year))
  )
}

# shallow ----------------------------------------------------------------------
dsh <- filter(d, depth_bin == 1)
gridshallow <- filter(grid, depth_bin == 1)
meshsh <- make_survey_mesh(dsh, gridshallow)
plot(meshsh)

ms <- fit_bin(dsh, meshsh)
ms |> strip_sdmTMB() |>
  saveRDS("output/fit-depth-bins-dog-shallow.rds")
ms <- readRDS("output/fit-depth-bins-dog-shallow.rds")

# deep -------------------------------------------------------------------------
ddp <- filter(d, depth_bin == 2)
griddeep <- filter(grid, depth_bin == 2)
meshdp <- make_survey_mesh(ddp, griddeep)
plot(meshdp)

md <- fit_bin(ddp, meshdp)
md |> strip_sdmTMB() |>
  saveRDS("output/fit-depth-bins-dog-deep.rds")
md <- readRDS("output/fit-depth-bins-dog-deep.rds")

sanity(ms)
sanity(md)
ms$sd_report
md$sd_report

# indices ----------------------------------------------------------------------
index <- get_index(ms, newdata = gridshallow, offset = rep(0, nrow(gridshallow)),
  derived_link = "log", bias_correct = TRUE)
index$loc <- "shallow"
index$surveyed <- index$year %in% dsh$year
saveRDS(index, "output/ind-depth-bin-dog-shallow.rds")

index_d <- get_index(md, newdata = griddeep, offset = rep(0, nrow(griddeep)),
  derived_link = "log", bias_correct = TRUE)
index_d$loc <- "deep"
index_d$surveyed <- index_d$year %in% ddp$year
saveRDS(index_d, "output/ind-depth-bin-dog-deep.rds")

# plots: see 01d-plot-indices.R
