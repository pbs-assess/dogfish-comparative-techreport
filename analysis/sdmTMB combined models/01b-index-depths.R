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

meshsh <- make_survey_mesh(dsh, gridshallow, max_edge = 6)
plot(meshsh)

ddp <- filter(d, depth_bin == 2)
griddeep <- filter(grid, depth_bin == 2)

meshdp <- make_survey_mesh(ddp, griddeep, max_edge = 6)
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

saveRDS(ms, "output/fit-depth-bins-shallow.rds")
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
saveRDS(md, "output/fit-depth-bins-deep.rds")
md <- readRDS("output/fit-depth-bins-deep.rds")

sanity(ms)
sanity(md)
ms$sd_report
md$sd_report


index <- get_index(ms, newdata = gridshallow, offset = rep(0, nrow(gridshallow)),
  derived_link = "log", bias_correct = TRUE)
index$loc <- "shallow"
index <- index |>
  mutate(mean = mean(est), meanlwr = mean(lwr), meanupr = mean(upr))

saveRDS(index, "output/ind-depth-bin-shallow.rds")

index_d <- get_index(md, newdata = griddeep, offset = rep(0, nrow(griddeep)),
  derived_link = "log", bias_correct = TRUE)
index_d$loc <- "deep"
index_d <- index_d |>
  mutate(mean = mean(est), meanlwr = mean(lwr), meanupr = mean(upr))

saveRDS(index_d, "output/ind-depth-bin-deep.rds")

index2 <- bind_rows(index, index_d)

index2 |>
  group_by(loc) |>
  #mutate(est2 = est - est[depth_m == min(depth_m)]) |>
  mutate(est2 = est - mean) |>
  mutate(est2 = est - mean, lwrc = lwr - mean, uprc = upr - mean) |>
  ggplot() +
  geom_line(aes(year, (est2), group = loc, colour = loc)) +
  geom_pointrange(mapping = aes(x = year, y = est2, ymin = (lwrc), ymax = (uprc),  group = loc, colour = loc), size = 0.2, pch = 5, alpha = 0.6, position = position_dodge(width = 0.5)) +
  scale_colour_viridis_d() +
  theme_classic()

index2 |>
  group_by(loc) |>
  #mutate(est2 = est - est[depth_m == min(depth_m)]) |>
  mutate(est2 = est - mean, lwrc = lwr - mean, uprc = upr - mean) |>
  ggplot(aes(year, (est2), ymin = (lwrc), ymax = (uprc), group = loc, colour = loc)) +
  geom_line(aes(year, (est2), group = loc, colour = loc)) +
  geom_pointrange(mapping = aes(x = year - 0.25), size = 1, pch = 5, alpha = 0.6, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("#d8b365", "#5ab4ac")) +
  guides(colour = guide_legend(title = "Depth bin")) +
  theme_classic() +
  ylab("Abundance index") +
  #guides(colour = "Depth bin") +
  xlab("Year")
ggsave("figures/depth-bin-figure.jpg", width = 4, height = 3)


