# Plot the indices saved by 01-index-generation.R, 01b-index-depths.R and
# 01e-index-depths-dog.R. Kept separate from the
# fitting loop so the plots are quick to tweak.
library(ggplot2)
library(tidyverse)
library(ggsidekick)

models <- c("hblldog", "hbll-n-s", "dog-j-circle")

ind <- map_dfr(models, \(m) readRDS(paste0("data-generated/ind-sog-intonly", m, ".rds")))
ind$modelloc <- factor(ind$modelloc, levels = models)

# all years, by model
g1 <- ggplot(ind, aes(as.factor(year), est, ymin = lwr, ymax = upr)) +
  geom_pointrange(position = position_dodge(width = 0.25)) +
  facet_wrap(~modelloc, ncol = 1, scales = "free_y") +
  theme_classic()
g1
ggsave("figures/index-all-years.pdf", g1, width = 5, height = 7)

# surveyed years only
g2 <- ind |>
  filter(model == "yrs_surved") |>
  ggplot(aes(as.factor(year), est, ymin = lwr, ymax = upr)) +
  geom_pointrange(position = position_dodge(width = 0.25)) +
  facet_wrap(~modelloc, ncol = 1, scales = "free_y") +
  theme_classic()
g2
ggsave("figures/index-surveyed-years.pdf", g2, width = 5, height = 7)

# depth-bin indices (01b-index-depths.R)--------------------------------
# each bin's index is divided by its geometric mean
ind_depth <- bind_rows(
  readRDS("output/ind-depth-bin-shallow.rds"),
  readRDS("output/ind-depth-bin-deep.rds")
) |>
  group_by(loc) |>
  mutate(gmean = exp(mean(log_est))) |>
  ungroup()

dodge <- position_dodge(width = 0.5)
depth_cols <- c(shallow = "#d8b365", deep = "#5ab4ac") # deep = darker
g_depth <- ind_depth |>
  mutate(est2 = est / gmean, lwrc = lwr / gmean, uprc = upr / gmean) |>
  ggplot(aes(year, est2, ymin = lwrc, ymax = uprc, colour = loc, fill = loc)) +
  geom_linerange(position = dodge) +
  geom_line(aes(group = loc), position = dodge) +
  geom_point(shape = 21, colour = "white", position = dodge, size = 2) +
  scale_colour_manual(values = depth_cols, breaks = c("shallow", "deep")) +
  scale_fill_manual(values = depth_cols, breaks = c("shallow", "deep")) +
  labs(x = "Year", y = "Abundance index", colour = "Depth bin", fill = "Depth bin") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_sleek() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.3)) +
  scale_x_continuous(breaks = seq(2000, 2030, 5)) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85))
g_depth
ggsave("figures/depth-bin-figure.pdf", g_depth, width = 5, height = 4)

# dogfish survey depth-bin indices (01e-index-depths-dog.R)----------------
# same style as above; surveyed years only, each bin divided by its
# geometric mean over those years
ind_depth_dog <- bind_rows(
  readRDS("output/ind-depth-bin-dog-shallow.rds"),
  readRDS("output/ind-depth-bin-dog-deep.rds")
) |>
  filter(surveyed) |>
  group_by(loc) |>
  mutate(gmean = exp(mean(log_est))) |>
  ungroup()

g_depth_dog <- ind_depth_dog |>
  mutate(est2 = est / gmean, lwrc = lwr / gmean, uprc = upr / gmean) |>
  ggplot(aes(year, est2, ymin = lwrc, ymax = uprc, colour = loc, fill = loc)) +
  geom_linerange(position = dodge) +
  geom_line(aes(group = loc), position = dodge, alpha = 0.4) +
  geom_point(shape = 21, colour = "white", position = dodge, size = 2) +
  scale_colour_manual(values = depth_cols, breaks = c("shallow", "deep")) +
  scale_fill_manual(values = depth_cols, breaks = c("shallow", "deep")) +
  labs(x = "Year", y = "Abundance index", colour = "Depth bin", fill = "Depth bin") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_sleek() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.3)) +
  scale_x_continuous(breaks = seq(1985, 2025, 5)) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85))
g_depth_dog
ggsave("figures/depth-bin-figure-dog.pdf", g_depth_dog, width = 5, height = 4)
ggsave("figures/depth-bin-figure-dog.png", g_depth_dog, width = 5, height = 4)
