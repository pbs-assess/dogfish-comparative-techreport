# Plot the indices saved by 01-index-generation.R and 01b-index-depths.R. Kept separate from the
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
g_depth <- ind_depth |>
  mutate(est2 = est / gmean, lwrc = lwr / gmean, uprc = upr / gmean) |>
  ggplot(aes(year, est2, ymin = lwrc, ymax = uprc, colour = loc, fill = loc)) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.3) +
  geom_linerange(linewidth = 0.4, alpha = 0.7, position = dodge) +
  geom_line(aes(group = loc), linewidth = 0.4, alpha = 0.7, position = dodge) +
  geom_point(size = 1.4, shape = 21, colour = "white", stroke = 0.3, position = dodge) +
  scale_colour_manual(values = c("#d8b365", "#5ab4ac")) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac")) +
  labs(x = "Year", y = "Abundance index (relative to geometric mean)", colour = "Depth bin", fill = "Depth bin") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_sleek(base_size = 9) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85))
g_depth
ggsave("figures/depth-bin-figure.pdf", g_depth, width = 4, height = 3)
