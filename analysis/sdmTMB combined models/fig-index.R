# combined indices and centered and scale indices figures
library(tidyverse)
library(ggplot2)

index1 <- readRDS(file = paste0("data-generated/ind-sog-intonly", "hblldog", ".rds")) |>
  mutate(gmean = exp(mean(log_est))) |> # geometric mean of the index
  mutate(est_c = est / gmean, lwr_c = lwr / gmean, upr_c = upr / gmean)

index2 <- readRDS(file = paste0("data-generated/ind-sog-intonly", "hbll-n-s", ".rds")) |>
  mutate(gmean = exp(mean(log_est))) |> # geometric mean of the index
  mutate(est_c = est / gmean, lwr_c = lwr / gmean, upr_c = upr / gmean)

index3 <- readRDS(file = paste0("data-generated/ind-sog-intonly", "dog-j-circle", ".rds")) |>
  mutate(gmean = exp(mean(log_est))) |> # geometric mean of the index
  mutate(est_c = est / gmean, lwr_c = lwr / gmean, upr_c = upr / gmean)

ind <- bind_rows(index1, index2, index3)
ind$group <- paste0(ind$modelloc, ind$type)
ind$modelloc <- factor(ind$modelloc, levels = c("dog-j-circle", "hbll-n-s", "hblldog"))

x <- palette.colors(palette = "Okabe-Ito")

gg <- ind %>%
  # mutate(dyear = year %in% year_dogfish) %>%
  # ggplot(aes(year, est, ymin = lwr, ymax = upr, group = Survey, colour = Survey)) +
  ggplot(aes(year, est_c, ymin = lwr_c, ymax = upr_c, group = modelloc, colour = modelloc)) +
  geom_line(aes(group = modelloc, colour = modelloc), linewidth = 1) +
  geom_point(aes(group = modelloc, colour = modelloc), size = 2) +
  geom_ribbon(aes(year, est_c, ymin = lwr_c, ymax = upr_c, fill = modelloc), alpha = 0.5) +
  # geom_linerange() +
  # facet_wrap(vars(Survey), ncol = 1, scales = "free_y") +
  expand_limits(y = 0) +
  scale_colour_manual(values = x[c(4, 7, 2)]) +
  scale_fill_manual(values = x[c(4, 7, 2)], labels = c("SoG dog circle/J-hook", "HBLL", "HBLL & SoG circle/J-hook")) +
  labs(x = "Year", y = "Relative abundance index") +
  theme_classic() +
  labs(fill = "Surveys incl.") +
  guides(colour = "none")

gg

ggsave("figures/stitched_indexes.jpg", gg, width = 5, height = 4)


# drop interpolated (grey) years before the HBLL starts; they carry no data
first_hbll <- min(filter(ind, modelloc == "hbll-n-s", model == "yrs_surved")$year)
ind_panel <- ind |> filter(!(model == "yrs_interp" & year < first_hbll))

gg2 <- ind_panel %>%
  # filter(model == "yrs_surved") |>
  ggplot(aes(year, est_c, ymin = lwr_c, ymax = upr_c, group = modelloc, colour = modelloc)) +
  # geom_line(data = filter(ind_panel, model == "yrs_surved"), aes(group = modelloc, colour = modelloc), linewidth = 1) +
  geom_point(data = filter(ind_panel, model == "yrs_surved"), aes(group = modelloc, colour = modelloc), size = 2) +
  geom_point(data = filter(ind_panel, model == "yrs_interp"), aes(group = modelloc), colour = "grey80", size = 2) +
  geom_pointrange(data = filter(ind_panel, model == "yrs_interp"), position = position_dodge(width = 0.25), colour = "grey80") +
  geom_pointrange(data = filter(ind_panel, model == "yrs_surved"), position = position_dodge(width = 0.25)) +
  expand_limits(y = 0) +
  scale_fill_manual(values = x[c(4, 7, 2)]) +
  scale_colour_manual(values = x[c(4, 7, 2)]) +
  labs(x = "Year", y = "Relative abundance index") +
  facet_wrap(~modelloc, nrow = 3, ncol = 1,
    labeller = labeller(
      modelloc = c(
        "hbll-n-s" = "HBLL",
        "dog-j-circle" = "SoG dog circle/J-hook",
        "hblldog" = "HBLL & SoG circle/J-hook"
      )
    )
  ) +
  guides(colour = "none") +
  ggsidekick::theme_sleek()

gg2

ggsave("figures/stitched_indexes_panel.pdf", gg2, width = 4, height = 9)
