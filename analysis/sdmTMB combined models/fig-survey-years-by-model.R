# Which survey data go into each spatiotemporal model, by year.
# Mirrors the data filtering in 01-index-generation.R.
library(tidyverse)

df <- readRDS("data-raw/wrangled-hbll-dog-sets.rds") |>
  drop_na(catch_count, offset, julian) |>
  filter(soak < 1000)

# 12/0 hook comparison sets in 2023-2024 are dropped from the dogfish models
id_remove <- df |>
  filter(
    grepl("COMPARISON", activity_desc) & hooksize_desc == "12/0" &
      year %in% c(2023, 2024)
  ) |>
  pull(fishing_event_id)

d_hblldog <- filter(df, !fishing_event_id %in% id_remove)
d_hbll <- filter(df, survey_abbrev %in% c("HBLL INS S", "HBLL INS N"))
d_dog <- filter(
  df,
  survey_sep %in% c("dog", "dog comp", "dog-jhook"),
  !fishing_event_id %in% id_remove
)

model_labels <- c(
  hblldog = "(a) HBLL & SoG dogfish circle/J-hook",
  hbll = "(b) HBLL",
  dog = "(c) SoG dogfish circle/J-hook"
)

survey_labels <- c(
  "HBLL INS N" = "HBLL inside north",
  "HBLL INS S" = "HBLL inside south",
  "hbll comp" = "HBLL comparison sets",
  "dog" = "Dogfish (circle hook)",
  "dog comp" = "Dogfish comparison sets",
  "dog-jhook" = "Dogfish (J-hook)"
)

# survey_lumped is the factor that gets its own catchability in the models
# (as.factor(survey_lumped)); HBLL is the reference level
catch_labels <- c(
  "hbll" = "HBLL (hook 13/0)",
  "dog" = "Dogfish circle hook (14/0)",
  "dog-jhook" = "Dogfish J-hook (12/0)"
)

dat <- bind_rows(
  mutate(d_hblldog, model = "hblldog"),
  mutate(d_hbll, model = "hbll"),
  mutate(d_dog, model = "dog")
) |>
  count(model, survey_sep, survey_lumped, year, name = "n_sets") |>
  mutate(
    model = factor(model, levels = names(model_labels), labels = model_labels),
    survey_sep = factor(survey_sep, levels = names(survey_labels), labels = survey_labels),
    survey_lumped = factor(survey_lumped, levels = names(catch_labels), labels = catch_labels)
  )

# Okabe-Ito, colour-blind safe
cols <- setNames(c("#0072B2", "#D55E00", "#CC79A7"), catch_labels)

g <- ggplot(dat, aes(year, survey_sep, colour = survey_lumped)) +
  # translucent fill underneath, opaque outline on top
  geom_point(aes(size = n_sets), shape = 16, alpha = 0.2) +
  geom_point(aes(size = n_sets), shape = 1, stroke = 1) +
  scale_size_area(max_size = 5, breaks = c(10, 25, 50, 75)) +
  ggforce::facet_col(~model, scales = "free_y", space = "free") +
  scale_colour_manual(values = cols, guide = guide_legend(override.aes = list(size = 4)))  +
  scale_x_continuous(breaks = seq(1985, 2025, 5), limits = range(df$year), minor_breaks = seq(1985, 2026, 1)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Year", y = NULL, colour = "Catchability group", size = "Number of sets") +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title.position = "top",
    strip.text = element_text(hjust = 0, colour = "grey10"),
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(colour = "grey88"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92")
  )

g
ggsave("figures/survey-years-by-model.pdf", g, width = 7, height = 5.5)
