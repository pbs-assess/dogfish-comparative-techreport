# Compare candidate meshes for the index models: the current
# make_mesh(cutoff = 5) against fmesher meshes built by make_survey_mesh()
# at several inner max.edge values. One figure per data set.
#
# Run from the repo root.

library(tidyverse)
library(sdmTMB)
source("analysis/sdmTMB combined models/utils-mesh.R")

dir.create("figures", showWarnings = FALSE)

max_edges <- c(8, 6, 5) # km; inner max.edge options to compare
outer_edge <- 30 # km; max.edge in the outer extension

# data and grids (as in 01-index-generation.R) -------------------------------

df <- readRDS("data-raw/wrangled-hbll-dog-sets.rds") |>
  drop_na(catch_count, offset, julian, soak, lglsp_hook_count)

grid_hbll <- rbind(gfplot::hbll_inside_n_grid$grid, gfplot::hbll_inside_s_grid$grid) |>
  add_utm_columns(utm_crs = 32609, ll_names = c("X", "Y"),
    utm_names = c("UTM.lon", "UTM.lat"))
grid_dog <- gfplot::dogfish_grid$grid |>
  mutate(UTM.lon = X / 1e3, UTM.lat = Y / 1e3) |>
  distinct()

configs <- list(
  hblldog = list(data = df, grid = grid_hbll),
  `hbll-n-s` = list(
    data = filter(df, survey_abbrev %in% c("HBLL INS S", "HBLL INS N")),
    grid = grid_hbll),
  `dog-j-circle` = list(
    data = filter(df, survey_sep %in% c("dog", "dog comp", "dog-jhook")),
    grid = grid_dog)
)

# coastline in UTM km --------------------------------------------------------

coast <- rnaturalearth::ne_countries(scale = 10,
  country = c("Canada", "United States of America"), returnclass = "sf") |>
  sf::st_geometry() |>
  sf::st_transform(32609) |>
  sf::st_crop(sf::st_bbox(c(xmin = 450e3, xmax = 1050e3,
    ymin = 5300e3, ymax = 5750e3), crs = sf::st_crs(32609)))
coast <- sf::st_set_crs(coast, NA) / 1000

# plot -----------------------------------------------------------------------

plot_mesh <- function(mesh, d, title, lims) {
  ggplot() +
    geom_sf(data = coast, fill = "grey85", colour = NA) +
    inlabru::gg(mesh$mesh, edge.color = "grey40", edge.linewidth = 0.15,
      int.color = "blue", ext.color = "red", int.linewidth = 0.5,
      ext.linewidth = 0.5) +
    geom_point(aes(UTM.lon, UTM.lat), data = d, size = 0.25,
      colour = "darkorange", alpha = 0.7) +
    coord_sf(xlim = lims$x, ylim = lims$y, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    ggsidekick::theme_sleek() +
    theme(plot.title = element_text(size = 9),
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA))
}

stats <- list()
for (m in names(configs)) {
  cfg <- configs[[m]]
  meshes <- c(
    list(current = make_mesh(cfg$data, c("UTM.lon", "UTM.lat"), cutoff = 5)),
    setNames(lapply(max_edges, \(me) make_survey_mesh(cfg$data, cfg$grid,
      max_edge = me, outer_edge = outer_edge)), paste0("max_edge_", max_edges))
  )
  s <- imap(meshes, \(x, nm) mutate(mesh_edge_summary(x, cfg$data),
    model = m, mesh = nm, .before = 1)) |> bind_rows()
  stats[[m]] <- s

  # zoom to the inner domain (data + grid) plus a margin
  xy <- rbind(cfg$data[, c("UTM.lon", "UTM.lat")], cfg$grid[, c("UTM.lon", "UTM.lat")])
  lims <- list(x = range(xy$UTM.lon) + c(-20, 20), y = range(xy$UTM.lat) + c(-20, 20))

  titles <- c(current = "current: make_mesh(cutoff = 5)",
    setNames(paste0("max_edge = c(", max_edges, ", ", outer_edge, "), cutoff = ",
      round(max_edges / 3, 1)), paste0("max_edge_", max_edges)))
  titles <- paste0(titles, "\n", s$n_nodes, " nodes; median edge at data ",
    round(s$edge_median, 1), " km (max ", round(s$edge_max, 1), ")")

  g <- cowplot::plot_grid(plotlist = map2(meshes, titles,
    \(x, t) plot_mesh(x, cfg$data, t, lims)), ncol = 2)
  ggsave(paste0("figures/mesh-options-", m, ".png"), g, width = 10, height = 8,
    bg = "white")
}

stats <- bind_rows(stats)
print(as.data.frame(stats), digits = 3)
