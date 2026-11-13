library(cowplot)

hblldog <- readRDS(paste0("figures/sog-model-mesh-", "hblldog", ".rds"))
hbll <- readRDS(paste0("figures/sog-model-mesh-", "hbll-n-s", ".rds"))
dog <- readRDS (paste0("figures/sog-model-mesh-", "dog-j-circle", ".rds"))

# coastline in UTM km (as in fig-mesh-options.R)
coast <- rnaturalearth::ne_countries(scale = 10,
  country = c("Canada", "United States of America"), returnclass = "sf") |>
  sf::st_geometry() |>
  sf::st_transform(32609) |>
  sf::st_crop(sf::st_bbox(c(xmin = 450e3, xmax = 1050e3,
    ymin = 5300e3, ymax = 5750e3), crs = sf::st_crs(32609)))
coast <- sf::st_set_crs(coast, NA) / 1000

# draw the coast beneath the existing layers, keeping the current plot extent
add_coast <- function(g) {
  pp <- ggplot2::ggplot_build(g)$layout$panel_params[[1]]
  # small semi-transparent orange points (as in fig-mesh-options.R)
  g$layers <- lapply(g$layers, function(l) {
    if (inherits(l$geom, "GeomPoint")) {
      ggplot2::geom_point(
        ggplot2::aes(UTM.lon, UTM.lat), data = l$data,
        size = 0.5, colour = "darkorange", alpha = 1, pch = 16
      )
    } else l
  })
  g <- g + ggplot2::geom_sf(data = coast, fill = "grey85", colour = NA, inherit.aes = FALSE)
  # lighten/thin the mesh edges so they don't swamp the data
  g$layers <- lapply(g$layers, function(l) {
    if (identical(l$aes_params$colour, "grey60")) {
      l$aes_params$colour <- "grey40"
      l$aes_params$linewidth <- 0.15
    }
    l
  })
  # order: coast, mesh, then points on top
  n <- length(g$layers)
  is_pt <- vapply(g$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))
  is_sf <- vapply(g$layers, function(l) inherits(l$geom, "GeomSf"), logical(1))
  g$layers <- g$layers[c(which(is_sf), which(!is_pt & !is_sf), which(is_pt))]
  suppressMessages(g + ggplot2::coord_sf(
    xlim = pp$x.range, ylim = pp$y.range, expand = FALSE, datum = NULL
  ))
}
hblldog <- add_coast(hblldog)
hbll <- add_coast(hbll)
dog <- add_coast(dog)

# panel titles are added to the plots themselves (aligned with the panel, which
# is fixed-aspect and so narrower than its grid cell) rather than as cowplot labels
add_title <- function(g, title) {
  g +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title.position = "panel",
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0)
    )
}
hblldog <- add_title(hblldog, "(a) HBLL & SoG dogfish circle/J-hook")
hbll <- add_title(hbll, "(b) HBLL")
dog <- add_title(dog, "(c) SoG dogfish circle/J-hook")

# 2 x 2 grid with a blank bottom-right cell
cv2 <- cowplot::plot_grid(
  hblldog, hbll, dog, NULL,
  ncol = 2,
  nrow = 2
)
cv2

ggplot2::ggsave("figures/mesh_points.png", cv2, height = 7, width = 8.5, dpi = 350, bg = "white")
