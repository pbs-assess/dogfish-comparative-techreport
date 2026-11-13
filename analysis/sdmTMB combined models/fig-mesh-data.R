library(cowplot)

hblldog <- readRDS(paste0("figures/sog-model-mesh-", "hblldog", ".rds"))
hbll <- readRDS(paste0("figures/sog-model-mesh-", "hbll-n-s", ".rds"))
dog <- readRDS (paste0("figures/sog-model-mesh-", "dog-j-circle", ".rds"))

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

cv <- cowplot::plot_grid(
  hblldog, hbll,
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 1)
)

cv2 <- cowplot::plot_grid(
  cv, dog,
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 0.56)
)
cv2

ggplot2::ggsave("figures/mesh_points.png", cv2, height = 10, width = 7, dpi = 200, bg = "white")
