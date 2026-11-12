# Mesh construction for the spatiotemporal index models.
#
# Mesh resolution is set from the smallest spatial range we want to resolve
# (inner max.edge ~ range / 3-5) rather than from set spacing. The inner
# boundary is a non-convex hull around the sets and the prediction grid, so
# the grid always falls in the fine part of the mesh, and a coarse outer
# extension about as wide as the largest plausible range limits SPDE boundary
# effects. The outer triangles only carry that buffer, so they can be much
# coarser (outer_edge) than the inner ones. Coordinates are UTM km.

make_survey_mesh <- function(d, grid, max_edge = 6, outer_edge = 30,
                             outer = 40, convex = 4,
                             xy_cols = c("UTM.lon", "UTM.lat")) {
  xy <- as.matrix(d[, xy_cols])
  gx <- unique(as.matrix(grid[, xy_cols])) # grids are often year-replicated
  inner <- fmesher::fm_nonconvex_hull(rbind(xy, gx), convex = convex,
    resolution = c(200, 200))
  m <- fmesher::fm_mesh_2d_inla(
    loc = xy, boundary = list(inner, NULL),
    max.edge = c(max_edge, outer_edge), # fine inside, coarse outside
    offset = c(-0.01, outer), # outer extension (km)
    cutoff = max_edge / 3
  )
  sdmTMB::make_mesh(d, xy_cols, mesh = m)
}

# Edge lengths of the triangles that contain data, as a check that the mesh
# is fine enough where it matters.
mesh_edge_summary <- function(mesh, d, xy_cols = c("UTM.lon", "UTM.lat")) {
  m <- mesh$mesh
  loc <- m$loc[, 1:2]
  tv <- m$graph$tv
  tri <- unique(fmesher::fm_bary(m, as.matrix(d[, xy_cols]))$index)
  e <- c(
    sqrt(rowSums((loc[tv[tri, 1], ] - loc[tv[tri, 2], ])^2)),
    sqrt(rowSums((loc[tv[tri, 2], ] - loc[tv[tri, 3], ])^2)),
    sqrt(rowSums((loc[tv[tri, 3], ] - loc[tv[tri, 1], ])^2))
  )
  data.frame(n_nodes = m$n, edge_median = median(e),
    edge_q90 = unname(quantile(e, 0.9)), edge_max = max(e))
}

