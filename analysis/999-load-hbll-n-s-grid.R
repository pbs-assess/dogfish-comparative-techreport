# HBLL inside north + south prediction grid with depth strata, from gfdata.
# Sourced by analysis/sdmTMB combined models/01b-index-depths.R.
#
# Produces `grid`: one row per active HBLL INS N/S survey block with
#   UTM.lon, UTM.lat  block centroid, UTM zone 9 (EPSG 32609), km
#   area              block area (km^2), clipped to land
#   depth_m           mean block depth from GFBio
#   depth_bin         1 = 40-70 m stratum, 2 = 71-100 m stratum (from the
#                     block's GFBio grouping_code, matching the data's strata)
#   survey_abbrev, grouping_code, block_id
#
# depth_bin uses the survey design strata rather than block depth_m because
# block depth_m spans 0-586 m and doesn't line up with the strata.

grid <- gfdata::load_survey_blocks("hbll_ins", type = "XY", active_only = TRUE) |>
  dplyr::filter(!is.na(grouping_code)) |>
  dplyr::mutate(
    # odd codes 279, 281, ... 297 = 40-70 m; even codes 280, ... 298 = 71-100 m
    depth_bin = ifelse(grouping_code %% 2 == 1, 1, 2),
    UTM.lon = X,
    UTM.lat = Y
  ) |>
  dplyr::select(survey_abbrev, block_id, grouping_code, depth_bin, depth_m, area, UTM.lon, UTM.lat) |>
  as.data.frame()
