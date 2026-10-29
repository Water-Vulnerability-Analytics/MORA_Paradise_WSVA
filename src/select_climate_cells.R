#' Keep the two GridMET cells selected for the Paradise WBM
#'
#' The cell search is documented in Reports/Explore_MORA_Climate_Grid.Rmd.
#' Six GridMET cells cover the buffered contributing area. The three southern
#' cells were dropped after comparing elevation and land cover, leaving the
#' cell containing the watershed and its western and eastern neighbors. The
#' eastern cell tracked cumulative precipitation like the other two and the
#' Paradise SNOTEL station, but sits lower, so more precipitation falls as rain
#' and its cumulative SWE runs well below the others; it was also dropped.
#' The retained cells are the one containing the watershed and the one to the
#' west.
#'
#' @param df Data frame or sf object with GridMET cell-center coordinates in
#'   columns `x` (longitude) and `y` (latitude), EPSG:4326.
#' @return `df` filtered to the two selected cells.
select_climate_cells <- function(df) {
  df %>%
    dplyr::filter(y > 46.81,    # drops the three southern cells
                  x < -121.7)   # drops the eastern cell
}
