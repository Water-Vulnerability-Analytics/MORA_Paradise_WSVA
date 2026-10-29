#' Extract Continuous Mainstem Flowline and Longitudinal Profile
#'
#' @description
#' Traverses a spatial stream vector network to extract a single continuous mainstem 
#' flowline from the maximum headwater elevation down to the basin outlet. The method 
#' converts vector flowlines into a spatial network graph, identifies the global 
#' elevation extremes using a coupled Digital Elevation Model (DEM), and solves 
#' Dijkstra's shortest topological path between the endpoints. The resulting mainstem 
#' is sampled to calculate downstream cumulative distance, elevation, and channel slope.
#'
#' @details
#' Standard spatial linestring merging (`sf::st_line_merge`) typically fails across complex 
#' river networks due to topological degree branch conflicts (confluences with degree > 2) 
#' and micro-gaps or digitization direction mismatches. An earlier version of this function
#' used that naive merge-and-pick-longest-piece approach despite this docstring already
#' describing the graph method below, which meant it silently returned the longest reach
#' BETWEEN confluences rather than the true mainstem. Steps 1-4 below are now actually
#' implemented in \code{\link{extract_mainstem_path}} (\code{src/extract_mainstem_path.R}),
#' which this function calls; step 5 remains here.
#' 
#' This function resolves network connectivity issues using a multi-step graph approach:
#' \enumerate{
#'   \item \strong{Network Cleaning & Graph Construction}: Merges input vector flowlines, 
#'         casts them to uniform \code{LINESTRING} geometries, and constructs an undirected 
#'         spatial network using \code{\link[sfnetworks]{as_sfnetwork}}.
#'   \item \strong{Connected Component Isolation}: Filters the network to retain only the 
#'         largest connected subgraph component using \code{\link[tidygraph]{group_components}}
#'         (a \code{\link[igraph]{components}} wrapper), removing isolated tributary segments.
#'   \item \strong{DEM Elevation Extremes Anchor}: Extracts surface elevations from the 
#'         provided DEM at all graph nodes. Identifies the global \emph{headwater} node 
#'         (maximum elevation) and \emph{outlet} node (minimum elevation).
#'   \item \strong{Graph Path Traversal}: Computes the edge-length weighted shortest path 
#'         from headwater to outlet using Dijkstra's algorithm (\code{\link[sfnetworks]{st_network_paths}}).
#'   \item \strong{Continuous Longitudinal Profiling}: Stitches path edges into a unified 
#'         linestring, samples nodes at regular distance intervals, enforces headwater-to-outlet 
#'         orientation, and calculates local gradient (\eqn{dz/dx}).
#' }
#'
#' @param flowlines An \code{\link[sf]{sf}} object containing 
#'   \code{LINESTRING} or \code{MULTILINESTRING} stream network geometries.
#' @param dem A \code{\link[terra]{SpatRaster}} digital elevation model covering the 
#'   extent of the stream network. Used for node elevation extraction and endpoint identification.
#' @param sample_dist Numeric. Spacing distance in meters between profile sampling points 
#'   along the extracted mainstem. Default is \code{10}.
#' @param min_reach_length Numeric. Minimum required length in meters for the extracted 
#'   mainstem line. Throws an error if the path length falls below this value. Default is \code{50}.
#'
#' @return A named \code{\link[base]{list}} containing three components:
#' \describe{
#'   \item{\code{nodes}}{An \code{\link[sf]{sf}} \code{POINT} dataset containing sampled 
#'     profile points with attributes: \code{node_id}, \code{cum_dist_m}, \code{elevation_m}, 
#'     \code{dz}, \code{dx}, \code{slope_m_per_m}, and \code{slope_pct}.}
#'   \item{\code{mainstem_flowline}}{An \code{\link[sf]{sf}} \code{LINESTRING} feature representing 
#'     the stitched continuous mainstem channel.}
#'   \item{\code{profile_plot}}{A \code{\link[ggplot2]{ggplot}} object depicting the longitudinal 
#'     elevation profile colored by channel slope percentage.}
#' }
#'
#' @references
#' Dijkstra, E. W. (1959). A note on two problems in connexion with graphs. 
#' \emph{Numerische Mathematik}, 1(1), 269–271. \doi{10.1007/BF01386390}
#' 
#' van der Walt, V., & Morgan, D. J. (2020). Spatial Networks in R: The sfnetworks Package. 
#' \emph{R Journal}, 12(2), 120–139.
#' 
#' Hack, J. T. (1957). Studies of longitudinal stream profiles in Virginia and Maryland. 
#' \emph{U.S. Geological Survey Professional Paper}, 294-B, 45–97.
#'
#' @seealso \code{\link{extract_mainstem_path}}, \code{\link[sfnetworks]{as_sfnetwork}}, 
#'   \code{\link[sfnetworks]{st_network_paths}}, \code{\link[terra]{extract}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # Load stream vector network and elevation raster
#' streams <- sf::st_read("flowlines.gpkg")
#' dem <- terra::rast("elevation.tif")
#'
#' # Extract mainstem profile sampled every 25 meters
#' profile_results <- extract_stream_profile(
#'   flowlines = streams,
#'   dem = dem,
#'   sample_dist = 25,
#'   min_reach_length = 100
#' )
#'
#' # Plot longitudinal profile
#' print(profile_results$profile_plot)
#'
#' # Access sampled nodes
#' head(profile_results$nodes)
#' }
extract_stream_profile <- function(flowlines, dem, sample_dist = 10, min_reach_length = 50) {
  
  # 1. Extract the headwater-to-outlet mainstem via graph shortest path
  #    (see src/extract_mainstem_path.R). This replaces a naive
  #    st_combine()+st_line_merge()+pick-longest-piece approach, which only
  #    merges through degree-2 nodes and therefore fractures -- and picks the
  #    wrong "longest" piece -- at every confluence. See that file's roxygen
  #    @details for the full explanation. extract_mainstem_path() handles the
  #    flowlines/dem CRS match internally.
  mainstem_result <- extract_mainstem_path(flowlines, dem)
  mainstem_line <- mainstem_result$mainstem_line
  
  if (mainstem_result$length_m < min_reach_length) {
    stop("Extracted mainstem path is shorter than min_reach_length threshold (",
         round(mainstem_result$length_m, 1), " m < ", min_reach_length, " m).")
  }
  
  # 2. Sample nodes along the continuous mainstem LINESTRING
  pts <- sf::st_line_sample(mainstem_line, density = 1 / sample_dist) %>%
    sf::st_cast("POINT") %>%
    sf::st_as_sf()
  
  # 3. Calculate continuous cumulative distance (meters)
  coords <- sf::st_coordinates(pts)
  dist_intervals <- c(0, sqrt(diff(coords[, 1])^2 + diff(coords[, 2])^2))
  pts$cum_dist_m <- cumsum(dist_intervals)
  
  # 4. Extract DEM elevation at each node
  pts$elevation_m <- terra::extract(dem, terra::vect(pts))[[2]]
  
  # 5. Orient from Headwater (high elevation) to Outlet (low elevation)
  # If the line was digitized downstream-to-upstream, reverse distance so 0 = headwater
  if (pts$elevation_m[1] < pts$elevation_m[nrow(pts)]) {
    pts$cum_dist_m <- max(pts$cum_dist_m) - pts$cum_dist_m
    pts <- pts[order(pts$cum_dist_m), ]
  }
  
  # 6. Calculate Reach Slope (m/m and %) along mainstem
  nodes <- pts %>%
    dplyr::mutate(
      node_id = dplyr::row_number(),
      dz = elevation_m - dplyr::lead(elevation_m),
      dx = dplyr::lead(cum_dist_m) - cum_dist_m,
      slope_m_per_m = dz / dx,
      slope_pct = slope_m_per_m * 100
    )
  
  # 7. Generate Clean Profile Plot
  p <- ggplot2::ggplot(nodes, ggplot2::aes(x = cum_dist_m, y = elevation_m)) +
    ggplot2::geom_line(color = "royalblue", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(color = slope_pct), size = 1.8) +
    ggplot2::scale_color_viridis_c(option = "magma", name = "Slope (%)", na.value = "grey50") +
    ggplot2::labs(
      title = "Mainstem Longitudinal Profile",
      subtitle = paste0("Total Length: ", round(max(nodes$cum_dist_m), 1), " m"),
      x = "Cumulative Distance Downstream (m)",
      y = "Elevation (m)"
    ) +
    ggplot2::theme_minimal()
  
  return(list(nodes = nodes, mainstem_flowline = mainstem_line, profile_plot = p))
}
