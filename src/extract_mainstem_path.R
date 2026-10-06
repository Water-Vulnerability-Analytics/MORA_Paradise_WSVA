#' Extract the Longest Connected Headwater-to-Outlet Mainstem Path from a Stream Network
#'
#' @description
#' Given a (possibly branching) vector stream network and a DEM, builds a spatial
#' network graph from the flowlines, restricts it to the single largest connected
#' component, identifies the highest- and lowest-elevation nodes within that
#' component (headwater and outlet), and returns the edge-length-weighted
#' shortest path (Dijkstra's algorithm) between them as a single continuous
#' LINESTRING.
#'
#' @details
#' This replaces a naive `sf::st_combine() + sf::st_line_merge() + pick-longest-
#' piece-after-st_cast()` approach previously used in both
#' `extract_stream_profile()` and `delineate_streamlines(mainline_only = TRUE)`.
#' `st_line_merge()` only stitches lines together through nodes of degree 2
#' (simple pass-through points); at any confluence (degree > 2 -- i.e. anywhere
#' a tributary joins the network) it cannot decide which branch continues the
#' line, so it leaves the network split there. `st_cast(..., "LINESTRING")`
#' then turns every one of those splits into a separate feature, and picking
#' the longest one just picks the longest uninterrupted run BETWEEN two
#' confluences -- frequently a short tributary or a short mid-network reach,
#' not the true headwater-to-outlet mainstem. This is why both functions were
#' observed clipping short segments instead of the full drainage.
#'
#' Using a graph shortest path instead correctly traverses every confluence
#' regardless of how many branches meet there, because Dijkstra's algorithm
#' operates on the full graph topology rather than relying on geometric line
#' merging.
#'
#' Note: this assumes the input `flowlines` network is vertex-exact at
#' confluences, which is true for whitebox-derived raster-to-vector stream
#' networks -- the only source currently feeding this function in this
#' codebase (via `delineate_streamlines()`'s `wbt_raster_streams_to_vector()`
#' output). A network with small digitizing gaps at junctions (e.g. some raw
#' NHD products) may need vertex snapping (e.g. `sf::st_snap()`) or
#' `sfnetworks::to_spatial_subdivision()` applied first -- not added here
#' since it isn't needed for the current data source and hasn't been tested
#' against this codebase's data. `stitch_path_edges()` below will emit a
#' warning if it detects a gap larger than `snap_tol`, which is meant to catch
#' exactly this situation if the function is ever pointed at a different
#' flowline source.
#'
#' @param flowlines sf LINESTRING/MULTILINESTRING network.
#' @param dem terra SpatRaster covering the flowline extent, used to identify
#'   the headwater (max elevation) and outlet (min elevation) nodes.
#'   `flowlines` is reprojected internally to match `dem`'s CRS if needed.
#' @param snap_tol Numeric, CRS units (typically meters). Diagnostic tolerance
#'   used when stitching path edges into a single linestring -- if consecutive
#'   edges in the shortest path don't share an endpoint within this tolerance,
#'   a warning is issued (see Details). Default 1.
#'
#' @return A list with `mainstem_line` (sfc LINESTRING, single feature, CRS =
#'   `dem`'s CRS) and `length_m` (numeric, meters).
#'
#' @seealso extract_stream_profile, delineate_streamlines
#' @export
extract_mainstem_path <- function(flowlines, dem, snap_tol = 1) {

  if (!terra::same.crs(flowlines, dem)) {
    flowlines <- sf::st_transform(flowlines, terra::crs(dem))
  }

  # Normalize to simple LINESTRING features before building the network.
  flowlines_clean <- flowlines %>%
    sf::st_geometry() %>%
    sf::st_cast("MULTILINESTRING") %>%
    sf::st_cast("LINESTRING") %>%
    sf::st_as_sf()

  # 1. Build an undirected spatial network graph from the flowlines.
  net <- sfnetworks::as_sfnetwork(flowlines_clean, directed = FALSE)

  # 2. Restrict to the single largest connected component -- drops isolated/
  #    disconnected tributary fragments that can't reach the rest of the
  #    network (and that would otherwise break the headwater/outlet search
  #    below if the DEM's global max/min happened to fall on one of them).
  net <- net %>%
    tidygraph::activate("nodes") %>%
    dplyr::mutate(.component = tidygraph::group_components()) %>%
    dplyr::filter(.component == names(which.max(table(.component)))) %>%
    dplyr::select(-".component")

  # 3. Edge weights = geographic edge length, for Dijkstra's algorithm below.
  net <- net %>%
    tidygraph::activate("edges") %>%
    dplyr::mutate(.weight = sfnetworks::edge_length())

  # 4. Identify headwater (max elevation) and outlet (min elevation) nodes
  #    within this connected component, from the DEM.
  node_pts <- net %>% tidygraph::activate("nodes") %>% sf::st_as_sf()
  node_elev <- terra::extract(dem, terra::vect(node_pts))[[2]]

  if (all(is.na(node_elev))) {
    stop("extract_mainstem_path(): could not extract any DEM elevation at ",
         "network nodes -- check that `dem` covers the extent of `flowlines`.")
  }

  headwater_idx <- which.max(ifelse(is.na(node_elev), -Inf, node_elev))
  outlet_idx    <- which.min(ifelse(is.na(node_elev),  Inf, node_elev))

  if (headwater_idx == outlet_idx) {
    stop("extract_mainstem_path(): headwater and outlet nodes resolved to the ",
         "same network node -- check DEM coverage/resolution relative to the ",
         "flowline network.")
  }

  # 5. Dijkstra's shortest (edge-length-weighted) path from headwater to outlet.
  path <- sfnetworks::st_network_paths(
    net,
    from = headwater_idx,
    to = outlet_idx,
    weights = ".weight",
    type = "shortest"
  )

  edge_idx <- path$edge_paths[[1]]
  node_idx <- path$node_paths[[1]]

  if (length(edge_idx) == 0) {
    stop("extract_mainstem_path(): no connected path found between the ",
         "headwater and outlet nodes -- the network may be disconnected ",
         "between the DEM's highest and lowest points within this component.")
  }

  path_edges <- net %>%
    tidygraph::activate("edges") %>%
    sf::st_as_sf() %>%
    dplyr::slice(edge_idx)

  start_coord <- sf::st_coordinates(node_pts[node_idx[1], ])[1, c("X", "Y")]

  mainstem_line <- stitch_path_edges(path_edges, start_coord, tol = snap_tol)
  mainstem_line <- sf::st_sfc(mainstem_line, crs = sf::st_crs(flowlines_clean))

  list(
    mainstem_line = mainstem_line,
    length_m = as.numeric(sf::st_length(mainstem_line))
  )
}


#' Stitch an Ordered Sequence of Network Path Edges into One Continuous LINESTRING
#'
#' @description
#' Internal helper for `extract_mainstem_path()`. Takes the edges of a graph
#' shortest path, already in path order (as returned by
#' `sfnetworks::st_network_paths()`), and walks them into a single ordered
#' LINESTRING -- reversing any edge whose digitized direction runs opposite to
#' the direction of travel, and warning if consecutive edges don't share an
#' endpoint within `tol` (a sign of a small digitizing gap in the source
#' network).
#'
#' @param edges_ordered sf object of LINESTRING edges, one row per path edge,
#'   already in path order (e.g. `edges[path$edge_paths[[1]], ]`).
#' @param start_coord Numeric length-2 vector (X, Y): the coordinate of the
#'   path's starting (headwater) node, used to orient the first edge.
#' @param tol Numeric. Endpoint-matching tolerance in CRS units. Default 1.
#'
#' @return A single sfg LINESTRING.
#' @keywords internal
stitch_path_edges <- function(edges_ordered, start_coord, tol = 1) {

  n <- nrow(edges_ordered)
  coords_list <- vector("list", n)
  current_pt <- as.numeric(start_coord)

  for (i in seq_len(n)) {
    mat <- sf::st_coordinates(sf::st_geometry(edges_ordered)[[i]])[, c("X", "Y"), drop = FALSE]

    d_start <- sqrt(sum((mat[1, ] - current_pt)^2))
    d_end   <- sqrt(sum((mat[nrow(mat), ] - current_pt)^2))

    if (min(d_start, d_end) > tol) {
      warning("stitch_path_edges(): edge ", i, " does not share an endpoint ",
              "with the running path within tol = ", tol, " (closest gap = ",
              round(min(d_start, d_end), 2), " CRS units). Likely a small ",
              "digitizing gap in the source network -- consider snapping ",
              "vertices (e.g. sf::st_snap()) before extraction.")
    }

    if (d_end < d_start) {
      mat <- mat[rev(seq_len(nrow(mat))), , drop = FALSE]
    }

    if (i > 1) {
      mat <- mat[-1, , drop = FALSE]  # drop the vertex shared with the previous edge
    }

    coords_list[[i]] <- mat
    current_pt <- mat[nrow(mat), ]
  }

  all_coords <- do.call(rbind, coords_list)
  sf::st_linestring(all_coords)
}

