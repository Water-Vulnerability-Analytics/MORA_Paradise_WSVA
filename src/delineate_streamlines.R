#' Delineate Stream Networks using WhiteboxTools and terra
#'
#' @param dem SpatRaster. Input Digital Elevation Model. Can be omitted if resume = TRUE and intermediate DEM exists on disk.
#' @param aoi sf or SpatVector object (optional). Area of interest/watershed polygon to restrict processing and stream output.
#' @param output_dir Character. Directory to save files and look for intermediate cached inputs.
#' @param file_prefix Character. Prefix for generated file names.
#' @param threshold Numeric. Minimum flow accumulation value (in cells) to initiate a stream channel.
#' @param breach_dist Numeric. Maximum breach distance for depression breaching (least cost).
#' @param aoi_buffer_dist Numeric. Buffer distance (in CRS units) applied around AOI for flow processing to avoid edge effects.
#' @param mainline_only Logical. If TRUE, isolates the continuous primary mainstem channel
#'   using a graph shortest-path (headwater to outlet) via
#'   \code{\link{extract_mainstem_path}} -- see that function's docs for why a naive
#'   line-merge approach was replaced. Default FALSE.
#' @param resume Logical. If TRUE, looks for existing processed DEM raster in output_dir and skips reprojection/cropping. Default TRUE.
#' @param output_option Character. Controls disk export behavior: "all" (saves all intermediate rasters & vector), "flowlines_only" (saves vector flowlines and deletes intermediate rasters), or "none" (returns sf object without keeping files on disk).
#' @param write_gpkg Logical. If TRUE and saving files, exports vector streams as GeoPackage.
#' @param overwrite Logical. If TRUE, overwrites stream output files.
#' @param crs_out Character. EPSG code or PROJ string for output vector projection (default EPSG:5070).
#'
#' @return An sf object containing the vector streamlines.
#' @export
delineate_streamlines <- function(dem = NULL,
                                  aoi = NULL,
                                  output_dir,
                                  file_prefix = "streams",
                                  threshold = 1000,
                                  breach_dist = 10,
                                  aoi_buffer_dist = 500,
                                  mainline_only = FALSE,
                                  resume = TRUE,
                                  output_option = c("all", "flowlines_only", "none"),
                                  write_gpkg = TRUE,
                                  overwrite = TRUE,
                                  crs_out = "EPSG:5070") {
  
  output_option <- match.arg(output_option)
  
  if (is.null(output_dir)) {
    stop("An 'output_dir' must be specified to store or read intermediate inputs.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # File paths
  dem_path         <- file.path(output_dir, paste0(file_prefix, "_dem.tif"))
  breached_path    <- file.path(output_dir, paste0(file_prefix, "_breached.tif"))
  flow_dir_path    <- file.path(output_dir, paste0(file_prefix, "_flow_dir.tif"))
  flow_accum_path  <- file.path(output_dir, paste0(file_prefix, "_flow_accum.tif"))
  stream_raster_p  <- file.path(output_dir, paste0(file_prefix, "_streams.tif"))
  stream_vector_p  <- file.path(output_dir, paste0(file_prefix, "_streams.shp"))
  
  # --------------------------------------------------------------------------
  # 1. DEM Import & Preprocessing Logic
  # --------------------------------------------------------------------------
  if (resume && file.exists(dem_path)) {
    message("--> Resuming from saved DEM: ", dem_path)
    dem <- terra::rast(dem_path)
  } else {
    if (is.null(dem)) {
      stop("No 'dem' provided and no existing DEM raster found at: ", dem_path)
    }
    
    message("--> Processing original DEM (reprojecting/cropping)...")
    
    if (!terra::same.crs(dem, crs_out)) {
      dem <- terra::project(dem, crs_out)
    }
    
    if (!is.null(aoi)) {
      aoi_v <- terra::vect(aoi) %>% 
        terra::project(terra::crs(dem))
      
      aoi_buf <- terra::buffer(aoi_v, width = aoi_buffer_dist)
      
      dem <- terra::crop(dem, aoi_buf) %>% 
        terra::mask(aoi_buf)
    }
    
    terra::writeRaster(dem, dem_path, overwrite = TRUE)
  }
  
  # --------------------------------------------------------------------------
  # 2. Whitebox Hydro Processing Workflow
  # --------------------------------------------------------------------------
  if (!resume || !file.exists(breached_path)) {
    message("--> Breaching depressions...")
    whitebox::wbt_breach_depressions_least_cost(
      dem = dem_path,
      output = breached_path,
      dist = breach_dist
    )
  }
  
  if (!resume || !file.exists(flow_dir_path) || !file.exists(flow_accum_path)) {
    message("--> Computing D8 flow pointer and accumulation...")
    whitebox::wbt_d8_pointer(
      dem = breached_path,
      output = flow_dir_path
    )
    
    whitebox::wbt_d8_flow_accumulation(
      input = breached_path,
      output = flow_accum_path,
      out_type = "cells"
    )
  }
  
  # --------------------------------------------------------------------------
  # 3. Stream Channel Extraction
  # --------------------------------------------------------------------------
  message("--> Extracting streams with cell threshold = ", threshold)
  whitebox::wbt_extract_streams(
    flow_accum = flow_accum_path,
    output = stream_raster_p,
    threshold = threshold
  )
  
  whitebox::wbt_raster_streams_to_vector(
    streams = stream_raster_p,
    d8_pntr = flow_dir_path,
    output = stream_vector_p
  )
  
  # Read vector streams into R
  streams_sf <- sf::st_read(stream_vector_p, quiet = TRUE) %>%
    sf::st_set_crs(terra::crs(dem)) %>%
    sf::st_transform(crs_out)
  
  # Clip precisely to exact AOI polygon boundary
  if (!is.null(aoi)) {
    aoi_sf <- sf::st_as_sf(aoi) %>% sf::st_transform(crs_out)
    streams_sf <- sf::st_intersection(streams_sf, sf::st_union(aoi_sf))
  }
  
  # --------------------------------------------------------------------------
  # 4. Mainline Filtering (If mainline_only = TRUE)
  # --------------------------------------------------------------------------
  if (mainline_only) {
    message("--> Isolating longest connected mainstem flowpath (graph shortest-path)...")
    
    # Graph shortest-path mainstem extraction (see src/extract_mainstem_path.R).
    # Replaces a naive st_combine()+st_line_merge()+pick-longest-piece approach,
    # which only merges through degree-2 nodes and therefore fractures -- and
    # picks the wrong "longest" piece -- at every confluence (tributary
    # junction). `dem` here is still the unbreached, (re)projected raster from
    # step 1 above -- whitebox's breach/flow-accum/streams outputs are written
    # to disk by file path and never reassigned into this R object.
    mainstem_result <- extract_mainstem_path(streams_sf, dem)
    mainstem_line <- sf::st_transform(mainstem_result$mainstem_line, crs_out)
    
    streams_sf <- sf::st_sf(
      geometry = mainstem_line,
      crs = crs_out
    ) %>% 
      dplyr::mutate(segment_id = 1, length_m = mainstem_result$length_m)
  }
  
  # --------------------------------------------------------------------------
  # 5. Cleanup & File Saving Logic
  # --------------------------------------------------------------------------
  if (output_option == "flowlines_only") {
    rasters_to_remove <- c(dem_path, breached_path, flow_dir_path, flow_accum_path, stream_raster_p)
    file.remove(rasters_to_remove[file.exists(rasters_to_remove)])
  }
  
  if (output_option %in% c("all", "flowlines_only")) {
    if (write_gpkg) {
      gpkg_path <- file.path(output_dir, paste0(file_prefix, "_streams.gpkg"))
      sf::st_write(
        obj = streams_sf,
        dsn = gpkg_path,
        delete_dsn = overwrite,
        quiet = TRUE
      )
      shp_files <- list.files(output_dir, pattern = paste0("^", file_prefix, "_streams\\."), full.names = TRUE)
      shp_files <- shp_files[!grepl("\\.gpkg$", shp_files)]
      file.remove(shp_files)
    }
  }
  
  return(streams_sf)
}
