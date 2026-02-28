
#' Remap an object to the ROMS grid. 
#' 
#' Find the nearest-neighbour coordinates of `x` in the coordinate arrays of `coords`. 
#' 
#' The input `coords` is assumed to be a 2-layer SpatRaster
#' and using [FNN::get.knnx()] the nearest matching position of the coordinates
#' of `x` is found in the grid space of `coords`. The motivating use-case is the
#' curvilinear longitude and latitude arrays of ROMS model output.
#' 
#' No account is made for the details of a ROMS cell, though this may be included in future. We tested only with the "lon_u" and "lat_u"
#' arrays. 
#' @param x object to transform to the grid space, e.g. an sf object, SpatVector, matrix of coordinates, or a wk-handleable geometry
#' @param coords romscoords SpatRaster
#' @param crop logical, if `TRUE` crop x to the extent of the boundary of the values in coords
#' @param lonlat logical, if `TRUE` check for need to back-transform to longitude/latitude and do it
#' @param ... unused
#' @note Do not use this for extraction purposes without checking the output, this is best used for exploration
#' and visualization. Re-mapping ROMS data is better done by looking up the `coords_points` within spatial objects, 
#' and transferring via the grid index. 
#' @return data.frame with columns `x_` and `y_` in grid index space, plus any feature identifiers
#' @export
#' @examples 
#' ant_idx <- romsmap(antarctica, ice_coords)
#' plot(ice_fake, col = c("transparent", hcl.colors(24)), asp = 1)
#' ## ant_idx is a data.frame of index-space coordinates
#' ## use wk to reconstruct if needed
romsmap <- function(x, ...) {
  UseMethod("romsmap")
}

#' @rdname romsmap
#' @export
romsmap.SpatVector <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  romsmap(wk::as_wk_wkb(x), coords = coords, crop = crop, lonlat = lonlat, ...)
}

#' @rdname romsmap
#' @export
romsmap.sf <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  romsmap(wk::as_wk_wkb(x), coords = coords, crop = crop, lonlat = lonlat, ...)
}

#' @rdname romsmap
#' @export
romsmap.sfc <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  romsmap(wk::as_wk_wkb(x), coords = coords, crop = crop, lonlat = lonlat, ...)
}

#' @rdname romsmap
#' @export
romsmap.SpatialPolygonsDataFrame <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  romsmap(wk::as_wk_wkb(x), coords = coords, crop = crop, lonlat = lonlat, ...)
}

#' @rdname romsmap
#' @export
romsmap.SpatialLinesDataFrame <- romsmap.SpatialPolygonsDataFrame

#' @rdname romsmap
#' @export
romsmap.SpatialPointsDataFrame <- romsmap.SpatialPolygonsDataFrame


#' @rdname romsmap
#' @export
romsmap.default <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    if (ncol(x) < 2) stop("need at least 2 columns for coordinates")
    ## treat as lon,lat points
    xy <- x[, 1:2, drop = FALSE]
  } else {
    ## try the wk path for anything else
    x <- wk::as_wk_wkb(x)
    return(romsmap(x, coords = coords, crop = crop, lonlat = lonlat, ...))
  }
  
  coord_xy <- cbind(terra::values(coords[[1]]), terra::values(coords[[2]]))
  kd <- FNN::get.knnx(coord_xy, xy, k = 1)
  index <- arrayInd(kd$nn.index[, 1L], .dim = c(nrow(coords), ncol(coords)))
  ## flip row index because terra is top-down, ROMS index space is bottom-up
  data.frame(x_ = index[, 2L], y_ = nrow(coords) - index[, 1L])
}


#' @rdname romsmap
#' @export
romsmap.wk_wkb <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  ## extract all coordinates from the geometry
  crd <- wk::wk_coords(x)
  xy <- cbind(crd$x, crd$y)
  
  ## look up nearest grid cell for each vertex
  coord_xy <- cbind(terra::values(coords[[1]]), terra::values(coords[[2]]))
  kd <- FNN::get.knnx(coord_xy, xy, k = 1)
  index <- arrayInd(kd$nn.index[, 1L], .dim = c(nrow(coords), ncol(coords)))
  
  crd$x <- index[, 2L]
  crd$y <- nrow(coords) - index[, 1L]
  crd
}

#' @rdname romsmap
#' @export
romsmap.wk_wkt <- romsmap.wk_wkb


#' Boundary polygon from raster of coordinates. 
#' 
#' Create a boundary polygon by tracking around coordinates stored in a SpatRaster. 
#' 
#' The first layer is treated as the X coordinate, second as Y. 
#' @param cds two-layer SpatRaster
#'
#' @export
#' @examples 
#' ice_grid_boundary <- romsboundary(ice_coords)
#' ## boundary is a wk_polygon in lon/lat space
#' ice_grid_boundary
romsboundary <- function(cds) {
  nc <- ncol(cds)
  nr <- nrow(cds)
  
  left   <- terra::cellFromCol(cds, 1)
  bottom <- terra::cellFromRow(cds, nr)
  right  <- rev(terra::cellFromCol(cds, nc))
  top    <- rev(terra::cellFromRow(cds, 1))
  
  cells <- unique(c(left, bottom, right, top))
  
  xvals <- terra::values(cds[[1]])[cells]
  yvals <- terra::values(cds[[2]])[cells]
  
  ## close the ring
  coords <- cbind(c(xvals, xvals[1]), c(yvals, yvals[1]))
  wk::wk_polygon(wk::xy(coords[, 1], coords[, 2]), feature_id = rep(1L, nrow(coords)))
}

#' Create SpatialPoints from coordinate arrays (deprecated)
#'
#' Use `coords_points` to get a matrix of coordinates from a two-layer coord SpatRaster. Returns a 
#' two column matrix with columns x (longitude) and y (latitude).
#' @param x two layer SpatRaster with longitude and latitude values
#' @param ... ignored
#'
#' @return matrix of coordinates with columns x, y
#' @export
#' @examples
#' pts <- coords_points(ice_coords)
#' head(pts)
coords_points <- function(x, ...) {
  cbind(x = terra::values(x[[1]]), y = terra::values(x[[2]]))
}


#' @param x a thing
#' @param mask logical
#' @param ... dots
#' 
#' @name romsboundary
#' @export
databoundary <- function(x, mask = NULL, ...) {
  stop("databoundary is defunct")
}

## put any raster into xy-index space (0, nc, 0, nr)
set_indextent <- function(x) {
  terra::ext(x) <- terra::ext(0, ncol(x), 0, nrow(x))
  x
}

#' Crop a ROMS layer
#' 
#' Crop a ROMS data layer from `romsdata` with a raster extent. 
#' 
#' The spatial crop is performed in the coordinate space of roms data. 
#' @param x ROMS xy- coordinates, see `romscoords`
#' @param ext `terra::ext` in the coordinate system of `x`
#' @param ... ignored
#'
#' @export
#' @examples
#' ## notice that extent is in long-lat, but ice_local is in the grid
#' ## space of ice_coords
#' ice_local <- croproms(ice_coords, terra::ext(100, 120, -75, -60))
#' plot(ice_coords[[2]], col = grey(seq(0, 1, length  = 20)))
#' plot(terra::crop(ice_fake, ice_local), add = TRUE)
croproms <- function(x, ext, ...) {
  xy <- cbind(terra::values(x[[1]]), terra::values(x[[2]]))
  incells <- which(xy[, 1] >= ext[1] & xy[, 1] <= ext[2] &
                     xy[, 2] >= ext[3] & xy[, 2] <= ext[4])
  x1 <- x[[1]]
  terra::values(x1) <- NA_real_
  x1[incells] <- 0
  terra::ext(terra::trim(x1))
}
