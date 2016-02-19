#' Title
#'
#' @param x 
#' @param ... 
#' @param coords romscoords RasterStack
#'
#' @return
#' @export
#'
#' @examples
romsmap <- function(x, coords, ...) {
  UseMethod("romsmap")
}

#' @rdname romsmap
#' @export
#' @importFrom spbabel sptable spFromTable
#' @importFrom nabor knn
#' @importFrom raster intersect as.matrix
romsmap.SpatialPolygonsDataFrame <- function(x, coords, ...) {
  ## first get the intersection
  op <- options(warn = -1)
  x <- raster::intersect(x, oms:::boundary(coords))
  options(op)
  
  tab <- spbabel::sptable(x)
  xy <- as.matrix(coords)
  kd <- nabor::knn(xy, raster::as.matrix(tab[, c("x", "y")]), k = 1, eps = 0)
  index <- expand.grid(x = seq(ncol(coords)), y = rev(seq(nrow(coords))))[kd$nn.idx, ]
  tab$x <- index$x
  tab$y <- index$y
  spbabel::spFromTable(tab, crs = projection(x))
}

## this is from rastermesh
boundary <- function(cds) {
  left <- cellFromCol(cds, 1)
  bottom <- cellFromRow(cds, nrow(cds))
  right <- rev(cellFromCol(cds, ncol(cds)))
  top <- rev(cellFromRow(cds, 1))
  ## need XYFromCell method
  SpatialPolygons(list(Polygons(list(Polygon(raster::as.matrix(cds)[unique(c(left, bottom, right, top)), ])), "1")))
}

#' @rdname romsmap
#' @export
romsmap.SpatialPointsDataFrame <- function(x, coords, ...) {
  ## first get the intersection
  op <- options(warn = -1)
  x <- raster::intersect(x, oms:::boundary(coords))
  options(op)
  
  tab <- sptable::sptable(x)
  xy <- as.matrix(coords)
  kd <- nabor::knn(xy, raster::as.matrix(tab[, c("x", "y")]), k = 1, eps = 0)
  index <- expand.grid(x = seq(ncol(coords)), y = rev(seq(nrow(coords))))[kd$nn.idx, ]
  tab$x <- index$x
  tab$y <- index$y
  sptable::spFromTable(tab, crs = projection(x))
}

#' Extract coordinate arrays from ROMS. 
#' 
#' Returns a RasterStack of the given variable names. 
#'
#' @param x ROMS file name
#' @param spatial names of coordinate variables (e.g. lon_u, lat_u) 
#'
#' @return \code{\link[raster]{RasterStack}}
#' @export 
#'
#' @examples
#' \dontrun{
#'   coord <- romscoord("roms.nc")
#' }
romscoords <- function(x, spatial = c("lon_u", "lat_u")) {
  l <- vector("list", length(spatial))
  for (i in seq_along(l)) l[[i]] <- raster(x, varname = spatial[i])
  stack(l)
}

#' Extract a data lyaer from ROMS by name and slice. 
#'
#' @param x ROMS file name
#' @param varname name of ROMS variable 
#' @param slice index in w and t (depth and time), defaults to first encountered
#'
#' @return \code{\link[raster]{RasterLayer}}
#' @export
#'
#' @examples
romsdata <-function(x, varname, slice = c(1, 1)) {
  brick(x, level = slice[1L], varname = varname)[[slice[2L]]]
}


