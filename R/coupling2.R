
#' Create SpatialPoints. 
#'
#' Convenience wrapper around `SpatialPoints` for a two layer brick with longitude and latitude coordinate arrays.  
#' @param x two layer `RasterBrick` with longitude and latitude values
#' @param ... ignored
#'
#' @return `SpatialPoints`
#' @export
#' @importFrom sp SpatialPoints CRS
#' @importFrom raster values
#' @examples
#' ## library(raadtools)
#' ##coords_points(romscoords(cpolarfiles()$fullname[1]))
#' 
#' pts <- coords_points(ice_coords)
coords_points <- function(x, ...) {
  SpatialPoints(cbind(values(x[[1]]), raster::values(x[[2]])), 
                proj4string = CRS("+init=epsg:4326"))
}



