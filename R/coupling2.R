
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
coords_points <- function(x, ...) {
  SpatialPoints(cbind(values(x[[1]]), raster::values(x[[2]])), 
                proj4string = CRS("+init=epsg:4326"))
}


#' Convenience function to transform map projection . 
#'
#' Transform `x` to whatever the projection of `to` is. 
#' @param x  object to transform
#' @param to object with a map projection
#'
#' @return `x`, transformed
#' @export
#' @importFrom raster projection
#' @importFrom sp spTransform
project_to <- function(x, to) {
  spTransform(x, CRS(projection(to)))
}

