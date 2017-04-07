
#' Create SpatialPoints. 
#'
#' Convenience wrapper around `SpatialPoints` for a two layer brick with longitude and latitude coordinate arrays.  
#' @param x two layer `RasterBrick` with longitude and latitude values
#' @param ... ignored
#'
#' @return `SpatialPoints`
#' @export
#' @examples
#' coords_points(romscoords(cpolarfiles()$fullname[1]))
coords_points <- function(x, ...) {
  SpatialPoints(cbind(values(x[[1]]), values(x[[2]])), 
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
#'
project_to <- function(x, to) {
  spTransform(x, CRS(projection(to)))
}

#' Circumpolar ROMS files. 
#' 
#' 
#' @param ... ignored
#'
#' @return data frame of fullname, date
#' @export
#'
#' @examples
#' cpolarfiles()
cpolarfiles <- function(...) {
  dplyr::mutate(dplyr::filter(raadtools::allfiles(), 
      grepl("s_corney/cpolar", fullname)), 
      date = as.POSIXct(strptime(sprintf("197%s-01", substr(basename(fullname), 12, 14)), "%Y%m-%d")))
}