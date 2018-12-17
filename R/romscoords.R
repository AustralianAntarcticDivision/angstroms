
#' Extract coordinate arrays from ROMS. 
#' 
#' Returns a RasterStack of the given variable names. 
#' 
#' The two layers from the model output are used to define the real-world space. This is used to create a boundary `romsboundary`, to map real-world
#' objects into  grid space `romscoords` and to generate graticules for mapping into the grid space with `graphics::contour`. 
#' @param x ROMS file name
#' @param spatial names of coordinate variables (e.g. lon_u, lat_u) 
#' @param ncdf default to NetCDF no matter what file name
#' @param transpose the extents (ROMS is FALSE, Access is TRUE)
#' @param ... unused
#'
#' @return RasterStack with two layers of the 2D-variables
#' @export 
#'
#' @examples
#' \dontrun{
#'   coord <- romscoord("roms.nc")
#' }
#' ## with in-built fake data
#' plot(ice_fake, asp = 0.5)
#' contour(ice_coords[[1]], add = TRUE, levels = seq(-165, 165, by = 15))
#' contour(ice_coords[[2]], add = TRUE)
#' 
#' @importFrom raster brick values
#' @importFrom raster stack
#' 
romscoords <- function(x, spatial = c("lon_u", "lat_u"), ncdf = TRUE,  transpose = FALSE, ... ) {
  l <- vector("list", length(spatial))
  for (i in seq_along(l)) l[[i]] <- raster(x, varname = spatial[i], ncdf = TRUE, ..., stopIfNotEqualSpaced = FALSE)
  if (transpose) {
    l <- lapply(l, function(x) setExtent(x, extent(0, ncol(x), 0, nrow(x))))
  } else {
    l <- lapply(l, function(x) setExtent(x, extent(0, nrow(x), 0, ncol(x))))
  }
  raster::readAll(raster::stack(l))
}


