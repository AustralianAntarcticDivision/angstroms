
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
#' @param varname in desperate cases, specify the variable that these coordinate variables belong to
#' @param flip_y Y coordinates are assumed to be in top-down order, set to FALSE to assume down-up
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
romscoords <- function(x, spatial = c("lon_u", "lat_u"), ncdf = TRUE,  transpose = TRUE, ..., varname = "", flip_y = TRUE) {
  l <- vector("list", length(spatial))
  for (i in seq_along(l)) l[[i]] <- try(raster(x, varname = spatial[i], ncdf = TRUE, ..., stopIfNotEqualSpaced = FALSE), 
                                        silent = TRUE)
  if (inherits(l[[1]], "try-error")) {
    ## assume it's rectilinear
    X <- c(rawdata(x, varname = spatial[1]))
    Y <-  c(rawdata(x, varname = spatial[2]))
    if (flip_y) Y <- rev(Y)
    xy <- expand.grid(X, Y)
  
    template <- suppressWarnings(raster(x, ncdf = TRUE, stopIfNotEqualSpaced = FALSE, varname = varname))  ## hope for the best
    
    #template <- setExtent(template, extent(0, ncol(template), 0, nrow(template)))
    l[[1]] <- setValues(template, xy[[1]])
    l[[2]] <- setValues(template, xy[[2]])
    names(l[[1]]) <- spatial[1]
    names(l[[2]]) <- spatial[2]
  }
  if (transpose) {
    l <- lapply(l, function(x) setExtent(x, extent(0, ncol(x), 0, nrow(x))))
  } else {
    l <- lapply(l, function(x) setExtent(x, extent(0, nrow(x), 0, ncol(x))))
  }
  raster::readAll(raster::stack(l))
}


