
#' Extract coordinate arrays from ROMS. 
#' 
#' Returns a SpatRaster of the given variable names. 
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
#' @return SpatRaster with two layers of the 2D-variables
#' @export 
#'
#' @examples
#' \dontrun{
#'   coord <- romscoords("roms.nc")
#' }
#' ## with in-built fake data
#' plot(ice_fake, asp = 0.5)
#' terra::contour(ice_coords[[1]], add = TRUE, levels = seq(-165, 165, by = 15))
#' terra::contour(ice_coords[[2]], add = TRUE)
#' 
romscoords <- function(x, spatial = c("lon_u", "lat_u"), ncdf = TRUE,  transpose = TRUE, ..., varname = "", flip_y = TRUE) {
  l <- vector("list", length(spatial))
  for (i in seq_along(l)) {
    l[[i]] <- try(terra::rast(x, subds = spatial[i]), silent = TRUE)
  }
  if (inherits(l[[1]], "try-error")) {
    ## assume it's rectilinear
    X <- c(rawdata(x, varname = spatial[1]))
    Y <- c(rawdata(x, varname = spatial[2]))
    if (flip_y) Y <- rev(Y)
    xy <- expand.grid(X, Y)

    template <- suppressWarnings(terra::rast(x, subds = varname))
    l[[1]] <- terra::setValues(template, xy[[1]])
    l[[2]] <- terra::setValues(template, xy[[2]])
    names(l[[1]]) <- spatial[1]
    names(l[[2]]) <- spatial[2]
  }
  out <- c(l[[1]], l[[2]])
  if (transpose) {
    terra::ext(out) <- terra::ext(0, ncol(out), 0, nrow(out))
  } else {
    terra::ext(out) <- terra::ext(0, nrow(out), 0, ncol(out))
  }
  ## force values into memory
  terra::values(out) <- terra::values(out)
  out
}
