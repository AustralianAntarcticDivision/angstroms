#' Rotate 'ROMS' vectors
#'  
#' 
#' `romsrotate` performs the rotation on a 2-layer SpatRaster
#' `romsangle` reads the (likely) angle variable
#'
#' @param uv SpatRaster of two layers with 'u' ("east-west") and 'v' ("north-south") vector components
#' @param angle SpatRaster layer of the 'angle' variable from 'ROMS'
#' @param x  filename
#' @param varname variable name from file 'x'
#' @param ... ignored currently
#' @references  [ROMS website](https://www.myroms.org/forum/viewtopic.php?f=3&t=295)
#' @name romsrotate
#' @aliases romsangle
#' @export
romsrotate <- function(uv, angle, ...) {
  u <- uv[[1]] * cos(angle) - uv[[2]] * sin(angle)
  v <- uv[[1]] * sin(angle) + uv[[2]] * cos(angle)
  c(u, v)
}
#' @name romsrotate
#' @export
romsangle <- function(x, varname = "angle", ...) {
  romsdata2d(x, varname = varname, ...)
}
