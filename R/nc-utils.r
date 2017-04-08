


#' Extract a data layer from ROMS by name and slice. 
#' 
#' Maybe this replaced by rastergetslice??
#' Returns a single slice 2D layer
#'
#' @param x ROMS file name
#' @param varname name of ROMS variable 
#' @param slice index in w and t (depth and time), defaults to first encountered
#' @param transpose the extents (ROMS is FALSE, Access is TRUE)
#' @param ... unused
#' @param ncdf default to \code{TRUE}, set to \code{FALSE} to allow raster format detection brick
#' @importFrom raster brick 
#' @return RasterLayer
#' @export
#'
romsdata <- function (x, varname, slice = c(1, 1), ncdf = TRUE, transpose = FALSE, ...) 
{
   stopifnot(!missing(varname))
  if (is.null(x)) stop("x must be a valid file name")
   x0 <- try(brick(x, level = slice[1L], varname = varname, ncdf = ncdf, ...), silent = TRUE)
  if (inherits(x0, "try-error")) {
     ## 
    stop(sprintf("%s is not multi-dimensional/interpretable as a RasterLayer, try extracting in raw form with rawdata()", varname))
    
    }
  x <- x0[[slice[2L]]]
   if (transpose) {
    e <- extent(0, ncol(x), 0, nrow(x)) 
   } else {
   e <- extent(0, nrow(x), 0, ncol(x))
   }
  setExtent(x, e)
}

#' Read the variable as is
#' 
#' @param x netcdf file path
#' @param varname variable name
#'
#' @export
rawdata <- function(x, varname) {
  return(ncdf4::ncvar_get(ncdf4::nc_open(x), varname))
 }


#' @importFrom ncdf4 nc_open nc_close ncvar_get 
ncget <- function(x, varname) {
  nc <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, varname)
}

ncgetslice <- function(x, varname, start, count) {
  con <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(con))
  ncdf4::ncvar_get(con, varname, start = start, count = count)
}

#' @importFrom raster getValuesBlock raster setExtent extent nlayers
rastergetslice <- function(x, slice) {
  ## expect slice to be c(xindex, NA, NA) or c(NA, yindex, NA)
  ## all longitudes
  if (is.na(slice[1]))  x1 <-  setExtent(raster(getValuesBlock(x, row = slice[2], nrows = 1)), extent(0, ncol(x), 0, nlayers(x)))
  ## all latitudes
  if (is.na(slice[2]))  x1 <-  setExtent(raster(getValuesBlock(x, col = slice[1], ncols = 1, nrows = nrow(x))), extent(0, nrow(x), 0, nlayers(x)))
  x1
}




