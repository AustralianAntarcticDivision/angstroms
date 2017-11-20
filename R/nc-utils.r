
#' Read the variable as is
#'
#' @param x netcdf file path
#' @param varname variable name
#' @param ... dots (ignored)
#' @export
rawdata <- function(x, varname, ...) UseMethod("rawdata")
#' @name rawdata
#' @export
rawdata.character <- function(x, varname, ...) {
  return(ncdf4::ncvar_get(ncdf4::nc_open(x), varname))
}
#' @name rawdata
#' @export
rawdata.NetCDF <- function(x, varname, ...) {
  rawdata(x$file$filename[1L], varname = varname, ...)
}


#' @importFrom ncdf4 nc_open nc_close ncvar_get
ncget <- function(x, varname) {
  nc <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, varname)
}

ncgetslice <- function(x, varname, start = c(1L, 1L, 1L, 1L), count = c(-1L, -1L, -1L, -1L)) {
  con <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(con))
  ncdf4::ncvar_get(con, varname, start = start, count = count)
}

# @importFrom raster getValuesBlock raster setExtent extent nlayers
# rastergetslice <- function(x, slice) {
#   ## expect slice to be c(xindex, NA, NA) or c(NA, yindex, NA)
#   ## all longitudes
#   if (is.na(slice[1]))  x1 <-  setExtent(raster(getValuesBlock(x, row = slice[2], nrows = 1L)), extent(0, ncol(x), 0, nlayers(x)))
#   ## all latitudes
#   if (is.na(slice[2]))  x1 <-  setExtent(raster(getValuesBlock(x, col = slice[1], ncols = 1L, nrows = nrow(x))), extent(0, nrow(x), 0, nlayers(x)))
#   x1
# }




