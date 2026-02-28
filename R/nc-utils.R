
#' Read the variable as-is
#' 
#' Read a raw variable from a NetCDF source with no transformation or 
#' raster interpretation.
#'
#' @param x netcdf file path
#' @param varname variable name
#' @param ... dots (ignored)
#' @param native maintain values in internal storage mode (default is `FALSE`)
#' @export
rawdata <- function(x, varname, ..., native = FALSE) UseMethod("rawdata")
#' @name rawdata
#' @export
rawdata.character <- function(x, varname, ..., native = FALSE) {
  nc <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, varname, raw_datavals = native)
}
#' @name rawdata
#' @export
rawdata.NetCDF <- function(x, varname, ..., native = FALSE) {
  rawdata(x$file$filename[1L], varname = varname, ..., native = native)
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
