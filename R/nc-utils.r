raster_ispace <- function(x, transpose = TRUE) {
  x <- t(x[,ncol(x):1])
  if (transpose) {
    e <- extent(0, ncol(x), 0, nrow(x)) 
  } else {
    e <- extent(0, nrow(x), 0, ncol(x))
  }
  setExtent(raster(x), e)
}

#' @examples 
#' #x <- raadtools:::cpolarfiles()$fullname[1]
#' #plot(roms_xy(x, "u"))
#' #plot(roms_xz(x, "u", slice = c(392L,1L)), asp = NA)
#' #plot(roms_xt(x, "u", slice = c(392L,1L)), asp = NA)
#' 
#' #plot(roms_yz(x, "u"))
#' #plot(roms_yt(x, "u", slice = c(1L,1L)), asp = NA)
#' #plot(roms_zt(x, "u", slice = c(1L, 392L)), asp = NA)
#' @name romsdata
#' @export
roms_xy <- function(x, varname, slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(1L, 1L, slice)
  count <- c(-1L, -1L, 1L, 1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}
#' @name romsdata
#' @export
roms_xz <- function(x, varname, slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(1L, slice[1L], 1L, slice[2L])
  count <- c(-1L, 1L, -1L, 1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}
#' @name romsdata
#' @export
roms_xt <- function(x, varname, slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(1L, slice[1L], slice[2L], 1L)
  count <- c(-1L, 1L, 1L, -1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}


#' @name romsdata
#' @export
roms_yz <- function(x, varname, slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(slice[1L], 1L, 1L, slice[2L])
  count <- c(1L, -1L, -1L,  1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}
#' @name romsdata
#' @export
roms_yt <- function(x, varname, slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(slice[1L], 1L,  slice[2L], 1L)
  count <- c(1L, -1L,  1L, -1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}

#' @name romsdata
#' @export
roms_zt <- function(x, varname, slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(slice, 1L, 1L)
  count <- c(1L, 1L, -1L, -1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}


#'  ROMS single slice 2D layer
#'  
#'  Extract a data layer from ROMS by name and 4-D slice. 
#' 
#' `romsdata` always works in the first two dimensions (x-y), the more specialist functions will
#' work in the space indicated by their name `roms_xy`, `roms_xt` and so on. 
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
romsdata <- function (x, varname, slice = c(1L, 1L), ncdf = TRUE, transpose = TRUE, ...) 
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

ncgetslice <- function(x, varname, start = c(1L, 1L, 1L, 1L), count = c(-1L, -1L, -1L, -1L)) {
  con <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(con))
  ncdf4::ncvar_get(con, varname, start = start, count = count)
}

#' @importFrom raster getValuesBlock raster setExtent extent nlayers
rastergetslice <- function(x, slice) {
  ## expect slice to be c(xindex, NA, NA) or c(NA, yindex, NA)
  ## all longitudes
  if (is.na(slice[1]))  x1 <-  setExtent(raster(getValuesBlock(x, row = slice[2], nrows = 1L)), extent(0, ncol(x), 0, nlayers(x)))
  ## all latitudes
  if (is.na(slice[2]))  x1 <-  setExtent(raster(getValuesBlock(x, col = slice[1], ncols = 1L, nrows = nrow(x))), extent(0, nrow(x), 0, nlayers(x)))
  x1
}




