#' @importFrom terra rast ext<-
raster_ispace <- function(x, transpose = TRUE) {
  x <- t(x[, ncol(x):1])
  r <- terra::rast(x)
  if (transpose) {
    terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  } else {
    terra::ext(r) <- terra::ext(0, nrow(r), 0, ncol(r))
  }
  r
}

# convert the depth ramp Cs_r, h (bottom depth), and cell number
# to a correctly oriented layer of depth values
romscoords_z <- function(x, cell) {
  h <- terra::rast(x, subds = "h")
  Cs_r <- rawdata(x, "Cs_r")

  hvals <- terra::extract(h, cell)[, 1L]
  m <- matrix(rep(hvals, each = length(Cs_r)) * rep(Cs_r, length(cell)),
              nrow = length(Cs_r))
  out <- terra::flip(terra::rast(m), direction = "vertical")
  terra::ext(out) <- terra::ext(0, ncol(out), 0, nrow(out))
  out
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
roms_xy <- function(x, varname = "", slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(1L, 1L, slice)
  count <- c(-1L, -1L, 1L, 1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}
#' @name romsdata
#' @export
roms_xz <- function(x, varname = "", slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(1L, slice[1L], 1L, slice[2L])
  count <- c(-1L, 1L, -1L, 1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}
#' @name romsdata
#' @export
roms_xt <- function(x, varname = "", slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(1L, slice[1L], slice[2L], 1L)
  count <- c(-1L, 1L, 1L, -1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}


#' @name romsdata
#' @export
roms_yz <- function(x, varname = "", slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(slice[1L], 1L, 1L, slice[2L])
  count <- c(1L, -1L, -1L, 1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}
#' @name romsdata
#' @export
roms_yt <- function(x, varname = "", slice = c(1L, 1L), transpose = TRUE, ...) {
  start <- c(slice[1L], 1L, slice[2L], 1L)
  count <- c(1L, -1L, 1L, -1L)
  raster_ispace(ncgetslice(x, varname, start = start, count = count))
}

#' @name romsdata
#' @export
roms_zt <- function(x, varname = "", slice = c(1L, 1L), transpose = TRUE, ...) {
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
#' @param verbose be chatty
#' @param lvar passed to `terra::rast` to specify 3rd or 4th dimension
#' @return SpatRaster
#' @export
#'
romsdata <- function(x, varname = "", slice = c(1L, 1L), transpose = TRUE, ...) {
  romsdata3d(x, varname = varname, slice = slice[2L], transpose = transpose)[[slice[1L]]]
}
#' @name romsdata
#' @export romsdata2d
romsdata2d <- romsdata

#' for romsdata3d slice must be length 1, intended to get all depths
#' @name romsdata
#' @export
romsdata3d <- function(x, varname = "", slice = 1L, transpose = TRUE, verbose = TRUE, ..., lvar = 4L) {
  stopifnot(length(slice) == 1L)
  if (is.null(x)) stop("x must be a valid NetCDF source name")

  ## terra reads all bands/layers by default from a subdataset
  x0 <- try(terra::rast(x, subds = varname), silent = TRUE)

  if (inherits(x0, "try-error")) {
    message(sprintf("cannot read in this form, need varname = ' a 4D variable in this source:\n%s", x))
    nc <- ncdf4::nc_open(x)
    if (verbose) {
      message("printing summary of source ...")
      print(nc)
    }
    ncdf4::nc_close(nc)
    stop("%s is not multi-dimensional/interpretable as a SpatRaster, try extracting in raw form with rawdata()")
  }
  if (transpose) {
    terra::ext(x0) <- terra::ext(0, ncol(x0), 0, nrow(x0))
  } else {
    terra::ext(x0) <- terra::ext(0, nrow(x0), 0, ncol(x0))
  }
  x0
}
