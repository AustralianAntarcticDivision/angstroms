raster_ispace <- function(x, transpose = TRUE) {
  x <- t(x[,ncol(x):1])
  if (transpose) {
    e <- extent(0, ncol(x), 0, nrow(x)) 
  } else {
    e <- extent(0, nrow(x), 0, ncol(x))
  }
  setExtent(raster(x), e)
}

# convert the depth ramp Cs_r, h (bottom depth), and cell number
# to a correctly oriented layer of depth values
romscoords_z <- function(x, cell) {
  ## important to readAll here, else extract is very slow in the loop
  h <- raster::readAll(raster(x, varname = "h", ncdf = TRUE))
  ## Cs_r is the S-coord stretching
  Cs_r <- rawdata(x, "Cs_r")
  
  out <- flip(raster(matrix(rep(extract(h, cell), each = length(Cs_r)) *  rep(Cs_r, length(cell)), 
                            length(Cs_r))), "y")
  setExtent(out, extent(0, ncol(out), 0, nrow(out)))
}

## read the 180th (reading up) latitude
## which happens to cut the coast a few times
#zz <- roms_xz(f, "temp", slice = c(180, 1))

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
#' @param verbose be chatty
#' @param lvar passed to `raster::brick` to specify 3rd or 4th dimension
#' @importFrom raster brick 
#' @return RasterLayer
#' @export
#'
romsdata <- function (x, varname, slice = c(1L, 1L), transpose = TRUE, ...) 
{
  if (missing(varname)) {
    stop("no varname supplied")
  }
  romsdata3d(x, varname = varname, slice = slice[2L], transpose = transpose)[[slice[1L]]]
}
#' @name romsdata
#' @export
romsdata2d <- romsdata
#' for romsdata3d slice must be length 1, intended to get all depths
#' @name romsdata
#' @export
romsdata3d <- function (x, varname, slice = 1L, transpose = TRUE, verbose = TRUE,  ..., lvar = 4L) 
{
  stopifnot(length(slice) == 1L)
  if (is.null(x)) stop("x must be a valid NetCDF source name")
  ## why is ncdf = TRUE needed? (maybe if the filename is not *.nc ...)
  x0 <- try(brick(x, level = slice[1L], lvar = lvar, varname = varname, ncdf = TRUE, ..., stopIfNotEqualSpaced = FALSE), silent = TRUE)

  if (inherits(x0, "try-error")) {
    message(sprintf("cannot read in this form, need varname = ' a 4D variable in this source:\n%s", x))
    # tnc <- try(tidync::tidync(x))
    tnc <- ncdf4::nc_open(x)
    if (!inherits(tnc, "try-error") && verbose) {
      message("printing summary of source ...")
      print(tnc)
      
    }
    stop("%s is not multi-dimensional/interpretable as a RasterLayer, try extracting in raw form with rawdata()")
  } 
  if (transpose) {
    e <- extent(0, ncol(x0), 0, nrow(x0)) 
  } else {
    e <- extent(0, nrow(x0), 0, ncol(x0))
  }
  setExtent(x0, e)
}
