## these should all be in rancid?
## if so must be exported . . .


#' NetCDF variable dimension
#'
#' This belongs in rancid . . .
#'
#' @param varname variable name
#' @param x file
#'
#' @export
#'
#' @importFrom dplyr transmute
ncdim <- function(x, varname) {
  roms <- NetCDF(x)
  # ## still exploring neatest way to do this . . .
  vdim <- vars(roms) %>% 
    filter(name == varname) %>% 
    inner_join(roms$vardim, "id") %>% 
    dplyr::transmute(id = dimids) %>% 
    inner_join(dims(roms), "id") 
  vdim$len
}



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
#'
#' @return RasterLayer
#' @export
#'
romsdata <- function (x, varname, slice = c(1, 1), ncdf = TRUE, transpose = FALSE, ...) 
{
   stopifnot(!missing(varname))
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

rastergetslice <- function(x, slice) {
  ## expect slice to be c(xindex, NA, NA) or c(NA, yindex, NA)
  ## all longitudes
  if (is.na(slice[1]))  x1 <-  setExtent(raster(getValuesBlock(x, row = slice[2], nrows = 1)), extent(0, ncol(x), 0, nlayers(x)))
  ## all latitudes
  if (is.na(slice[2]))  x1 <-  setExtent(raster(getValuesBlock(x, col = slice[1], ncols = 1, nrows = nrow(x))), extent(0, nrow(x), 0, nlayers(x)))
  x1
}



#' Read an arbitrary 2D or 3D slice from NetCDF as a RasterBrick
#' 
#' @param x ROMS file name
#' @param varname variable name
#' @param slice index, specified with NA for the index to read all steps
#'
#' @export
ncraster <- function(x, varname, slice) {
  nc <- rancid::NetCDF(x)
  vd <- ## how is order controlled here?
    rancid::vars(nc) %>% filter(name == varname) %>% 
    inner_join(nc$vardim, "id") %>% transmute(vid = id, id = dimids) %>% 
    inner_join(dims(nc), "id")
  ## if slice is NA, we get all
  start <- ifelse(is.na(slice), 1, slice)
  count <- ifelse(is.na(slice), vd$len, 1)
 # print(start)
#  print(count)
  a <- ncgetslice(x, varname, start, count)
  if (length(dim(a)) == 2) {
    a <- a[, ncol(a):1 ]
    a <- setExtent(raster(t(a)), extent(0, nrow(a), 0, ncol(a)))
  } else {
    a <- a[,ncol(a):1,]
    a <- setExtent(brick(a,  transpose = TRUE)  , extent(0, nrow(a), 0, ncol(a)))
  }
  a
}
