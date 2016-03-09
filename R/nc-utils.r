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


#' Read an arbitrary 2D or 3D slice from NetCDF as a RasterBrick
#'
#' @param x ROMS file name
#' @param varname variable name
#' @param slice index, specified with NA for the index to read all steps
#'
#' @return
#' @export
#'
#' @examples
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
