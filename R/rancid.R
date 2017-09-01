#' Read an arbitrary 2D or 3D slice from NetCDF as a RasterBrick
#' 
#' @param x ROMS file name
#' @param varname variable name
#' @param slice index, specified with NA for the index to read all steps
#'
#' @noRd
# ncraster <- function(x, varname, slice) {
#   nc <- rancid::NetCDF(x)
#   vd <- ## how is order controlled here?
#     rancid::vars(nc) %>% filter(name == varname) %>% 
#     inner_join(nc$vardim, "id") %>% transmute(vid = id, id = dimids) %>% 
#     inner_join(dims(nc), "id")
#   ## if slice is NA, we get all
#   start <- ifelse(is.na(slice), 1, slice)
#   count <- ifelse(is.na(slice), vd$len, 1)
#   # print(start)
#   #  print(count)
#   a <- ncgetslice(x, varname, start, count)
#   if (length(dim(a)) == 2) {
#     a <- a[, ncol(a):1 ]
#     a <- setExtent(raster(t(a)), extent(0, nrow(a), 0, ncol(a)))
#   } else {
#     a <- a[,ncol(a):1,]
#     a <- setExtent(brick(a,  transpose = TRUE)  , extent(0, nrow(a), 0, ncol(a)))
#   }
#   a
# }


## these should all be in rancid?
## if so must be exported . . .


#' NetCDF variable dimension
#'
#' This belongs in rancid . . .
#'
#' @param varname variable name
#' @param x file
#'
#' @noRd
#'
#' @importFrom dplyr transmute
# ncdim <- function(x, varname) {
#   roms <- NetCDF(x)
#   # ## still exploring neatest way to do this . . .
#   vdim <- vars(roms) %>% 
#     filter(name == varname) %>% 
#     inner_join(roms$vardim, "id") %>% 
#     dplyr::transmute(id = dimids) %>% 
#     inner_join(dims(roms), "id") 
#   vdim$len
# }
