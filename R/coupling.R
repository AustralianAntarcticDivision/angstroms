# ## worker functions for ROMS netcdf and R raster
# mkget <- function(filename) {
#   function(varname, start = NA, count = NA) {
#     on.exit(nc_close(handle))
#     handle <- nc_open(filename)
#      ncvar_get(handle, varname, start, count)
#   }
# }
# 
# applyget <- function(vars) {
#   a <- vector("list", length(vars))
#   names(a) <- vars
#   for (i in seq_along(vars)) a[[vars[i]]] <- getr(vars[i])
#   a
# }
# 
# 
# setIndexExt <- function(x) {
#   ex <- extent(1, ncol(x), 1, nrow(x))
#   setExtent(x, ex)
# }
# nctor <- function(x) {
#   if (length(dim(x)) == 2) {
#     x <- flip(raster(t(x)), "y")
#   } else {
#     if (length(dim(x)) == 3) x <- flip(brick(x, transpose = TRUE), "y")
#   } 
#   setIndexExt(x)
# }
