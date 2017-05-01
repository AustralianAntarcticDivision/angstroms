# ##gris:::quadmeshFromRaster
# 
# quadmesh <- 
# function (x, z = NULL, na.rm = FALSE) 
# {
#   
#   rc <- rhocoords(x)
#   uc <- ucoords(x)
#   vc <- vcoords(x)
#   psi <- psicoords(x)
#   
#   x <- x[[1]]
#   exy <- edgesXY(x)
#   ind <- apply(prs(seq(ncol(x) + 1)), 1, p4, nc = ncol(x) + 
#                  1)
#   ind0 <- as.vector(ind) + rep(seq(0, length = nrow(x), by = ncol(x) + 
#                                      1), each = 4 * ncol(x))
#   if (na.rm) {
#     ind1 <- matrix(ind0, nrow = 4)
#     ind0 <- ind1[, !is.na(values(x))]
#   }
#   q3d <- NULL
#   data(q3d)
#   ob <- q3d
#   if (!is.null(z)) 
#     z <- extract(z, exy, method = "bilinear")
#   else z <- 0
#   ob$vb <- t(cbind(exy, z, 1))
#   ob$ib <- matrix(ind0, nrow = 4)
#   ob
# }

#' Crop a ROMS layer
#' 
#' Crop a ROMS data layer from `romsdata` with a raster extent. 
#' 
#' The spatial crop is performed in the coordinate space of roms data. 
#' @param x ROMS xy- coordinates, see `romscoords`
#' @param ext `raster::extent` in the coordinate system of `x`
#' @param ... ignored
#'
#' @export
#' @examples
#' ## notice that extent is in long-lat, but ice_local is in the grid
#' ## space of ice_coords
#' ice_local <- croproms(ice_coords, extent(100, 120, -75, -60))
#' plot(ice_coords[[2]], col = grey(seq(0, 1, length  = 20)))
#' plot(crop(ice_fake, ice_local), add = TRUE)
croproms <- function(x, ext, ...) {
  xy <- as.matrix(x)
  incells <- which(xy[,1] >= xmin(ext) & xy[,1] <= xmax(ext) &
                     xy[,2] >= ymin(ext) & xy[,2] <= ymax(ext))
  x1 <- setValues(x[[1]], NA_real_)
  x1[incells] <- 0
  extent(trim(x1))
}
## https://www.myroms.org/wiki/easygrid
# The rho-grid represents the grid centers (red dots)
# The u-grid represents the grid East-West sides (blue triangles)
# The v-grid represents the grid North-South sides (green triangles)
# The psi-grid represents the grid corners (purple crosses)

#' @importFrom raster crop 
plot_cgrid <- function(x, ex = extent(0, 15, 0, 20), 
                       include = c("u", "v", "rho", "psi"), cell = TRUE,  ...) {
  
  rc <- crop(rhocoords(x), ex, snap = "out")
  uc <- crop(ucoords(x), ex, snap = "out")
  vc <- crop(vcoords(x), ex, snap = "out")
  psi <- crop(psicoords(x), ex, snap = "out")
  
  plot(as.matrix(rc), type = "n", ...) 
  if ("rho" %in% include) {
    points(as.matrix(rc), col = "firebrick", pch = 19, cex = 0.4)
  }
  
  if ("u" %in% include) {
    points(as.matrix(uc), pch = 17, col = "blue", cex = 0.6)
  }
  if ("v" %in% include) {
    points(as.matrix(vc), pch = 17, col = "green3", cex = 0.6)
  }
  if (cell) {
    for (i in seq(ncol(uc))) lines(extract(uc, cellFromCol(uc, i)))
    for (j in seq(nrow(vc))) lines(extract(vc, cellFromRow(vc, j)))
 }
 invisible(NULL) 
} 
ucoords <- function(x, ...) {
  s <- stack(raster(x, varname = "lon_u"), raster(x, varname = "lat_u"))
  setExtent(s, extent(0, ncol(s), 0, nrow(s)))
}
vcoords <- function(x, ...) {
  s <- stack(raster(x, varname = "lon_v"), raster(x, varname = "lat_v"))
  setExtent(s, extent(0, ncol(s), 0, nrow(s)))
}
rhocoords <- function(x, ...) {
  s <- stack(raster(x, varname = "lon_rho"), raster(x, varname = "lat_rho"))
  setExtent(s, extent(0, ncol(s), 0, nrow(s)))
}
psicoords <- function(x, ...) {
  # Draw vectors at PSI-points (cell's corners):
  #Upsi(i,j) = 0.5 * ( U(i,j-1) + U(i,j) )
  #Vpsi(i,j) = 0.5 * ( V(i-1,j) + V(i,j) )
#  uc <- ucoords(x)
  vc <- vcoords(x)
  rc <- rhocoords(x)
  ## probably should drop to raw matrix to do this stuff . . .
 # ugex <- extent(uc) + c(0, 0, 0, -1)
  vgex <- extent(vc) + c(0, -1, 0, 0)
 # crop(uc, ugex)
  # these are the same, more or less
 # Upsi <- 0.5 * (crop(uc, ugex) + setExtent(crop(uc, extent(uc) + c(0, 0, 1, 0)), ugex))
  #   
  #Vpsi <- 0.5 * (crop(vc, vgex) + setExtent(crop(vc, extent(vc) + c(1, 0, 0, 0)), vgex))
  
  0.5 * (crop(vc, vgex) + setExtent(crop(vc, extent(vc) + c(1, 0, 0, 0)), vgex))
}


# rhocoords <- function(x, ...) {
#   # Draw vectors at interior RHO-points (cell's center):
#   #Urho(i,j) = 0.5 * ( U(i,j) + U(i+1,j) )
#   #Vrho(i,j) = 0.5 * ( V(i,j) + V(i,j+1) )
# }

