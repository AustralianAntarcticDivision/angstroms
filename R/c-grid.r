
## https://www.myroms.org/wiki/easygrid
# The rho-grid represents the grid centers (red dots)
# The u-grid represents the grid East-West sides (blue triangles)
# The v-grid represents the grid North-South sides (green triangles)
# The psi-grid represents the grid corners (purple crosses)

plot_cgrid <- function(x, ex = extent(0, 15, 0, 20)) {
  uc <- crop(ucoords(x), ex)
  vc <- crop(vcoords(x), ex)
  rc <- crop(rhocoords(x), ex)
  psi <- crop(psicoords(x), ex)
  plot(as.matrix(rc), col = "firebrick", pch = 19, cex = 0.4)
  for (i in seq(ncol(uc))) lines(extract(uc, cellFromCol(uc, i)))
  for (j in seq(nrow(vc))) lines(extract(vc, cellFromRow(vc, j)))
  points(as.matrix(uc), pch = 17, col = "blue", cex = 0.6)
  points(as.matrix(vc), pch = 17, col = "green3", cex = 0.6)
  
  
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

