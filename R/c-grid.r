ucoords <- function(x, ...) {
  s <- stack(raster(x, varname = "lon_u"), raster(x, varname = "lat_u"))
  setExtent(s, extent(0, ncol(s), 0, nrow(s)))
}
vcoords <- function(x, ...) {
  s <- stack(raster(x, varname = "lon_v"), raster(x, varname = "lat_v"))
  setExtent(s, extent(0, ncol(s), 0, nrow(s)))
}
psicoords <- function(x, ...) {
  # Draw vectors at PSI-points (cell's corners):
  #Upsi(i,j) = 0.5 * ( U(i,j-1) + U(i,j) )
  #Vpsi(i,j) = 0.5 * ( V(i-1,j) + V(i,j) )
  uc <- ucoords(x)
  vc <- vcoords(x)
  ## probably should drop to raw matrix to do this stuff . . .
  gex <- extent(uc) + c(0, 0, 0, -1)
  Upsi <- 0.5 * (crop(uc, gex) + setExtent(crop(uc, extent(uc) + c(0, 0, 1, 0)), gex))
  gex <- extent(vc) + c(0, -1, 0, 0)
  Vpsi <- 0.5 * (crop(vc, gex) + setExtent(crop(vc, extent(vc) + c(1, 0, 0, 0)), gex))
  
}

rhocoords <- function(x, ...) {
  # Draw vectors at interior RHO-points (cell's center):
  #Urho(i,j) = 0.5 * ( U(i,j) + U(i+1,j) )
  #Vrho(i,j) = 0.5 * ( V(i,j) + V(i,j+1) )
}

