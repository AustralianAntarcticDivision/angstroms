#' @name plot_cgrid
#' @param x ROMS file
#' @param ex extent in index space
#' @param include grid elements to plot (all by default, uv, rho, psi)
#' @param cell draw the cells (defaults to TRUE)
#' @param ... arguments passed to plot
#' @export
plot_cgrid <- function(x, ex = terra::ext(0, 15, 0, 20), 
                       include = c("u", "v", "rho", "psi"), cell = TRUE, ...) {
  
  rc <- terra::crop(rhocoords(x), ex, snap = "out")
  uc <- terra::crop(ucoords(x), ex, snap = "out")
  vc <- terra::crop(vcoords(x), ex, snap = "out")
  psi <- terra::crop(psicoords(x), ex, snap = "out")
  
  xy_rc <- cbind(terra::values(rc[[1]]), terra::values(rc[[2]]))
  plot(xy_rc, type = "n", ...) 
  if ("rho" %in% include) {
    points(xy_rc, col = "firebrick", pch = 19, cex = 0.4)
  }
  
  if ("u" %in% include) {
    xy_uc <- cbind(terra::values(uc[[1]]), terra::values(uc[[2]]))
    points(xy_uc, pch = 17, col = "blue", cex = 0.6)
  }
  if ("v" %in% include) {
    xy_vc <- cbind(terra::values(vc[[1]]), terra::values(vc[[2]]))
    points(xy_vc, pch = 17, col = "green3", cex = 0.6)
  }
  if (cell) {
    for (i in seq(ncol(uc))) {
      cells <- terra::cellFromCol(uc, i)
      vals <- cbind(terra::values(uc[[1]])[cells], terra::values(uc[[2]])[cells])
      lines(vals)
    }
    for (j in seq(nrow(vc))) {
      cells <- terra::cellFromRow(vc, j)
      vals <- cbind(terra::values(vc[[1]])[cells], terra::values(vc[[2]])[cells])
      lines(vals)
    }
  }
  invisible(NULL) 
} 
ucoords <- function(x, ...) {
  s <- c(terra::rast(x, subds = "lon_u"), terra::rast(x, subds = "lat_u"))
  terra::ext(s) <- terra::ext(0, ncol(s), 0, nrow(s))
  s
}
vcoords <- function(x, ...) {
  s <- c(terra::rast(x, subds = "lon_v"), terra::rast(x, subds = "lat_v"))
  terra::ext(s) <- terra::ext(0, ncol(s), 0, nrow(s))
  s
}
rhocoords <- function(x, ...) {
  s <- c(terra::rast(x, subds = "lon_rho"), terra::rast(x, subds = "lat_rho"))
  terra::ext(s) <- terra::ext(0, ncol(s), 0, nrow(s))
  s
}
psicoords <- function(x, ...) {
  vc <- vcoords(x)
  rc <- rhocoords(x)
  vgex <- terra::ext(vc) + c(0, -1, 0, 0)
  vc_crop1 <- terra::crop(vc, vgex)
  vc_crop2 <- terra::crop(vc, terra::ext(vc) + c(1, 0, 0, 0))
  terra::ext(vc_crop2) <- vgex
  0.5 * (vc_crop1 + vc_crop2)
}
