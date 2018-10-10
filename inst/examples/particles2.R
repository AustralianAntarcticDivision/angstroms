library(raadtools)
library(angstroms)

## 3 years of monthly files
## each file is 31 (or thereabouts) days
## x, y, z, t  (z is also length 31)
files <- cpolarfiles()

afile <- files$fullname[1]
u <- romsdata3d(afile, varname = "u")
v <- romsdata3d(afile, varname = "v")
ll <- romscoords(afile)
h <- romsdepth(afile, grid_type = "u")

keep_valid <- function(xy) {
  ## xy is lon/lat
  xy[!is.na(extract(u[[1]], romsmap(xy, ll))), ]
}
library(dplyr)
particles <- sample_n(tibble::as_tibble(values(ll)) %>% 
                        setNames(c("lon", "lat")), 1000) %>% 
  keep_valid()

cols <- viridis::viridis(1000)
mxy <- function(x, y) sqrt(x[[1]]^2 + y[[1]]^2)
xyz <- cbind(coordinates(romsmap(particles, ll)), 0)

par(mar = rep(0, 4))
plot(mxy(u, v), col = cols, asp = NA)
points(xyz)

library(nabor)
## give every *valid* lon, lat, z to the nabor engine
xyz <- cbind(values(ll[[1]]), values(ll[[2]]))
xyz <- cbind(xyz[, 1], xyz[, 2], as.vector(extract(h, romsmap(xyz[,1:2], ll))))
lookup <- WKNNF(xyz)


lookup$query(cbind(150, -43, -10), k = 1, eps = 0)



extract_uvh <- function(particles, z = 0) {
  ## we need romsmap to be non-spatial
  ## https://github.com/AustralianAntarcticDivision/angstroms/issues/12
  xyindex <- coordinates(romsmap(particles, ll))
  hindex <- extract(h, xyindex)
  ## subtract z by expanding into this matrix and
  ## find the first intersect
  image(hindex)
  contour(hindex, levels = -10, add = T)
contour(hindex, levels = -100, add = T)
contour(hindex, levels = -1000, add = T)
rowM
}

