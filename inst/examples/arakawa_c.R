library(raster)
library(gris)
library(rgl)
library(RTriangle)

f <- "data/mer_his_1992_01.nc"
ulon <- raster(f, varname = "lon_u")
ulat <- raster(f, varname = "lat_u")
vlon <- raster(f, varname = "lon_v")
vlat <- raster(f, varname = "lat_v")
rlon <- raster(f, varname = "lon_rho")
rlat <- raster(f, varname = "lat_rho")

.crop <- function(x, ext = NULL) {
  if (!is.null(ext)) x <- crop(x, ext)
  x
}
mkget <- function(x) {
  function(varname = "lon_u", xylim = NULL) {
    .crop(raster(x, varname = varname), xylim)
  }
}
mkvget <- function(getter, xylim) {
  function(varname) {
  values(getter(varname, xylim))
  }
}
mkrvget <- function(getter) {
  function(varname, row)  {
    dat <- getter(varname)
    extract(dat, cellFromRow(dat, row))
  }
}
mkcvget <- function(getter) {
  function(varname, col) {
    dat <- getter(varname)
    extract(dat, cellFromCol(dat, col))
  }
}

ext <- extent(340, 360, 1, 20)

get <- mkget(f)
vget <- mkvget(get, ext)
rowget <- mkrvget(get)
colget <- mkcvget(get)

crho <- "firebrick"
cu <- "dodgerblue"
cv <- "green"
par(pch = 16, cex = 0.5)
plot(vget("lon_u"), vget("lat_u"), xlim = range(vget("lon_u"), vget("lon_v"), vget("lon_rho")), col = cu)
points(vget("lon_v"), vget("lat_v"), col = cv)
points(vget("lon_rho"), vget("lat_rho"), col = crho)

for (i in 1:359) {
lines(colget("lon_u", i), colget("lat_u", i), col = "red")
}

for (i in 1:199) {
  lines(rowget("lon_v", i), rowget("lat_v", i), col = "red")
}

l <- vector("list", 46)
cnt <- 0
for (i in 1:359) {
  cnt <- cnt + 1
  l[[cnt]] <- cbind(colget("lon_u", i), colget("lat_u", i))
}

for (i in 1:199) {
  cnt <- cnt + 1
  l[[cnt]] <- cbind(rowget("lon_v", i), rowget("lat_v", i))
}

library(sp)
ll <- SpatialLines(list(Lines(lapply(l, Line), "1")))
library(rgeos)
g <- gPolygonize(gNode(ll))

