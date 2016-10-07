romsfile

h <- romsdata(romsfile, "h")
Cs_r <- ncvar_get(nc_open(romsfile), "Cs_r")
c1 <- 172
slice <- ncvar_get(nc_open(romsfile), "temp", 
                   start = c(c1, 1, 1, 1), count = c(1, 392, 31, 1))

u <- ncvar_get(nc_open(romsfile), "u", 
                   start = c(c1, 1, 1, 1), count = c(1, 391, 31, 1))
v <- ncvar_get(nc_open(romsfile), "v", 
                    start = c(c1, 1, 1, 1), count = c(1, 391, 31, 1))

slice <- sqrt(u * u + v * v)

rslice <- raster(rbind(slice[391, ], slice[391:1, 31:1]))

#ys <- getValuesBlock(coords[[2]], col = c1, ncols = 1, nrows = 392)
#hs <- getValuesBlock(h, col = c1, ncols = 1, nrows = 392)
ln <- spLines(cbind(c1, 0:392))
hs <- extract(h, ln)[[1]]
ys <- extract(coords[[2]], ln)[[1]]

plot(raster(outer(hs, Cs_r)))
ln <- spLines(cbind(c1, 0:392))
#plot(raster(slice[nrow(slice):1, ncol(slice):1]), col = palr::sstPal(100))
h_s <- t(raster(outer(hs, Cs_r)[,31:1 ]))
y_s <- setValues(raster(h_s), rep(ys, 31))
bd <- tabularaster::boundary(brick(y_s, h_s))

bd_tab <- raster::geom(bd)
library(dplyr)
tab <- data_frame(x = bd_tab[, "x"], y = bd_tab[, "y"])
hy <- cbind(ys, values(h_s))
uxy <- data_frame(x = hy[,1], y = hy[,2], id = seq(nrow(hy)), z = values(rslice)) %>% distinct(x, y)
bndID <- tab %>% inner_join(uxy)
library(RTriangle)
ps <- pslg(as.matrix(uxy[, c("x", "y")]), 
           S = head(matrix(bndID$id, ncol = 2, nrow = nrow(bndID) + 1), -1), 
           PA = matrix(uxy$z))


tri <- triangulate(ps)
save(tri, uxy, file = "triroms.rdata")
g <- gris:::as.gris.triangulation(tri)
g$v$z <- uxy$z
save(g, file = "groms.rdata")


library(rgl)
library(gris)
load("triroms.rdata")

o <- tetrahedron3d()
o$vb <- t(cbind(45, tri$P, 1))
o$it <- t(tri$T)
uxy$z2 <- zoo::na.locf(uxy$z, na.rm = F)
scl <- function(x) (x  - min(x))/diff(range(x))
shade3d(o, col = palr::sstPal(100)[scl(uxy$z2[o$it]) * 99 + 1], specular = "black"); aspect3d(1, .2, 0.1)






