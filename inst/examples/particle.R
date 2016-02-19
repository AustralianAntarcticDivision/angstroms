## TODO:
##   consider use of polar-projected coordinates for kdtree search 
##   need an out of bounds test, similar to the hit-the-bottom test
## PREP 
## a single file from BGF
## f <- "/ds/projects/mertz/mdl/mer015_1/mer_his_0574.nc"
## on katabatic
## single time slice from BGF (to be replaced by summer climatology)
## module load netcdf/3.6.2-gnu
## module load nco/3.9.5
## ncks -v u,v,w,Cs_w,h -d ocean_time,1 /ds/projects/mertz/mdl/mer015_1/mer_his_0574.nc ~/mer_his_2017Mar03.nc
## ftp krill (etc)


## ROMS slice
sf <- "http://staff.acecrc.org.au/~mdsumner/mer_his_2017Mar03.nc"
if (!file.exists(basename(sf))) download.file(sf, basename(sf), mode = "wb")
f <- basename(sf)

## nabor should install with install.packages("nabor"), ncdf4 already on VM
library(ncdf4)
library(nabor)

nc <- nc_open(f)

## obtain relevant coordinate arrays (these are [lon lat])
## u,v,rho distinction is ignored below 
lon_u <- ncvar_get(nc, "lon_u")
lat_u <- ncvar_get(nc, "lat_u")

## h is bathymetry, we need Cs function of w to 
## determine every cell's depth
h0 <- ncvar_get(nc, "h")
Cs_w <- ncvar_get(nc, "Cs_w")

## we read u,v,w from our dummy file which has only one time slice
## (this might be in the loop later)
i_u <- ncvar_get(nc, "u", start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
i_v <- ncvar_get(nc, "v", start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
i_w <- ncvar_get(nc, "w", start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))

## be nice to kittens
nc_close(nc)

## top to bottom
Cs_wr <- rev(Cs_w)

## populate the depths in one big array [lon lat w]
## remember h is one row extra in the u dim
h <- h0[-nrow(h0),]
hh <- array(rep(as.vector(h), length(Cs_wr)) * rep(Cs_wr, each = prod(dim(h))), c(dim(h), length(Cs_wr)))


## we want one big table of all coordinates for kdtree search
## lon and lat get recycled
allxyz <- cbind(as.vector(lon_u), as.vector(lat_u), as.vector(hh))
## needs to be a separate lookup for detecting when we hit the bottom
xytop <- seq_len(prod(dim(lon_u)))

## parameters for run 
ntime <- 1200
time_step <- 30 * 60 
w_sink <- -0.001
## only keep a sample of the full trajectory
thin <- ntime %/% 10


## these are the starting points for particles
## (replace with chla-locations)
pts <- allxyz[xytop, ][sample(length(xytop), length(xytop) %/% 20), ]

## these are the particles! 
ptrack <- array(0, c(nrow(pts), 3, ntime/thin))

## this is our search tree engine, built once upfront
kdtree <- WKNND(allxyz)
## we store a second tree for the halting condition
kdxy <-  WKNND(allxyz[xytop, 1:2 ]) ## knn(allxyz[xytop, ], xyzt[,1:2,itime - 1], k = 1)

## halt
stopped <- rep(FALSE, nrow(pts))

## copies of the starting points for updating in the loop
plast <- pts
pnow <- plast


plotit <- TRUE
if (plotit) plot(pts, pch = ".")
system.time({
for (itime in seq_len(ntime)) {
  ## index 1st nearest neighbour of trace points to grid points
  dmap <- kdtree$query(plast, k = 1, eps = 0)
  ## extract component values from the vars
  thisu <- i_u[dmap$nn.idx]
  thisv <- i_v[dmap$nn.idx]
  thisw <- i_w[dmap$nn.idx]
  
  ## update this time step longitude, latitude, depth
  pnow[,1] <- plast[,1] + (thisu * time_step) / (1.852 * 60 * 1000 * cos(pnow[,2] * pi/180))
  pnow[,2] <- plast[,2] + (thisv * time_step) / (1.852 * 60 * 1000)
  pnow[,3] <- pmin(0, plast[,3])  + ((thisw + w_sink)* time_step )
  
  ## hit the bottom 
  stopped <- stopped | pnow[,3] <= -h[kdxy$query(pnow, k = 1, eps = 0)$nn.idx]
  stopped[is.na(stopped)] <- TRUE
  pnow[stopped,] <- plast[stopped,]
  if (itime %% thin == 0) {
    ptrack[,,itime/thin] <- pnow
    if ((itime/thin) > 1 & plotit) {
    arrows(ptrack[,1,(itime/thin) - 1], ptrack[,2,(itime/thin) - 1], 
           pnow[,1], pnow[,2], length = 0)
    points(plast[is.na(pnow[,1]), ], pch = 16, cex = 0.5, col = "firebrick")
    }
    
  }
  ##if (any(is.na(plast)) | any(is.na(pnow))) stop()
  plast <- pnow
  print(itime)
  if (all(stopped)) break; 
}
})



## there's no way to start rgl in RStudio Server, so save out the ptrack object to visualize at home
## save(ptrack, file = "ptrack.Rdata")
## library(rgl)
## plot3d(ptrack[,,dim(ptrack)[3]])
