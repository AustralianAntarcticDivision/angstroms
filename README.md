## oms - R for ROMS



```r
## filename of ROMS data
fp <- "slice3101.nc"
```

```r
library(oms)
## RasterStack of the lon_u/lat_u coords
coords <- romscoords(fp)

## a polygon data set
library(rworldxtra)
data(countriesHigh)

## extract a single z/time slice of ROMS data, on native grid 
## 4th level (depth), 3rd time step
slc <- c(4, 3)

# u and v components of velocity
u <- romsdata(fp, "u", slice = slc)
v <- romsdata(fp, "v", slice = slc)
   
# temperature
temp <- romsdata(fp, "temp", slice = slc)   

## translate a SpatialPolygons data set to the ROMS mesh
## (also clips to the extents)
map <- romsmap(countriesHigh, coords)

## plot temperature + current mag contours
plot(temp, col = rev(heat.colors(100)))
contour(sqrt(u ^ 2 + v ^ 2), add = TRUE)
## add polygons
plot(map, add = TRUE, col = "grey")
```
