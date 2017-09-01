#' Remap an object to the ROMS grid. 
#' 
#' Find the nearest-neighbour coordinates of `x` in the coordinate arrays of `coords`. 
#' 
#' The input `coords` is a assumed to be a 2-layer RasterStack or RasterBrick and
#' using `nabor::knn` the nearest matching position of the coordinates of `x` is found in the grid space of `coords`. The
#' motivating use-case is the curvilinear longitude and latitude arrays of ROMS model output. 
#' 
#' No account is made for the details of a ROMS cell, though this may be included in future. We tested only with the "lon_u" and "lat_u"
#' arrays. 
#' @param x object to transform to the grid space, e.g. a \code{\link[sp]{Spatial}} object
#' @param coords romscoords RasterStack
#' @param crop logical, if \code{TRUE} crop x to the extent of the boundary of the values in coords
#' @param lonlat logical, if \code{TRUE} check for need to back-transform to longitude/latitude and do it
#' @param ... unused
#' @note Do not use this for extraction purposes without checking the output, this is best used for exploration
#' and visualization. Re-mapping ROMS data is better done by looking up the `coords_points` within spatial objects, 
#' and transferring via the grid index. 
#' @return input object with coordinates transformed to space of the coords 
#' @export
#' @examples 
#' ant_ice_coords <- romsmap(antarctica, ice_coords)
#' plot(ice_fake, main = "sea ice in pure grid space")
#' plot(ant_ice_coords, add = TRUE)
#' 
#' 
romsmap <- function(x, ...) {
  UseMethod("romsmap")
}

#' @rdname romsmap
#' @export
#' @importFrom spbabel sptable sp
#' @importFrom nabor knn
#' @importFrom raster intersect as.matrix projection
#' @importFrom sp CRS
romsmap.SpatialPolygonsDataFrame <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  ## first get the intersection
  if (crop) {
  op <- options(warn = -1)
  x <- raster::intersect(x, romsboundary(coords))
  options(op)
  }
  ## do we need to invert projection?
  repro <- !raster::isLonLat(x)
  proj <- projection(x)
  tab <- spbabel::sptable(x)
  
  if (repro & !is.na(proj)) {
    llproj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
    xy <- proj4::ptransform(cbind(tab$x_, tab$y_, 0), src.proj = proj, dst.proj = llproj, silent = FALSE) * 180 / pi
    ## spbabel standard is attributes with underlines
    tab$x_ <- xy[,1]
    tab$y_ <- xy[,2]
    proj <- llproj
  }
  xy <- as.matrix(coords)
  kd <- nabor::knn(xy, raster::as.matrix(tab[, c("x_", "y_")]), k = 1, eps = 0)
  index <- expand.grid(x = seq(ncol(coords)), y = rev(seq(nrow(coords))))[kd$nn.idx, ]
 tab$x_ <- index$x
  tab$y_ <- index$y
  spbabel::sp(tab, attr_tab = as.data.frame(x), crs = NA_character_)
}

#' @rdname romsmap
#' @export
romsmap.SpatialLinesDataFrame <- romsmap.SpatialPolygonsDataFrame

#' @rdname romsmap
#' @export
romsmap.SpatialPointsDataFrame <- romsmap.SpatialPolygonsDataFrame

#' Boundary polygon from raster of coordinates. 
#' 
#' Create a boundary polygon by tracking around coordinates stored in a RasterStack. 
#' 
#' The first layer in the stack is treated as the X coordinate, second as Y. 
#' @param cds two-layer Raster
#'
#' @importFrom sp SpatialPolygons Polygons Polygon
#' @importFrom raster as.matrix cellFromRow cellFromCol xmin xmax ymin ymax trim setExtent setValues raster extract flip extent 
#' @export
#' @examples 
#' ice_grid_boundary <- romsboundary(ice_coords)
#' plot(antarctica)
#' ## does not make sense in this space
#' plot(ice_grid_boundary, add = TRUE, border = "grey")
#' 
#' ## ok in this one
#' #library(rgdal)
#'#   proj4string(ice_grid_boundary) <- CRS("+init=epsg:4326")
#'# pweird <- "+proj=laea +lon_0=147 +lat_0=-42 +ellps=WGS84"
#'#  laea_world <- spTransform(antarctica, pweird)
#'#  plot(extent(laea_world) + 8e6, type = "n", asp = 1)
#'#  plot(laea_world, add = TRUE)
#'#  plot(spTransform(ice_grid_boundary, pweird), add  = TRUE, border = "darkgrey")
romsboundary <- function(cds) {
  left <- cellFromCol(cds, 1)
  bottom <- cellFromRow(cds, nrow(cds))
  right <- rev(cellFromCol(cds, ncol(cds)))
  top <- rev(cellFromRow(cds, 1))
  ## need XYFromCell method
  SpatialPolygons(list(Polygons(list(Polygon(raster::as.matrix(cds)[unique(c(left, bottom, right, top)), ])), "1")))
}


## put any raster into xy-index space (0, nc, 0, nr)
set_indextent <- function(x) {
  setExtent(x, extent(0, ncol(x), 0, nrow(x)))
}


#' Extract coordinate arrays from ROMS. 
#' 
#' Returns a RasterStack of the given variable names. 
#' 
#' The two layers from the model output are used to define the real-world space. This is used to create a boundary `romsboundary`, to map real-world
#' objects into  grid space `romscoords` and to generate graticules for mapping into the grid space with `graphics::contour`. 
#' @param x ROMS file name
#' @param spatial names of coordinate variables (e.g. lon_u, lat_u) 
#' @param ncdf default to NetCDF no matter what file name
#' @param transpose the extents (ROMS is FALSE, Access is TRUE)
#' @param ... unused
#'
#' @return RasterStack with two layers of the 2D-variables
#' @export 
#'
#' @examples
#' \dontrun{
#'   coord <- romscoord("roms.nc")
#' }
#' ## with in-built fake data
#' plot(ice_fake, asp = 0.5)
#' contour(ice_coords[[1]], add = TRUE, levels = seq(-165, 165, by = 15))
#' contour(ice_coords[[2]], add = TRUE)
#' 
#' @importFrom raster brick values
#' @importFrom raster stack
#' 
romscoords <- function(x, spatial = c("lon_u", "lat_u"), ncdf = TRUE,  transpose = FALSE, ... ) {
  l <- vector("list", length(spatial))
  for (i in seq_along(l)) l[[i]] <- raster(x, varname = spatial[i], ncdf = TRUE, ...)
  if (transpose) {
   l <- lapply(l, function(x) setExtent(x, extent(0, ncol(x), 0, nrow(x))))
   } else {
  l <- lapply(l, function(x) setExtent(x, extent(0, nrow(x), 0, ncol(x))))
  }
  raster::stack(l)
}


#' Coordinates at depth
#' 
#' Extract the multi-layer 'h'eight grid with S-coordinate stretching applied
#' 
#'  \code{S} and \code{h} are the  names of the appropriate variables
#' @param x ROMS file name 
#' @param depth depth thing
#' @param S  of S-coordinate stretching curve at RHO-points
#'
#' @return RasterStack with a layer for every depth
#' @export
romshcoords <- function(x, S = "Cs_r", depth = "h"){
  h <- raster(x, varname = depth)
  Cs_r <- ncget(x, S)
  v <- values(h)
  setExtent(brick(array(rep(rev(Cs_r), each = length(v)) * v, c(ncol(h), nrow(h), length(Cs_r))), transpose = TRUE), 
            extent(0, ncol(h), 0, nrow(h)))
}




