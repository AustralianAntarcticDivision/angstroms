# romsmap.sf <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
#   ## first get the intersection
#   #  if (crop) {
#   #    op <- options(warn = -1)
#   #    x <- raster::intersect(x, romsboundary(coords))
#   #    options(op)
#   #  }
#   ## do we need to invert projection?
#   #  repro <- !raster::isLonLat(x)
#   #  proj <- projection(x)
#   # tab <- spbabel::sptable(x)
#   gm <- gibble::gibble(x)
#   vertex <- silicate::sc_coord(x)
#   xy <- as.matrix(coords)
#   kd <- nabor::knn(xy, raster::as.matrix(vertex[, c("x_", "y_")]), k = 1, eps = 0)
#   index <- expand.grid(x = seq(ncol(coords)), y = rev(seq(nrow(coords))))[kd$nn.idx, ]
#   vertex$x_ <- index$x
#   vertex$y_ <- index$y
#   #  spbabel::sp(tab, attr_tab = as.data.frame(x), crs = NA_character_)
#   out <- tibble::as_tibble(as.data.frame(x))
#   out[[attr(x, "sf_column")]]  <- silicate:::build_sf(gm, vertex, crs = NA)
#   out
# }

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
romsmap.default <- function(x, coords, crop = FALSE, lonlat = TRUE,...) {
  x <- as.matrix(x)
  obj <- sp::SpatialPointsDataFrame(SpatialPoints(x, proj4string = sp::CRS("+init=epsg:4326")), data.frame(a = 1:nrow(x)))
  romsmap(obj, coords)
}

#' @rdname romsmap
#' @export
romsmap.SpatialLinesDataFrame <- romsmap.SpatialPolygonsDataFrame

#' @rdname romsmap
#' @export
romsmap.SpatialPointsDataFrame <- romsmap.SpatialPolygonsDataFrame

keeponlymostcomplexline <- function (x) 
{
  for (iObj in seq_len(nrow(x))) {
    if (inherits(x, "SpatialLinesDataFrame")) {
      wmax <- which.max(sapply(x[iObj, ]@lines[[1]]@Lines, 
                               function(x) nrow(x@coords)))
      x@lines[[iObj]]@Lines <- x@lines[[iObj]]@Lines[wmax]
    }
  }
  x
}
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
  polys <- sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(raster::as.matrix(cds)[unique(c(left, bottom, right, top)), ])), "1")), 
                                                    proj4string = sp::CRS("+init=epsg:4326")), 
                                    data.frame(boundary = "coords", stringsAsFactors = FALSE))
}

#' @name romsboundary
#' @export
databoundary <- function(x, mask = NULL, ...) {
  if (is.null(mask)) {
    ## longshot, put an NA band around the mask and contour on that
    maskboundary <- rgeos::gPolygonize(methods::as(keeponlymostcomplexline(
      raster::rasterToContour(is.na(raster::extend(x, c(2, 2), value = NA)), level = 1)), 
      "SpatialLinesDataFrame"))
    maskboundary <- SpatialPolygonsDataFrame(maskboundary, data.frame(boundary = "mask", stringsAsFactors = FALSE))
  }
  maskboundary
}

## put any raster into xy-index space (0, nc, 0, nr)
set_indextent <- function(x) {
  setExtent(x, extent(0, ncol(x), 0, nrow(x)))
}

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