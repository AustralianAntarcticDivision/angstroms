#' Remap an object to the space defined by coordinate arrays. 
#' 
#' Find the nearest-neighbour coordinates of `x` in the coordinate arrays of `coords`. 
#' 
#' The input `coords` is a assumed to be a 2-layer RasterStack or RasterBrick and
#' using `nabor::knn` the nearest matching position of the coordinates of `x` is found in the grid space of `coords`. The
#' motivating use-case is the curvilinear longitude and latitude arrays of ROMS model output. 
#' 
#' Cropping is complicated more details . . .
#' No account is made for the details of a ROMS cell, though this may be included in future. We tested only with the "lon_u" and "lat_u"
#' arrays. 
#' @param x object to transform to the grid space, e.g. a \code{\link[sp]{Spatial}} object
#' @param coords romscoords RasterStack
#' @param crop logical, if \code{TRUE} crop x to the extent of the boundary of the values in coords
#' @param lonlat logical, if \code{TRUE} check for need to back-transform to longitude/latitude and do it
#' @param ... unused
#'
#' @return input object with coordinates transformed to space of the coords 
#' @export
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
  x <- raster::intersect(x, boundary(coords))
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
  spbabel::sp(tab, attr_tab = as.data.frame(x), crs = proj)
}

#' @rdname romsmap
#' @export
romsmap.SpatialLinesDataFrame <- romsmap.SpatialPolygonsDataFrame

#' @rdname romsmap
#' @export
romsmap.SpatialPointsDataFrame <- romsmap.SpatialPolygonsDataFrame



#' Extract coordinate arrays from ROMS. 
#' 
#' Returns a RasterStack of the given variable names. 
#'
#' @param x ROMS file name
#' @param spatial names of coordinate variables (e.g. lon_u, lat_u) 
#' @param ncdf default to NetCDF no matter what file name
#' @param ... unused
#'
#' @return RasterStack with two layers of the 2D-variables
#' @export 
#'
#' @examples
#' \dontrun{
#'   coord <- romscoord("roms.nc")
#' }
#' @importFrom raster stack
romscoords <- function(x, spatial = c("lon_u", "lat_u"), ncdf = TRUE,  ... ) {
  l <- vector("list", length(spatial))
  for (i in seq_along(l)) l[[i]] <- raster(x, varname = spatial[i], ncdf = TRUE, ...)
  l <- lapply(l, function(x) setExtent(x, extent(0, nrow(x), 0, ncol(x))))
  raster::stack(l)
}


#' Coordinates at depth
#' 
#' \code{S} and \code{h} are the  names of the appropriate variables
#'
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




