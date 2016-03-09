#' Remap an object to the space defined by coordinate arrays. 
#' 
#' Find the nearest-neighbour coordinates of `x` in the coordinate arrays of `coords`. 
#' 
#' The input `coords` is a assumed to be a 2-layer \code{\link[raster]{RasterStack}} or \code{\link[raster]{RasterBrick}} and
#' using `nabor::knn` the nearest matching position of the coordinates of `x` is found in the grid space of `coords`. The
#' motivating use-case is the curvilinear longitude and latitude arrays of ROMS model output. 
#' 
#' Cropping is complicated more details . . .
#' No account is made for the details of a ROMS cell, though this may be included in future. We tested only with the "lon_u" and "lat_u"
#' arrays. 
#' @param x object to transform to the grid space, e.g. a \code{\link[sp]{Spatial}} object
#' @param coords romscoords RasterStack
#' @param crop logical, if \code{TRUE} crop x to the extent of the boundary of the values in coords
#' @param longlat logical, if \code{TRUE} check for need to back-transform to longitude/latitude and do it
#' @param ... 
#'
#' @return input object with coordinates transformed to space of the coords 
#' @export
romsmap <- function(x, ...) {
  UseMethod("romsmap")
}

#' @rdname romsmap
#' @export
#' @importFrom spbabel sptable spFromTable
#' @importFrom nabor knn
#' @importFrom raster intersect as.matrix projection
#' @importFrom sp CRS
romsmap.SpatialPolygonsDataFrame <- function(x, coords, crop = FALSE, lonlat = TRUE, ...) {
  ## first get the intersection
  if (crop) {
  op <- options(warn = -1)
  x <- raster::intersect(x, oms:::boundary(coords))
  options(op)
  }
  ## do we need to invert projection?
  repro <- !raster::isLonLat(x)
  proj <- projection(x)
  tab <- spbabel::sptable(x)
  
  if (repro & !is.na(proj)) {
    llproj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
    xy <- proj4::ptransform(cbind(tab$x, tab$y, 0), src.proj = proj, dst.proj = llproj, silent = FALSE) * 180 / pi
    tab$x <- xy[,1]
    tab$y <- xy[,2]
    proj <- llproj
  }
  xy <- as.matrix(coords)
  kd <- nabor::knn(xy, raster::as.matrix(tab[, c("x", "y")]), k = 1, eps = 0)
  index <- expand.grid(x = seq(ncol(coords)), y = rev(seq(nrow(coords))))[kd$nn.idx, ]
  tab$x <- index$x
  tab$y <- index$y
  spbabel::spFromTable(tab, crs = NA_character_)
}

#' @rdname romsmap
#' @export
romsmap.SpatialLinesDataFrame <- romsmap.SpatialPolygonsDataFrame

#' @rdname romsmap
#' @export
romsmap.SpatialPointsDataFrame <- romsmap.SpatialPolygonsDataFrame

## this is from rastermesh
boundary <- function(cds) {
  left <- cellFromCol(cds, 1)
  bottom <- cellFromRow(cds, nrow(cds))
  right <- rev(cellFromCol(cds, ncol(cds)))
  top <- rev(cellFromRow(cds, 1))
  ## need XYFromCell method
  SpatialPolygons(list(Polygons(list(Polygon(raster::as.matrix(cds)[unique(c(left, bottom, right, top)), ])), "1")))
}


#' Extract coordinate arrays from ROMS. 
#' 
#' Returns a RasterStack of the given variable names. 
#'
#' @param x ROMS file name
#' @param spatial names of coordinate variables (e.g. lon_u, lat_u) 
#'
#' @return \code{\link[raster]{RasterStack}} with two layers of the 2D-variables
#' @export 
#'
#' @examples
#' \dontrun{
#'   coord <- romscoord("roms.nc")
#' }
romscoords <- function(x, spatial = c("lon_u", "lat_u")) {
  l <- vector("list", length(spatial))
  for (i in seq_along(l)) l[[i]] <- raster(x, varname = spatial[i])
  stack(l)
}


#' Coordinates at depth
#' \code{S} and \code{h} are the  names of the appropriate variables
#' @param x ROMS file name 
#' @param S  of S-coordinate stretching curve at RHO-points
#' @param h bathymetry at RHO-points
#' @return \code{\link[raster]{RasterStack}} with a layer for every depth
#' @export
romshcoords <- function(x, S = "Cs_r", depth = "h"){
  h <- raster(x, varname = depth)
  Cs_r <- ncget(x, S)
  v <- values(h)
  setExtent(brick(array(rep(rev(Cs_r), each = length(v)) * v, c(ncol(h), nrow(h), length(Cs_r))), transpose = TRUE), 
            extent(0, ncol(h), 0, nrow(h)))
}


ncget <- function(x, varname) {
  nc <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, varname)
}

#' Extract a data lyaer from ROMS by name and slice. 
#'
#' @param x ROMS file name
#' @param varname name of ROMS variable 
#' @param slice index in w and t (depth and time), defaults to first encountered
#'
#' @return \code{\link[raster]{RasterLayer}}
#' @export
#'
romsdata <-function(x, varname, slice = c(1, 1)) {
  brick(x, level = slice[1L], varname = varname)[[slice[2L]]]
}



#' Title
#'
#' @param varname 
#'
#' @return
#' @export
#'
#' @examples
ncdim <- function(x, varname) {
   library(rancid)
   roms <- NetCDF(x)
  # ## still exploring neatest way to do this . . .
   vdim <- vars(roms) %>% 
     filter(name == varname) %>% 
     inner_join(roms$vardim, "id") %>% 
     dplyr::transmute(id = dimids) %>% 
     inner_join(dims(roms), "id") 
   vdim$len
}


