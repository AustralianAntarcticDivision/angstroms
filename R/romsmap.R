# reverse and detect functions taken from spbabel
#' @importFrom sp Line Lines SpatialLinesDataFrame SpatialLines SpatialPolygonsDataFrame SpatialPointsDataFrame SpatialPointsDataFrame 
geomnames <- function () 
{
  list(SpatialPolygonsDataFrame = c("object_", "branch_", "island_", 
                                    "order_", "x_", "y_"), SpatialLinesDataFrame = c("object_", 
                                                                                     "branch_", "order_", "x_", "y_"), SpatialMultiPointsDataFrame = c("branch_", 
                                                                                                                                                       "object_", "x_", "y_"), SpatialPointsDataFrame = c("object_", 
                                                                                                                                                                                                          "x_", "y_"))
}
reverse_geomPoly <- function(x, d, proj) {
  objects <- split(x, x$object_)
  ## remove those columns used by reconstruction?
  d$branch_ <- d$object_ <- d$island_ <- d$order_ <- d$x_ <- d$y_ <- NULL
  if (ncol(d) < 1) d[["rownumber_"]] <- seq(nrow(d))
  if (ncol(d) < 1L) d$id_ <- seq(nrow(d))  ## we might end up with no attributes
  ## match.ID should be replaced by method to carry the original rownames somehow
  SpatialPolygonsDataFrame(SpatialPolygons(lapply(objects, loopBranchPoly), proj4string = CRS(proj)), d, match.ID = FALSE)
}
loopBranchPoly <- function(a) Polygons(lapply(dropZeroRowFromList(split(a, a$branch_)),
                                              function(b) Polygon(as.matrix(b[, c("x_", "y_")]), hole = !b$island_[1L])), as.character(a$object_[1L]))


reverse_geomLine <- function(x, d, proj) {
  objects <- split(x, x$object_)
  d$branch_ <- d$object_ <- d$order_ <- d$x_ <- d$y_ <- NULL
  if (ncol(d) < 1L) d$rownumber_ <- seq(nrow(d))  ## we might end up with no attributes
  SpatialLinesDataFrame(SpatialLines(lapply(objects, loopBranchLine), proj4string = CRS(proj)), d, match.ID = FALSE)
}
dropZeroRowFromList <- function(x) x[unlist(lapply(x, nrow), use.names = FALSE) > 0L]
loopBranchLine<- function(a) Lines(lapply(dropZeroRowFromList(split(a, a$branch_)), function(b) Line(as.matrix(b[, c("x_", "y_")]))), as.character(a$object_[1L]))

reverse_geomPoint <- function(a, d, proj) {
  # stop("not implemented")
  ## the decomposition is not yet applied for Multipoints . . .
  ## if (length(unique(a$object)) > 1) warning("no support for Multipoints yet")
  spts <- SpatialPoints(as.matrix(a[, c("x_", "y_")]))
  d$object_ <- d$x_ <- d$y_ <- NULL
  if (ncol(d) < 1) d[["rownumber_"]] <- seq(nrow(d))
  SpatialPointsDataFrame(spts, d, proj4string = CRS(proj))
}
#' @importFrom sp SpatialMultiPointsDataFrame SpatialMultiPoints
reverse_geomMultiPoint <- function(a, d, proj) {
  d$branch_ <- d$object_  <- d$x_ <- d$y_ <- NULL
  if (ncol(d) < 1L) d$rownumber_ <- seq(nrow(d))  ## we might end up with no attributes
  
  SpatialMultiPointsDataFrame(SpatialMultiPoints(lapply(split(a[, c("x_", "y_")], a$object_), as.matrix)), d, proj4string = CRS(proj))
}
detectSpClass <- function(x) {
  #if ("topol_" %in% names(x)) return(topol2sp(x$topol_))
  gn <- geomnames()
  for (i in seq_along(gn)) {
    if (all(gn[[i]] %in% names(x))) return(names(gn)[i])
  }
  
  #if (all(gn$SpatialPolygonsDataFrame %in% names(x))) return("SpatialPolygonsDataFrame")
  #if (all(gn$SpatialLinesDataFrame %in% names(x))) return("SpatialLinesDataFrame")
  #if (all(gn$SpatialPointsDataFrame %in% names(x))) return("SpatialPointsDataFrame")
  #if (all(gn$SpatialMultiPointsDataFrame %in% names(x))) return("SpatialMultiPointsDataFrame")
  # cat("cannot find matching topology type from these columns")
  #print(x)
  stop('cannot create Spatial* object from this input, matching topology type from these columns')
  
}


#' @importFrom rlang .data
.sp_from_table <- function (x, attr_tab = NULL, crs, ..., topol_ = NULL) 
{
  if (missing(crs)) 
    crs <- attr(x, "crs")
  if (is.null(crs)) 
    crs <- NA_character_
  if (is.null(topol_)) 
    target <- detectSpClass(x)
  else target <- topol_
  minc <- c(SpatialPolygonsDataFrame = 3, SpatialLinesDataFrame = 2, 
            SpatialMultiPointsDataFrame = 1, SpatialPointsDataFrame = 1)[target]
  if (nrow(x) < minc) 
    stop(sprintf("target is %s but input table has  %i %s", 
                 target, nrow(x), c("rows", "row")[(nrow(x) == 1) + 
                                                     1]))
  dat <- dplyr::distinct(x, .data$object_, .keep_all = TRUE)
  n_object <- nrow(dat)
  n_attribute <- nrow(attr_tab)
  if (is.null(n_attribute)) 
    n_attribute <- n_object
  if (!(n_attribute == n_object)) 
    stop("number of rows in attr must match distinct object in x")
  if (!is.null(attr_tab)) 
    dat <- dplyr::bind_cols(dat, attr_tab)
  gom <- switch(target, SpatialPolygonsDataFrame = reverse_geomPoly(x, 
                                                                    dat, crs), SpatialLinesDataFrame = reverse_geomLine(x, 
                                                                                                                        dat, crs), SpatialPointsDataFrame = reverse_geomPoint(x, 
                                                                                                                                                                              dat, crs), SpatialMultiPointsDataFrame = reverse_geomMultiPoint(x, 
                                                                                                                                                                                                                                              dat, crs))
  gom
}


#' Remap an object to the ROMS grid. 
#' 
#' Find the nearest-neighbour coordinates of `x` in the coordinate arrays of `coords`. 
#' 
#' The input `coords` is a assumed to be a 2-layer RasterStack or RasterBrick
#' and using [FNN::get.knnx()] the nearest matching position of the coordinates
#' of `x` is found in the grid space of `coords`. The motivating use-case is the
#' curvilinear longitude and latitude arrays of ROMS model output.
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
#' library(raster)
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
  sc <- silicate::PATH0(x)
  ind <- do.call(rbind, sc$object$path_)
  ind[["x_"]] <- silicate::sc_vertex(x)[["x_"]][ind[["vertex_"]]]
  ind[["y_"]] <- silicate::sc_vertex(x)[["y_"]][ind[["vertex_"]]]
  
  if (inherits(x, "SpatialPolygonsDataFrame")) {
  tab <- dplyr::transmute(ind, .data$x_, .data$y_, 
                          .data$object_, 
                          branch_ = .data$path_, 
                          island_ = TRUE, order_ = dplyr::row_number())
  #tab <- spbabel::sptable(x)
  }
  if (inherits(x, "SpatialLinesDataFrame")) {
    tab <- dplyr::transmute(ind, .data$x_, .data$y_, 
                            .data$object_, 
                            branch_ = .data$path_, 
                            order_ = dplyr::row_number())
    
  }
  if (inherits(x, "SpatialPointsDataFrame")) {
    tab <- dplyr::transmute(ind, .data$x_, .data$y_, 
                            .data$object_ 

                            )
    
  }
  
  if (repro & !is.na(proj)) {
    llproj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
    ## replace with reproj
    xy <- reproj::reproj(cbind(tab$x_, tab$y_, 0), target = llproj, source = proj)
  ## spbabel standard is attributes with underlines
    tab$x_ <- xy[,1]
    tab$y_ <- xy[,2]
    proj <- llproj
  }
  xy <- as.matrix(coords)
  kd <- FNN::get.knnx(xy, raster::as.matrix(tab[, c("x_", "y_")]), k = 1)
  index <- expand.grid(x = seq(ncol(coords)), y = rev(seq(nrow(coords))))[kd$nn.index[,1L,drop=TRUE], ]
 tab$x_ <- index$x
  tab$y_ <- index$y
  .sp_from_table(tab, attr_tab = as.data.frame(x), crs = NA_character_)
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

#' @param x a thing
#' @param mask logical
#' @param ... dots
#' 
#' @name romsboundary
#' @export
databoundary <- function(x, mask = NULL, ...) {
  stop("databoundary is defunct")
  # if (is.null(mask)) {
  #   ## longshot, put an NA band around the mask and contour on that
  #   maskboundary <- rgeos::gPolygonize(methods::as(keeponlymostcomplexline(
  #     raster::rasterToContour(is.na(raster::extend(x, c(2, 2), value = NA)), level = 1)), 
  #     "SpatialLinesDataFrame"))
  #   maskboundary <- sp::SpatialPolygonsDataFrame(maskboundary, data.frame(boundary = "mask", stringsAsFactors = FALSE))
  # }
  # maskboundary
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