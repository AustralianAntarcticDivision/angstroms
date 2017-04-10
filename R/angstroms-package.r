#' Tools for ROMS model output.
#' 
#'  
#' @name angstroms
#' @aliases angstroms-package
#' @docType package
#' @importFrom graphics lines plot points
#' @importFrom stats filter
NULL

#' Fake model data. 
#' 
#' `ice_coords` and `ice_fake` are generated from a projected map of southern Ocean sea ice data. 
#' 
#' The coords layer is the longitude and latitude values for the centres of the polar cells. This is veyr loosely
#' analogous to the coordinate arrays used by ROMS data, included here for working examples, illustration and code tests.
#' 
#' The proper metadata for these layers is 
#' "-3950000, 3950000, -3950000, 4350000  (xmin, xmax, ymin, ymax)"
#' 
#' "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs" 
#' @name ice_fake
#' @aliases ice_coords
#' @docType data
NULL

#' Antarctica simple coastline. 
#' 
#' Taken from "rnaturalearth::countries110"
#' @name antarctica
#' @docType data
NULL