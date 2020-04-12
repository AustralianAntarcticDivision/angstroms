#define dimensions
transport_dimensions <- function(bgm, ntime = 365 * 2, nlevels = 2L, 
                                 time_unit = "seconds since 1970-01-01 00:00:00") {
  
  
  timedim   <- ncdim_def("time", "", seq_len(ntime), unlim = TRUE, create_dimvar = FALSE)
  leveldim  <- ncdim_def("level", "", seq_len(nlevels), create_dimvar = FALSE)
  facesdim  <- ncdim_def("faces", "", seq_len(bgm$extra$nface), create_dimvar = FALSE)
  
  var.time  <- ncvar_def("time",  time_unit,      dim = timedim,   prec = "double")
  var.face  <- ncvar_def("faces", "",             dim = facesdim,  longname="Face IDs", prec='integer')
  var.lev   <- ncvar_def("level", "",              dim = leveldim,  longname="layer index; 1=near surface", prec="integer")
  var.trans <- ncvar_def("transport", "m3/s",    dim = list(leveldim,facesdim,timedim), 0, prec="float")
  var.destb <- ncvar_def("dest_boxid","id",      dim = facesdim, longname="ID of destination box", prec="integer")
  var.srcb  <- ncvar_def("source_boxid","id",  facesdim,     longname="ID of source box",prec="integer")
  
  var.pt1x  <- ncvar_def("pt1_x", "degree_east", facesdim, longname = "x-coord of pt 1 of face", prec='float')
  var.pt2x  <- ncvar_def("pt2_x", "degree_east", facesdim, longname = "x-coord of pt 2 of face", prec='float')
  var.pt1y  <- ncvar_def("pt1_y", "degree_north", facesdim, longname = "y-coord of pt 1 of face", prec='float')
  var.pt2y  <- ncvar_def("pt2_y", "degree_north", facesdim, longname = "y-coord of pt 2 of face", prec='float')
  
  
  list(var.time, var.face, var.lev, var.destb, var.srcb, var.pt1x, var.pt2x, var.pt1y, var.pt2y, var.trans)
}

mass_dimensions <- function(bgm, ntime = 365 * 2, nlevels = 2L, 
                            time_unit = "seconds since 1970-01-01 00:00:00") {
  timedim   <- ncdim_def("time", "", seq_len(ntime), unlim = TRUE, create_dimvar = FALSE)
  leveldim  <- ncdim_def("level", "", seq_len(nlevels), create_dimvar = FALSE)
  boxesdim  <- ncdim_def("boxes", "", seq_len(bgm$extra$nbox), create_dimvar = FALSE)
  
  var.time <- ncvar_def("time", time_unit, timedim, prec="double")
  var.box  <- ncvar_def("boxes", "", boxesdim, longname="Box IDs", prec='integer')
  var.lev  <- ncvar_def("level", "", leveldim, longname="layer index; 1=near surface; positive=down" ,
                        prec="integer")
  var.vertflux <- ncvar_def("verticalflux", "m3/s", list(leveldim, boxesdim, timedim), -999,
                            longname = "vertical flux averaged over floor of box", prec = "float")
  var.temp     <- ncvar_def("temperature", "degree_C", list(leveldim, boxesdim, timedim), -999,
                            longname = "temperature volume averaged", prec = "float")
  var.salt     <- ncvar_def("salinity", "psu", list(leveldim,boxesdim,timedim), -999,
                            longname = "salinity volume averaged", prec = "float")
  
  list(var.time,var.box, var.lev, var.salt, var.temp, var.vertflux)
}

#' Create "transport" file for Atlantis. 
#'
#' Transport is the hydrodynamic forcing for the Atlantis model from oceanographic models. These functions
#' aim to minimize the amount of manual handling of details, creating an NetCDF file that can be modified directly. 
#' All the Box Geometry Model details are determined from the bgm file itself. The output file is the precursor file
#' for input to HydroConstruct (as per section 4.1 of the Atlantis User Guide). 
#' 
#' Please use the original BGM file in its native map projection, there is no need to copy this to a longitude-latitude version
#' of the same file. 
#' 
#' TODO: check that bgmlevels should only be an index, check that time_steps must be regular ?( or just monotonic incresing?)
#' @param filename the NetCDF file to create
#' @param model_title name of the model
#' @param bgmfilepath path to BGM file
#' @param bgmlevels the Atlantis depth levels
#' @param time_steps actual time steps, need to be regularly space
#' @param transp_params optional details to put in the NetCDF notes
#' @param overwrite set to \code{TRUE} to clobber an existing file
#' @return the filename of the output (use ncdf4 to inspect, modify it)
#' @export
#' @importFrom ncdf4 nc_open ncatt_put ncvar_put nc_close ncdim_def ncvar_def
#' @examples
#' tfile <- file.path(tempdir(), "my_hydrofile.nc")
#' (ncfile <- create_transport(tfile, bgmfilepath = bgmfiles::bgmfiles("antarctica_99")))
#' ## use ncdf4::nc_open to inspect the file
create_transport <- function(filename,
                             model_title = "Transport file [placeholder]", 
                             bgmfilepath = "",
                             bgmlevels = NULL, 
                             time_steps = NULL,
                             transp_params = "", overwrite = FALSE) {
  
  
  bgm <- rbgm::bgmfile(bgmfilepath)
  
  if (is.null(bgmlevels)) {
    deltas <- rbgm::build_dz(as.numeric(bgm$extra$maxwcbotz))
    bgmlevels <- cumsum(c(bgm$extra$maxwcbotz, deltas))
    message(sprintf("'bgmlevels' not supplied, building levels at %s", paste(bgmlevels, collapse = ", ")))
    
  } 
  nlevels <- length(bgmlevels)
  if (is.null(time_steps)) {
    time_steps <- seq(as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), by = 12 * 3600, length = 10)
  }
  
  stopifnot(inherits(time_steps, "POSIXct"))
  stopifnot(length(unique(diff(unclass(time_steps)))) == 1)
  
  dimensions <- transport_dimensions(bgm, nlevels = nlevels)
  
  if (file.exists(filename) && !overwrite) {
    stop(sprintf("'filename' already exists, use 'overwrite = TRUE' or delete: \n%s", filename))
  }
  nc_transp <-  ncdf4::nc_create(filename, dimensions)
  on.exit(ncdf4::nc_close(nc_transp), add = TRUE)
  #assign global attributes to file
  ncdf4::ncatt_put(nc_transp, 0, "title",      model_title)
  ncatt_put(nc_transp, 0, "geometry",   basename(bgmfilepath))
  ncatt_put(nc_transp, 0, "parameters", transp_params)
  #assign attributes to variables
  ncatt_put(nc_transp, "time", "dt", diff(time_steps[1:2]), prec="double") 
  ncdf4::ncvar_put(nc_transp, "time", time_steps, count = length(time_steps))
  ncvar_put(nc_transp, "faces",  bgm$faces$.fx0)  ## or bgm$faces$.fx0 ??
  ncvar_put(nc_transp, "level",   seq_along(bgmlevels), count = length(bgmlevels)) 
  ## from RyanM 
  # positive flow is right to left across face from p1 to p2
  ncvar_put(nc_transp, "dest_boxid",    bgm$faces$left)
  ncvar_put(nc_transp, "source_boxid",  bgm$faces$right)
  ll <- bgm$vertices
  if (!raster::isLonLat(bgm$extra$projection)) {
    ll[c("x", "y")] <- reproj::reproj(as.matrix( ll[c("x", "y")]), source = bgm$extra$projection, 
                                      target = "+proj=longlat +datum=WGS84")[,1:2]
  }
  v1 <- bgm$facesXverts[bgm$facesXverts[[".p0"]] == 1, ] %>%  dplyr::inner_join(ll, ".vx0")
  v2 <- bgm$facesXverts[bgm$facesXverts[[".p0"]] == 2, ] %>%  dplyr::inner_join(ll, ".vx0")
  

  ncvar_put(nc_transp, "pt1_x",  v1$x)
  ncvar_put(nc_transp, "pt2_x",  v2$x)
  ncvar_put(nc_transp, "pt1_y",  v1$y)
  ncvar_put(nc_transp, "pt2_y",  v2$y)
  filename
}


create_mass <- function(filename,
                             model_title = "Box averaged properties file [placeholder]", 
                             bgmfilepath = "",
                             bgmlevels = NULL, 
                             time_steps = NULL,
                             mass_params = "", overwrite = FALSE) {
  
  
  bgm <- rbgm::bgmfile(bgmfilepath)
  
  if (is.null(bgmlevels)) {
    deltas <- rbgm::build_dz(as.numeric(bgm$extra$maxwcbotz))
    bgmlevels <- cumsum(c(bgm$extra$maxwcbotz, deltas))
    message(sprintf("'bgmlevels' not supplied, building levels at %s", paste(bgmlevels, collapse = ", ")))
    
  } 
  nlevels <- length(bgmlevels)
  if (is.null(time_steps)) {
    time_steps <- seq(as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), by = 12 * 3600, length = 10)
  }
  
  stopifnot(inherits(time_steps, "POSIXct"))
  stopifnot(length(unique(diff(unclass(time_steps)))) == 1)
  
  dimensions <- mass_dimensions(bgm, nlevels = nlevels)
  
  if (file.exists(filename) && !overwrite) {
    stop(sprintf("'filename' already exists, use 'overwrite = TRUE' or delete: \n%s", filename))
  }
  nc_varfile <-  ncdf4::nc_create(filename, dimensions)
  on.exit(ncdf4::nc_close(nc_varfile), add = TRUE)
  

  #assign global attributes to file
  ncatt_put(nc_varfile,0,"title", model_title)
  ncatt_put(nc_varfile, 0, "geometry", bgmfilepath)
  ncatt_put(nc_varfile, 0, "parameters", mass_params)

  #assign attributes to variables
  ncatt_put(nc_varfile, "time","dt",86400,prec="double")

  ncvar_put(nc_varfile, "time", time_steps, count = length(time_steps))
  ncvar_put(nc_varfile, "level", bgmlevels)
  ncvar_put(nc_varfile, "boxes", bgm$boxes$.bx0)
  filename
}