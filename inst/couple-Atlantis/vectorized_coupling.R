## helper functions for model output NetCDF using R raster package
library(angstroms) # remotes::install_github("hypertidy/angstroms")
## read and parse functions for Atlantis BGM files
library(rbgm) # remotes::install_github("AustralianAntarcticDivision/rbgm")
library(dplyr)

## TODO
## - deepest_depth should be botz, check


##################################################################
## Utility functions for this process script

#' Clean up function for raster cache (shouldn't be necessary)
delif <- function(x) {
  f <- filename(x)
  if (file.exists(f)) {
    ff <- gsub("grd$", "gri", f)
    if (file.exists(ff)) file.remove(ff)
  }
  invisible(NULL)
}

#' Wrap angles (degree) into -180,180
wrap180 <- function(x, lmin = -180) (x - lmin) %% 360 + lmin

#' sign function (never 0, so simple multiply)
sgn1 <- function(x) ifelse(x < 0, -1, 1)

#' Get raster values, explicitly pull from disk if necessary
rvalues <- function(x) {
  if (nchar(raster::filename(x)) > 0) x <- readAll(x)
  as.vector(raster::values(x))
}


####################################################################
## prep 

## get local ROMS files for his_31
cpolar <- raadtools:::cpolarfiles() %>% 
  dplyr::filter(stringr::str_detect(fullname, "his_31"))

## convert files list to time step form, filename, and "dim4_slice"
file_db <- purrr::map_df(cpolar$fullname, 
                         function(x) tibble::tibble(fullname = x,
                                                    time = tidync::tidync(x)$transforms$ocean_time$ocean_time) %>% 
                           dplyr::mutate(dim4_slice = row_number()))

## ignore the file time steps and assume they daily from the start value
## absolute time is not relevant to Atlantis (yet)
file_db$time_steps <-  seq(ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC") + 
                             file_db$time[1], by = "1 day", length.out = nrow(file_db))
## a constant
time_step <- 86400

## get a BGM and read it
## - box geometry model for Atlantis
## - it's an edge list, line segments (faces), defining boxes (2D) and faces (1D) that extend
## those boxes to multiple depths
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- rbgm:::read_bgm(bfile)

## we need a dummy file to get into things, build indexes
roms_file <- file_db$fullname[1]

## this is the 2D coordinates (longlat) of ROMS
roms_ll <- romscoords(roms_file, transpose = TRUE)
## boxes of BGM in sp form (a local LCC projection)
boxes <- boxSpatial(bgm)

## width and height of each box (not used currently but no cost)
e <- vector("list", nrow(boxes))
for (i in seq_along(e)) {ex <- raster::extent(boxes[i, ]); e[[i]] <- c(xmin = xmin(ex), xmax = xmax(ex), ymin = ymin(ex), ymax = ymax(ex))}
boxes$width <- apply(do.call(rbind, e)[,1:2], 1, diff)
boxes$height <- apply(do.call(rbind, e)[,3:4], 1, diff)

## faces of BGM in sp form 
faces <- faceSpatial(bgm)
## convert faces and boxes to longlat (so we can use angular trig in ROMS coordinates)
faces_ll <- spTransform(faces, "+init=epsg:4326")
boxes_ll <- spTransform(boxes, "+init=epsg:4326")
faces$.x0 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][1, 1])
faces$.x1 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][2, 1])
faces$.y0 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][1, 2])
faces$.y1 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][2, 2])

## map the bgm into ROMS grid (1:ncol, 1:nrow)
romsbox <- romsmap(boxes, roms_ll)
roms_ll <- crop(roms_ll , romsbox, snap = "out")
romsface <- romsmap(faces, romscoords(roms_file, transpose = TRUE))

## dummy grids from the ROMS (later we iterate over slice, which is time)
u <- crop(romsdata3d(roms_file, varname = "u", slice = 1), romsbox)
v <- crop(romsdata3d(roms_file, varname = "v", slice = 1), romsbox)
temp <- crop(romsdata3d(roms_file, varname = "temp", slice = 1), romsbox)
salt <- crop(romsdata3d(roms_file, varname = "salt", slice = 1), romsbox)
vert <- subset(crop(romsdata3d(roms_file, varname = "w", slice = 1), romsbox), 1:31)
h <- crop(romsdata(roms_file, varname = "h"), romsbox, snap = "out")

## every ROMS cell depth (though we don't use it currently, we want this finite element
## approach for comparison/completeness)
hhh <- crop(romsdepth(roms_file), romsbox, snap = "out")
## snap h onto the first layer (the bottom) and calculate delta_rho (dz)
delta_rho <- calc(brick(-h, hhh), fun = function(x) diff(x))


## get a reference depth for the bottom (should be botz not max(abs(h)) - let's see)
deepest_depth <- max(raster::extract(h, tabularaster::cellnumbers(h, sf::st_as_sf(romsbox))$cell_), na.rm = TRUE)
## here we must use the same as we used for the nc output files
dz <- rev(rbgm::build_dz(-deepest_depth, zlayers = c(-Inf,  -2000, -1000, -750, -400, -300, 
                                                     -200, -100, -50, -20, 0)))
atlantis_depths <- head(-cumsum(c(0, dz)), -1)

## now we are ready to build the index of ROMS cells
## (ignore dz here because we aren't using ROMS cells for their height, only their
## position in the layer stack)
romstab <- tabularaster::as_tibble(u, dim = FALSE) %>% 
  dplyr::rename(u = cellvalue) %>% dplyr::mutate(v = rvalues(v), 
                                                 temp = rvalues(temp), 
                                                 salt = rvalues(salt), 
                                                 vert = rvalues(vert),
                                                 depth = rvalues(hhh), 
                                                 #dz = rvalues(delta_rho),  ## not this, this is ROMS cell depth, not face delta
                                                # harea = rvalues(harea),
                                                 layer = findInterval(-depth, -atlantis_depths), 
                                                delta_layer = dz[layer])
delif(delta_rho)

## now the box and face indexes
box_index <- tabularaster::cellnumbers(u, romsbox)
box_index$boxid <- romsbox$box_id[box_index$object_]
## keep botz as we use it to filter too-low ROMS data
box_index$maxlayer  <- findInterval(-romsbox$botz[box_index$object_], -atlantis_depths)

## FIXME clean this up, it's a mess due to iterative fixes/experimentation
face_index <- tabularaster::cellnumbers(u, romsface)
face_index$faceid <- romsface$.fx0[face_index$object_]
face_index$object_ <- NULL
face_index <- inner_join(face_index, faces@data[c("length", ".fx0", ".x0", ".y0", ".x1", ".y1")], c("faceid" = ".fx0"))
face_index <- face_index %>% 
  mutate(angle_face = (180 * atan2(.x0 - .x1, .y0 - .y1))/pi)   # angle in radians of tan(y/x) / pi resulting in angle in degrees
face_index$.x0 <- NULL
face_index$.x1 <- NULL
face_index$.y0 <- NULL
face_index$.y1 <- NULL

### FOR TRANSPORT/MASS NC FILES
transp_filename <- sprintf("%s_transport.nc", tempfile())
mass_filename  <- sprintf("%s_mass.nc", tempfile())
transp_params <- ""
mass_params <- ""

## required angle for rotating ROMS vectors into longlat space
angle <- crop(angstroms:::romsangle(roms_file), romsbox)

## create netcdf files (do it here to ensure var shape is reliable)
## but FIXME: clean this up
angstroms:::create_transport(transp_filename, model_title = "Transport file Antarctica_28", bgmfilepath = bfile, 
                             bgmlevels = atlantis_depths, time_steps = file_db$time_steps)

angstroms:::create_mass(mass_filename, model_title = "Mass file Antarctica_28", bgmfilepath = bfile, 
                        bgmlevels = atlantis_depths, time_steps = file_db$time_steps)

## get variables from the files, so we know they reliable (a paranoid leftover FIXME)
vertical <- salinity <- temperature <- angstroms::rawdata(mass_filename, "temperature")
transport <- angstroms::rawdata(transp_filename, "transport")

# system.time(doitfun(1))
# user  system elapsed 
# 11.264   2.105  13.465 
# system.time(doitfun(300))
# user  system elapsed 
# 10.214   1.342  16.637 

## a trivial function to run a single iteration, all is a single time
doitfun <- function(itime) {
  roms_file <- file_db$fullname[itime]
  roms_slice <- file_db$dim4_slice[itime]
  u0 <- crop(romsdata3d(roms_file, varname = "u", slice = roms_slice), romsbox)
  v0 <- crop(romsdata3d(roms_file, varname = "v", slice = roms_slice), romsbox)
  
  ## apply required ROMs rotation  https://www.myroms.org/forum/viewtopic.php?f=3&t=295
  ul <- vector("list", nlayers(u0))
  vl <- vector("list", nlayers(u0))
  for (j in seq_along(ul)) {
    uv <- romsrotate(brick(u0[[j]], v0[[j]]), angle)
    ul[[j]] <- uv[[1]]
    vl[[j]] <- uv[[2]]
  }
  u0 <- brick(ul)
  v0 <- brick(vl)
  for (j in seq_along(ul)) {
    delif(ul[[j]])
    delif(vl[[j]])
  }
  rm(ul, vl)
  temp0 <- crop(romsdata3d(roms_file, varname = "temp", slice = roms_slice), romsbox)
  salt0 <- crop(romsdata3d(roms_file, varname = "salt", slice = roms_slice), romsbox)
  vert0 <- subset(crop(romsdata3d(roms_file, varname = "w", slice = roms_slice), romsbox), 1:31)
  
  romstab$u <- rvalues(u0)
  romstab$v <- rvalues(v0)
  romstab$temp <- rvalues(temp0)
  romstab$salt <- rvalues(salt0)
  romstab$vert <- rvalues(vert0)
  delif(u0)
  delif(v0)
  delif(temp0)
  delif(salt0)
  delif(vert0)
  
box_summ <- romstab %>% 
  inner_join(box_index[c("cell_", "boxid", "maxlayer")], c("cellindex" = "cell_")) %>% 
  dplyr::mutate(layer = pmin(layer, maxlayer)) %>%  ## push lower data into lowest layer
  group_by(boxid, layer) %>% summarize(temp = mean(temp, na.rm = TRUE), 
                                       salt = mean(salt, na.rm = TRUE), 
                                       vert = sum(vert, na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(boxid, layer)


face_summ <- romstab %>% dplyr::select(-temp, -salt) %>% 
  inner_join(face_index, c("cellindex" = "cell_")) %>% 
  inner_join(box_index[c("cell_", "maxlayer")], c("cellindex" = "cell_")) %>% 
  dplyr::mutate(layer = pmin(layer, maxlayer)) %>%  ## push lower data into lowest layer
 ## mutate(angle_face = (180 * atan2(.x0 - .x1, .y0 - .y1))/pi)  %>%   # angle in radians of tan(y/x) / pi resulting in angle in degrees
  mutate(angle_flow = (180 * atan2(u, v))/pi) %>%
  mutate(relative_angle = wrap180(angle_flow - angle_face), 
         velocity_magnitude = sqrt(u^2 + v^2) * sgn1(relative_angle)) %>% 
  group_by(faceid, layer) %>% 
  ## mean because (i)  https://github.com/AustralianAntarcticDivision/EA-Atlantis-dev/issues/5#issuecomment-457453852
  summarize(flow = mean(velocity_magnitude * length * dz[layer], na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(faceid, layer)

 list(face = face_summ[c("faceid", "layer", "flow")], 
      box = box_summ[c("boxid", "layer", "temp", "salt", "vert")])
}


library(future)
plan(multiprocess)
## here I tried distributing the slices different but it makes no difference
x <- future.apply::future_lapply(1:nrow(file_db), doitfun)

## clean up the output from the time loop
for (itime in seq_along(x)) {
  face_summ <- x[[itime]]$face
  box_summ <- x[[itime]]$box
  ## FIXME loop over these for itime
  fidx <- cbind(face_summ$layer, face_summ$faceid + 1, itime)
  bidx <- cbind( box_summ$layer, box_summ$boxid + 1, itime)
  temperature[bidx] <- box_summ$temp
  salinity[bidx] <- box_summ$salt
  vertical[bidx] <- box_summ$vert
  transport[fidx] <-face_summ$flow
}




## output the arrays to the files created earlier
library(ncdf4)


ncmass <- ncdf4::nc_open(mass_filename, write = TRUE)
nctran <- ncdf4::nc_open(transp_filename, write = TRUE)


ncvar_put(nctran, "transport", transport)
ncvar_put(ncmass, "temperature", temperature)
ncvar_put(ncmass, "salinity", salinity)
ncvar_put(ncmass, "verticalflux", vertical)
nc_close(nctran)
nc_close(ncmass)

## copy to local dir for local workflow
file.copy(mass_filename, ".")
file.copy(transp_filename, ".")

