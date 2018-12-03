##devtools::install_github(c("AustralianAntarcticDivision/angstroms"))
library(angstroms)
library(rbgm) ## read BGM
library(bgmfiles) ## archive of BGM files
library(raadtools)
library(tidync)
library(dplyr)


## TODO
## rotate u/v vectors on import (see ROMS details)
## check hyper diffusion scaling
## check depth interval (do we multiply by face height?)
## apply interpolation between ROMs layers (per pixel for the pre-stage)
## check vertical flux calcs


cpolar <- raadtools:::cpolarfiles()
  
## for each file we know about each separate time, and this is lvar = 4, level = [dim4_slice]
## there are 31 depth levels (surface water is the first one)
## files have 28, 30, or 31 time steps file_db %>% group_by(fullname) %>% tally() %>% distinct(n)
## so we can't just assume a constant

## we are expected to read  via romsdata3d(fullname[itime], slice =  dim4_slice[itime], lvar = 4)
## as fullname is replicated out for each time step ()
file_db <- purrr::map_df(cpolar$fullname, 
                         function(x) tibble::tibble(fullname = x,
                                                    time = tidync(x)$transforms$ocean_time$ocean_time) %>% 
                           dplyr::mutate(dim4_slice = row_number()))

time_step <- 3600 * 24
## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- rbgm:::read_bgm(bfile)


library(dplyr)


clamp <- function(x) {
  x[x < 1] <- 1
  x[x > 31] <- 31
  x
}


itime <- 1
roms_file <- file_db$fullname[itime]
roms_slice <- file_db$dim4_slice[itime]
roms_ll <- romscoords(roms_file, transpose = TRUE)
boxes <- boxSpatial(bgm)
## width and height of each box
e <- vector("list", nrow(boxes))
for (i in seq_along(e)) {ex <- raster::extent(boxes[i, ]); e[[i]] <- c(xmin = xmin(ex), xmax = xmax(ex), ymin = ymin(ex), ymax = ymax(ex))}
boxes$width <- apply(do.call(rbind, e)[,1:2], 1, diff)
boxes$height <- apply(do.call(rbind, e)[,3:4], 1, diff)
faces <- faceSpatial(bgm)
romsbox <- romsmap(boxes, roms_ll)
roms_ll <- crop(roms_ll , romsbox, snap = "out")
romsface <- romsmap(faces, romscoords(roms_file, transpose = TRUE))
## dummy grids
u <- crop(romsdata3d(roms_file, varname = "u", slice = 1), romsbox)
v <- crop(romsdata3d(roms_file, varname = "v", slice = 1), romsbox)
temp <- crop(romsdata3d(roms_file, varname = "temp", slice = 1), romsbox)
salt <- crop(romsdata3d(roms_file, varname = "salt", slice = 1), romsbox)

h <- crop(romsdata(roms_file, varname = "h"), romsbox, snap = "out")
hhh <- crop(romsdepth(roms_file), romsbox, snap = "out")

## Cs_r is the S-coord stretching
Cs_r <- rawdata(roms_file, "Cs_r")
#list_nc_z_index <- vector('list', nrow(box_roms_index))
deepest_depth <- max(raster::extract(h, tabularaster::cellnumbers(h, sf::st_as_sf(romsbox))$cell_), na.rm = TRUE)
## here we must use the same as we used for the nc output files
dz <- rev(rbgm::build_dz(-3000))
atlantis_depths <- -cumsum(c(0, dz))

## boxes
box_mass_index <- matrix(NA_integer_, nrow(hhh@data@values), length(atlantis_depths))
for (i in seq_len(nrow(box_mass_index))) {
  idx <- clamp(findInterval(atlantis_depths, hhh@data@values[i, ]))
  ## max out at ~200m beyond the deepest depth
  idx[(atlantis_depths - min(hhh@data@values[i, ])) < -20] <- NA_integer_
  box_mass_index[i, ] <- idx
  if (i %% 1000 == 0) print(i)
}

## faces
face_index <- matrix(NA_integer_, nrow(u@data@values), 
                     length(atlantis_depths))


for (i in seq_len(nrow(face_index))) {
  idx <- clamp(findInterval(atlantis_depths, hhh@data@values[i, ]))
  ## max out at ~200m beyond the deepest depth
  idx[(atlantis_depths - min(hhh@data@values[i, ])) < -20] <- NA_integer_
  face_index[i, ] <- idx
  if (i %% 1000 == 0) print(i)
}


## box properties
temp <- crop(romsdata3d(roms_file, varname = "temp", slice = 1), romsbox)
salt <- crop(romsdata3d(roms_file, varname = "salt", slice = 1), romsbox)
salt_atlantis <- temp_atlantis <- box_mass_index * NA_real_
for (i in seq_len(nrow(temp_atlantis))) {
  temp_atlantis[i, ] <- temp@data@values[i, box_mass_index[i, ]]
  salt_atlantis[i, ] <- salt@data@values[i, box_mass_index[i, ]]
  
  if (i %% 1000 == 0) print(i)
}





## create NetCDF templates for Atlantis hydro
library(ncdf4)

### FOR TRANSPORT NC FILE
transp_filename <- sprintf("%s_transport.nc", tempfile())
mass_filename  <- sprintf("%s_mass.nc", tempfile())
transp_params <- ""
mass_params <- ""

bgmfilepath <- bgmfiles::bgmfiles("antarctica_28")
library(angstroms)


## there are 1095 daily time steps

# ocean_time  Size:31   *** is unlimited ***
#   long_name: time since initialization
# units: seconds since 0001-01-01 00:00:00
# calendar: julian
# field: time, scalar, series

## GRRR 
#time_steps <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC") + file_db$time
## ignore the file time steps and assume they daily from the start value
time_steps <-  seq(ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC") + file_db$time[1], by = "1 day", length.out = nrow(file_db))


angstroms:::create_transport(transp_filename, model_title = "Transport file Antarctica_28", bgmfilepath = bgmfilepath, 
                 bgmlevels = atlantis_depths, time_steps = time_steps)

angstroms:::create_mass(mass_filename, model_title = "Mass file Antarctica_28", bgmfilepath = bgmfilepath, 
            bgmlevels = atlantis_depths, time_steps = time_steps
            )

ncmass <- ncdf4::nc_open(mass_filename, write = TRUE)
nctran <- ncdf4::nc_open(transp_filename, write = TRUE)
vertical <- salinity <- temperature <- ncvar_get(ncmass, "temperature")
transport <- ncvar_get(nctran, "transport")
cols <- palr::sstPal(2600)
delif <- function(x) {
  f <- filename(x)
  if (file.exists(f)) {
    ff <- gsub("grd$", "gri", f)
    if (file.exists(ff)) file.remove(ff)
  }
  invisible(NULL)
}
for (itime in seq_len(nrow(file_db))) {

  ## box properties
  temp <- crop(romsdata3d(file_db$fullname[itime], varname = "temp", slice = file_db$dim4_slice[itime]), romsbox)
  tempvals <- values(temp)
  delif(temp)
  rm(temp)
  salt <- crop(romsdata3d(file_db$fullname[itime], varname = "salt", slice = file_db$dim4_slice[itime]), romsbox)
  saltvals <- values(salt)
  delif(salt)
  rm(salt)
  salt_atlantis <- temp_atlantis <- box_mass_index * NA_real_
 

  ## face properties
  u <- crop(romsdata3d(file_db$fullname[itime], varname = "u", slice = file_db$dim4_slice[itime]), romsbox)
  uvals <- values(u)

  v <- crop(romsdata3d(file_db$fullname[itime], varname = "v", slice = file_db$dim4_slice[itime]), romsbox)
  vvals <- values(v)

  u_atlantis <- v_atlantis <- face_index * NA_real_

  for (i in seq_len(nrow(temp_atlantis))) {
    temp_atlantis[i, ] <- tempvals[i, box_mass_index[i, ]]
    salt_atlantis[i, ] <- saltvals[i, box_mass_index[i, ]]
    u_atlantis[i, ] <- uvals[i, face_index[i, ]]
    v_atlantis[i, ] <- vvals[i, face_index[i, ]]
  }
  temp0  <- setValues(subset(hhh, 1:ncol(temp_atlantis)), temp_atlantis)
  salt0 <- setValues(subset(hhh, 1:ncol(temp_atlantis)), salt_atlantis)
  names(temp0) <- sprintf("depth_%i", as.integer(abs(atlantis_depths)))
  names(salt0) <- names(temp0)
  u0  <- setValues(subset(u * NA_real_, 1:ncol(u_atlantis)), u_atlantis)
  v0  <- setValues(subset(v * NA_real_, 1:ncol(v_atlantis)), v_atlantis)
  delif(u)
  delif(v)
  rm(u, v)
  names(u0) <- sprintf("depth_%i", as.integer(abs(atlantis_depths)))
  names(v0) <- names(temp0)
  
  if (itime == 1) {
    mass_cells <- tabularaster::cellnumbers(temp0[[1]], sf::st_as_sf(romsbox))
    face_cells <- tabularaster::cellnumbers(u0[[1]], sf::st_as_sf(romsface))
    #face_cells$left <- boxes$box_id[match(faces$left[face_cells$object_], boxes$box_id)]
    #face_cells$right <- boxes$box_id[match(faces$right[face_cells$object_], boxes$box_id)]
    face_cells$box_width <- boxes$width[match(faces$left[face_cells$object_], boxes$box_id)]
    face_cells$box_height <- boxes$height[match(faces$left[face_cells$object_], boxes$box_id)]
    face_cells$length <- faces$length[face_cells$object_]
  }
 
  for (ilayer in seq_len(nrow(temperature))) {
    mass_cells$temperature <- raster::extract(temp0[[ilayer]], mass_cells$cell_)
    mass_cells$salinity <- raster::extract(salt0[[ilayer]], mass_cells$cell_)
    mass_cells$vertical <- raster::extract(sqrt(u0[[ilayer]] * u0[[ilayer]] + v0[[ilayer]] * v0[[ilayer]]), mass_cells$cell_)
    mass_summ <- mass_cells %>% dplyr::group_by(object_) %>% dplyr::summarize(temperature = mean(temperature, na.rm = TRUE), 
                                                                           salinity = mean(salinity, na.rm = TRUE), 
                                                                           vertical = mean(vertical, na.rm = TRUE))
    temperature[ilayer, mass_summ$object_, itime] <- mass_summ$temperature
    salinity[ilayer, mass_summ$object_, itime] <- mass_summ$salinity
    vertical[ilayer, mass_summ$object_, itime] <- mass_summ$vertical
    
    face_cells$u <- raster::extract(u0[[ilayer]], face_cells$cell_)
    face_cells$v <- raster::extract(v0[[ilayer]], face_cells$cell_)
  
    face_summ <- face_cells %>% 
      dplyr::group_by(object_) %>% 
      dplyr::summarize(u = sum(u / length, na.rm = TRUE), v = sum(v/length, na.rm = TRUE), 
                       transport = sqrt(u * u + v * v) * time_step * c(dz, tail(dz, 1))[ilayer])
    transport[ilayer, face_summ$object_, itime] <- face_summ$transport
  }
  delif(temp0)
  delif(u0)
  delif(salt0)
  delif(v0)
  rm(temp0, u0, salt0, v0)
  print(itime)
  gc()
  if (itime %% 10 == 0) {
    image(transport[,,itime], col = cols, zlim = c(0, 100), useRaster = TRUE)
   
  }
}

ncvar_put(nctran, "transport", transport)
ncvar_put(ncmass, "temperature", temperature)
ncvar_put(ncmass, "salinity", salinity)
ncvar_put(ncmass, "verticalflux", vertical)
nc_close(nctran)
nc_close(ncmass)

# 
# ## Plot
# range_raster <- function(x) {
#   c(min = min(cellStats(x, min)), max = max(cellStats(x, max)))
# }
# temp0  <- setValues(subset(hhh, 1:ncol(temp_atlantis)), temp_atlantis)
# salt0 <- setValues(subset(hhh, 1:ncol(temp_atlantis)), salt_atlantis)
# names(temp0) <- sprintf("depth_%i", as.integer(abs(atlantis_depths)))
# names(salt0) <- names(temp0)
# 
# u0  <- setValues(subset(u * NA_real_, 1:ncol(u_atlantis)), u_atlantis)
# v0  <- setValues(subset(v * NA_real_, 1:ncol(v_atlantis)), v_atlantis)
# names(u0) <- sprintf("depth_%i", as.integer(abs(atlantis_depths)))
# names(v0) <- names(temp0)
# 
# ix <- c(1, 4, 7, 10)
# library(quadmesh)
# mesh_plot(temp0[[2]], crs = projection(boxes), coords = roms_ll)
# plot(subset(temp0, ix), zlim = range_raster(temp0), col = palr::sstPal(26))
# plot(subset(salt0, ix), zlim = range_raster(salt0), col = palr::sstPal(26))
# 
# plot(subset(u0, ix), zlim = range_raster(u0), col = palr::sstPal(26))
# plot(subset(v0, ix), zlim = range_raster(v0), col = palr::sstPal(26))
# 
# get_coords <- memoise::memoise(function(f) {
#     cd <- romscoords(f, transpose = TRUE)
#     # raster::crop(cd, romsmap(geom, cd))
#     cd
#     }
# )
# plot_roms <- function(file, varname = "temp", map = NULL, slice = c(1, 1)) {
#   coords <- get_coords(file)
#   ro_map <- NULL
#   if (!is.null(map)) ro_map <- romsmap(map, coords)
#   if (varname == "[uv]") {
#     u <- romsdata(file, varname = "u", slice = slice)
#     v <- romsdata(file, varname = "v", slice = slice)
#     dat <- sqrt(u * u + v* v)
#   } else {
#     dat <- romsdata(file, varname = varname, slice = slice)
#   }
#   if (!is.null(ro_map)) {
#     dat <- raster::crop(dat, ro_map)
#     coords <- raster::crop(coords, ro_map)
#   }
#   crs <- NULL
#   if (!is.null(map)) crs <- raster::projection(map)
#   
#   quadmesh::mesh_plot(dat, crs = crs, coords = coords)
#   if (!is.null(map)) plot(map, add = TRUE)
#   invisible(NULL)
# }
# plot_roms(roms_file, varname = "[uv]", map = boxes, slice = c(31, 12))
