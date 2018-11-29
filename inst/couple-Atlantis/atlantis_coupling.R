##devtools::install_github(c("AustralianAntarcticDivision/angstroms"))
library(angstroms)
library(rbgm) ## read BGM
library(bgmfiles) ## archive of BGM files
library(raadtools)
library(tidync)
library(dplyr)

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
faces <- faceSpatial(bgm)
romsbox <- romsmap(boxes, roms_ll)
roms_ll <- crop(roms_ll , romsbox, snap = "out")

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
atlantis_depths <- -cumsum(c(0, rev(rbgm::build_dz(-deepest_depth))))

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


## face properties

u <- crop(romsdata3d(roms_file, varname = "u", slice = 1), romsbox)
v <- crop(romsdata3d(roms_file, varname = "v", slice = 1), romsbox)
u_atlantis <- v_atlantis <- face_index * NA_real_
for (i in seq_len(nrow(u_atlantis))) {
  u_atlantis[i, ] <- u@data@values[i, face_index[i, ]]
  v_atlantis[i, ] <- v@data@values[i, face_index[i, ]]
  
  if (i %% 1000 == 0) print(i)
}



## Plot
range_raster <- function(x) {
  c(min = min(cellStats(x, min)), max = max(cellStats(x, max)))
}
temp0  <- setValues(subset(hhh, 1:ncol(temp_atlantis)), temp_atlantis)
salt0 <- setValues(subset(hhh, 1:ncol(temp_atlantis)), salt_atlantis)
names(temp0) <- sprintf("depth_%i", as.integer(abs(atlantis_depths)))
names(salt0) <- names(temp0)

u0  <- setValues(subset(u * NA_real_, 1:ncol(u_atlantis)), u_atlantis)
v0  <- setValues(subset(v * NA_real_, 1:ncol(v_atlantis)), v_atlantis)
names(u0) <- sprintf("depth_%i", as.integer(abs(atlantis_depths)))
names(v0) <- names(temp0)

ix <- c(1, 4, 7, 10)
library(quadmesh)
mesh_plot(temp0[[2]], crs = projection(boxes), coords = roms_ll)
plot(subset(temp0, ix), zlim = range_raster(temp0), col = palr::sstPal(26))
plot(subset(salt0, ix), zlim = range_raster(salt0), col = palr::sstPal(26))

plot(subset(u0, ix), zlim = range_raster(u0), col = palr::sstPal(26))
plot(subset(v0, ix), zlim = range_raster(v0), col = palr::sstPal(26))

get_coords <- memoise::memoise(function(f) {
    cd <- romscoords(f, transpose = TRUE)
    # raster::crop(cd, romsmap(geom, cd))
    cd
    }
)
plot_roms <- function(file, varname = "temp", map = NULL, slice = c(1, 1)) {
  coords <- get_coords(file, map)
  ro_map <- NULL
  if (!is.null(map)) ro_map <- romsmap(map, coords)
  if (varname == "[uv]") {
    u <- romsdata(file, varname = "u", slice = slice)
    v <- romsdata(file, varname = "v", slice = slice)
    dat <- sqrt(u * u + v* v)
  } else {
    dat <- romsdata(file, varname = varname, slice = slice)
  }
  if (!is.null(ro_map)) {
    dat <- raster::crop(dat, ro_map)
    coords <- raster::crop(coords, ro_map)
  }
  crs <- NULL
  if (!is.null(map)) crs <- raster::projection(map)
  
  quadmesh::mesh_plot(dat, crs = crs, coords = coords)
  if (!is.null(map)) plot(map, add = TRUE)
  invisible(NULL)
}
plot_roms(roms_file, varname = "[uv]", map = boxes, slice = c(31, 12))
