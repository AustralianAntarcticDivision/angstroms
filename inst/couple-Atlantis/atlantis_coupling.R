
## crop to region with cell safety 
## really characterize different levels by box, and overall


#' Convenience function to transform map projection . 
#'
#' Transform `x` to whatever the projection of `to` is. 
#' @param x  object to transform
#' @param to object with a map projection
#'
#' @return `x`, transformed
#' @export
#' @importFrom raster projection
#' @importFrom sp spTransform
project_to <- function(x, to) {
  spTransform(x, CRS(projection(to)))
}

## which box does each point fall in
index_box <- function(box_sp, roms_ll) {
  ind <- sp::over(project_to(coords_points(roms_ll), box_sp) , as(box_sp, "SpatialPolygons"))
  #  tibble(box = ind,
  tibble(box = box_sp$label[ind], 
         cell = seq_len(ncell(roms_ll))) %>% 
    filter(!is.na(box))
}

## 
index_face <- function(face_sp, roms_ll) {
  roms_face <- romsmap(project_to(face_sp, "+init=epsg:4326"), roms_ll)
  ind_face <- cellFromLine(roms_ll, roms_face)
  tibble(face = roms_face$label[rep(seq_len(nrow(roms_face)), lengths(ind_face))], 
         cell = unlist(ind_face))
}



##devtools::install_github(c("mdsumner/angstroms"))
library(angstroms)
library(rbgm) ## read BGM
library(bgmfiles) ## archive of BGM files
library(raadtools)
library(tidync)
library(dplyr)

cpolar <- raadtools:::cpolarfiles()
  
## for each file we know about each separate time, and this is lvar = 4, level = [band_level]
## there are 31 depth levels (surface water is the first one)
## files have 28, 30, or 31 time steps file_db %>% group_by(fullname) %>% tally() %>% distinct(n)
## so we can't just assume a constant
file_db <- purrr::map_df(cpolar$fullname[2], 
                         function(x) tibble::tibble(fullname = x,
                                                    time = tidync(x)$transforms$ocean_time$ocean_time) %>% 
                           dplyr::mutate(band_level = row_number()))


## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- rbgm:::read_bgm(bfile)


library(dplyr)


clamp <- function(x) {
  x[x < 1] <- 1
  x[x > 31] <- 31
  x
}
roms_file <- file_db$fullname[1]
roms_ll <- romscoords(roms_file, transpose = TRUE)
boxes <- boxSpatial(bgm)
faces <- faceSpatial(bgm)
romsbox <- romsmap(boxes, roms_ll)
roms_ll <- crop(roms_ll , romsbox, snap = "out")

h <- crop(romsdata(roms_file, varname = "h"), romsbox, snap = "out")
hhh <- crop(romsdepth(roms_file), romsbox, snap = "out")

## Cs_r is the S-coord stretching
Cs_r <- rawdata(roms_file, "Cs_r")
#list_nc_z_index <- vector('list', nrow(box_roms_index))
deepest_depth <- max(raster::extract(h, tabularaster::cellnumbers(h, sf::st_as_sf(romsbox))$cell_), na.rm = TRUE)
## here we must use the same as we used for the nc output files
atlantis_depths <- -cumsum(c(0, rev(rbgm::build_dz(-deepest_depth))))


box_mass_index <- matrix(NA_integer_, nrow(hhh@data@values), length(atlantis_depths))
for (i in seq_len(nrow(box_mass_index))) {
  idx <- clamp(findInterval(atlantis_depths, hhh@data@values[i, ]))
  ## max out at ~200m beyond the deepest depth
  idx[(atlantis_depths - min(hhh@data@values[i, ])) < -20] <- NA_integer_
  box_mass_index[i, ] <- idx
  if (i %% 1000 == 0) print(i)
}
temp <- crop(romsdata3d(roms_file, varname = "temp", slice = 1), romsbox)
salt <- crop(romsdata3d(roms_file, varname = "salt", slice = 1), romsbox)
salt_atlantis <- temp_atlantis <- box_mass_index * NA_real_
for (i in seq_len(nrow(temp_atlantis))) {
  temp_atlantis[i, ] <- temp@data@values[i, box_mass_index[i, ]]
  salt_atlantis[i, ] <- salt@data@values[i, box_mass_index[i, ]]
  
  if (i %% 1000 == 0) print(i)
}
range_raster <- function(x) {
  c(min = min(cellStats(x, min)), max = max(cellStats(x, max)))
}
temp0  <- setValues(subset(hhh, 1:ncol(temp_atlantis)), temp_atlantis)
salt0 <- setValues(subset(hhh, 1:ncol(temp_atlantis)), salt_atlantis)
ix <- c(1, 4, 7, 10)
names(temp0) <- sprintf("depth_%i", as.integer(abs(atlantis_depths)))
names(salt0) <- names(temp0)
plot(subset(temp0, ix), zlim = range_raster(temp0), col = palr::sstPal(26))
plot(subset(salt0, ix), zlim = range_raster(salt0), col = palr::sstPal(26))


