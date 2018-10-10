## ROMS Atlantis coupling MDSumner, Jessica Melbourne-Thomas April 2017

## NOTES

#https://confluence.csiro.au/display/Atlantis/Hydro+FAQ

#https://github.com/cecilieha/NoBA/blob/master/create_flux_ncdf_atlantis.R

##corners <- read.delim("https://raw.githubusercontent.com/cecilieha/NoBA/master/corners_neighbours_nordic.txt")

## curvlinear vector rotation: https://www.myroms.org/forum/viewtopic.php?f=3&t=295

## hydro stuff
# https://confluence.csiro.au/display/Atlantis/Forcing+Files
# https://confluence.csiro.au/display/Atlantis/Current+Forcing+File+Structure
# https://confluence.csiro.au/display/Atlantis/HydroConstruct
# https://confluence.csiro.au/display/Atlantis/Matlab+-+Hydro

## dummy hydro
# sohydro <- "~/projects/AtlantisEastAntartica/SOhydrodummy.nc"
# library(ncdump)
# nc <- NetCDF(sohydro)

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

#' return the ramp of positive depths from surface down 
#' (so that the order is native to the NetCDF order)
roms_level <- function(Cs_r, h, cell) {
  as.vector(extract(h, cell)) * Cs_r
}

roms_z <- function(Cs_r, h, cell) {
  out <- flip(raster(matrix(rep(extract(h, cell), each = length(Cs_r)) *  rep(Cs_r, length(cell)), 
                            length(Cs_r))), "y")
  setExtent(out, extent(0, ncol(out), 0, nrow(out)))
}

## new sf approach for rbgm
# make_boxes <- function(tab) {
#   tab <- tab[grepl("box", tab$vert), ]
#   lapply(split(tab[, c("X", "Y")], tab$vert), function(x) sf::st_polygon(list(as.matrix(x))))
# }
# make_faces <- function(tab) {
#   lapply(split(tab[, c("x1", "y1", "x2", "y2")], tab$face), function(x) sf::st_linestring(matrix(unlist(x), ncol = 2, byrow = TRUE)))
# }

##devtools::install_github(c("mdsumner/angstroms"))
library(angstroms)
library(rbgm) ## read BGM
library(bgmfiles) ## archive of BGM files
library(raadtools)
library(tidync)
library(dplyr)

cpolar <- raadtools:::cpolarfiles()
  
file_db <- purrr::map_df(cpolar$fullname, 
                            function(x) tibble::tibble(fullname = x,
                                                       time = tidync(x)$transforms$ocean_time$ocean_time) %>% 
                              dplyr::mutate(band_level = row_number()))

## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- rbgm:::read_bgm(bfile)


## get the longitude/latitude arrays
roms_ll <- romscoords(file_db$fullname[1], transpose = TRUE)

## we need the unsullied boxes to identify points inside them
# boxes <- as(sf::st_sf(label = bgm$box$box, 
#                       geometry =  sf::st_sfc(make_boxes(bgm$vertices), 
#                                              crs = bgm$meta$projection), stringsAsFactors = FALSE), "Spatial")
boxes <- boxSpatial(bgm)
faces <- faceSpatial(bgm) 
#roms_face <- romsmap(project_to(faceSpatial(bgm), "+init=epsg:4326"), roms_ll)
roms_face <- romsmap(project_to(faces, "+init=epsg:4326"), roms_ll)




## build the index for each box to the ROMS cells it contains
## and each face for the ROMS cells it traverses
library(dplyr)
box_roms_index <- index_box(boxes, roms_ll)
ind_face <- cellFromLine(romsdata(file_db$fullname[1], "u"), roms_face)
face_roms_index <- tibble(face = roms_face$label[rep(seq_len(nrow(roms_face)), lengths(ind_face))], 
                          cell = unlist(ind_face))


## important to readAll here, else extract is very slow in the loop
h <- brick(readAll(raster(file_db$fullname[1], varname = "h")))
## Cs_r is the S-coord stretching
Cs_r <- rawdata(file_db$fullname[1], "Cs_r")


## build the level index between Atlantis and ROMS

list_nc_z_index <- vector('list', nrow(box_roms_index))
max_depth <- max(extract(h, unique(box_roms_index$cell)))
atlantis_depths <- -rev(cumsum(c(0, rev(rbgm::build_dz(-max_depth)))))
for (i in seq_len(nrow(box_roms_index))) {
  rl <- roms_level(Cs_r, h, box_roms_index$cell[i])
  ## implicit 0 at the surface, and implicit bottom based on ROMS
  list_nc_z_index[[i]] <- length(atlantis_depths) - findInterval(rl, atlantis_depths) + 1
if (i %% 1000 == 0) print(i)
}


# 
# ll <- extract(roms_ll, box_roms_index$cell)
# box_roms_index$lon <- ll[, 1]
# box_roms_index$lat <- ll[, 2]
# library(ggplot2)
# ggplot(box_roms_index, aes(lon, lat, colour = box)) + geom_point(pch = ".")
# 

## join the box-xy-index to the level index
box_z_index <- bind_rows(lapply(list_nc_z_index, 
                 function(x) tibble(atlantis_level = pmax(1, x), roms_level = seq_along(x))), 
          .id = "cell_index") %>% 
  inner_join(mutate(box_roms_index, cell_index = as.character(row_number()))) %>% 
  select(-cell_index)


# ll <- extract(roms_ll, box_z_index$cell)
# box_z_index$lon <- ll[, 1]
# box_z_index$lat <- ll[, 2]
# library(ggplot2)
# ggplot(box_z_index, aes(lon, lat, colour = roms_level)) + geom_point(pch = ".") + 
#   facet_wrap(~atlantis_level)


for (i in seq_len(nrow(face_roms_index))) {
  rl <- roms_level(Cs_r, h, face_roms_index$cell[i])
  ## implicit 0 at the surface, and implicit bottom based on ROMS
  list_nc_z_index[[i]] <- length(atlantis_depths) -findInterval(rl, atlantis_depths) + 1
  if (i %% 1000 == 0) print(i)
}
## join the face-xy-index to the level index
face_z_index <- bind_rows(lapply(list_nc_z_index, 
                                function(x) tibble(atlantis_level = x, roms_level = seq_along(x))), 
                         .id = "cell_index") %>% 
  inner_join(mutate(face_roms_index, cell_index = as.character(row_number()))) %>% 
  select(-cell_index)



#z_index$atlantis_depth <- c(0, rev(rbgm::build_dz(box_roms_index$botz[i])))[z_index$atlantis_level]
#z_index$roms_depth <- rl[z_index$roms_level]
## driver function to loop over x raster by levels
## matching $cell for the right group
extract_at_level <- function(x, cell_level) {
  ulevel <- unique(cell_level$level)
  values <- numeric(nrow(cell_level))
  for (ul in seq_along(ulevel)) {
    asub <- cell_level$level == ulevel[ul]
    values[asub] <- extract(x[[ulevel[ul]]], 
            cell_level$cell[asub])
  }
  values
}

#file_db <- file_db[1, ]
box_props <- face_props <- vector("list", nrow(file_db))
i_timeslice <- 1
set_indextent <- tabularaster:::set_indextent
for (i_timeslice in seq(nrow(file_db))) {
  roms_file <- file_db$fullname[i_timeslice]
  level <- file_db$band_level[i_timeslice]
  ru <- set_indextent(brick(roms_file, varname = "u", lvar = 4, level = level))
  face_z_index$u <- extract_at_level(ru, rename(face_z_index, level = roms_level, cell = cell))
  rv <- set_indextent(brick(roms_file, varname = "v", lvar = 4, level = level))
  face_z_index$v <- extract_at_level(rv, rename(face_z_index, level = roms_level, cell = cell))
  rdata <- set_indextent(brick(roms_file, varname = "temp", lvar = 4, level = level))
  box_z_index$temp <- extract_at_level(rdata, rename(box_z_index, level = roms_level, cell = cell))
  rdata <- set_indextent(brick(roms_file, varname = "salt", lvar = 4, level = level))
  box_z_index$salt <- extract_at_level(rdata, rename(box_z_index, level = roms_level, cell = cell))
  box_props[[i_timeslice]] <- box_z_index %>% group_by(atlantis_level, box) %>% 
    summarize(temp = mean(temp, na.rm = TRUE), salt = mean(salt ,na.rm = TRUE)) %>% 
    mutate(band_level = level)
  face_props[[i_timeslice]] <-  face_z_index %>% group_by(atlantis_level, face) %>% 
    summarize(flux = mean(sqrt(u * u + v * v), na.rm = TRUE) * 24 * 3600) %>% 
    mutate(band_level = level)
}

box_props <- bind_rows(box_props)
face_props <- bind_rows(face_props)
library(ggplot2)
ggplot(box_props, aes(x = temp, y = salt, colour = atlantis_level)) + geom_point() + facet_wrap(~box)
ggplot(face_props, aes(x = atlantis_level, y = flux, colour = atlantis_level)) + geom_point() 

as_tib_slab <- function(trans, slab, activ) {
  tib <- list()
  tib[[activ]] <- as.vector(slab)
  tib <- as_tibble(tib)
  prod_dims <- 1
  total_prod <- prod(dim(slab))
  
  for (i in seq_along(trans)) {
    nm <- names(trans)[i]
    nr <- nrow(trans[[i]])
    tib[[nm]] <- rep(trans[[nm]][[nm]], each = prod_dims, length.out = total_prod)
    prod_dims <- prod_dims * nr
  }
  tib
}
library(ncdump)
library(dplyr)
ncfile <- "~/Git/shinyrAtlantis/egfiles/SOtempdummy.nc"
con <- ncdf4::nc_open(ncfile)
trans <- ncdump::NetCDF(ncfile) %>% activate("temperature") %>% filtrate() 
hslab <- bind_rows(lapply(trans, function(x) tibble(name = x$name[1], start = min(x$step), count = length(x$step))))
slab <- ncdf4::ncvar_get(con, "temperature", 
                         start = hslab$start, 
                         count = hslab$count)
temp <- as_tib_slab(trans, slab, activ = "temperature")

ggplot(box_props %>% filter(box %in% c("Box1", "Box2", "Box3", "Box4")) %>% 
         group_by(box, atlantis_level, band_level) %>% 
         summarize(temp  = mean(temp)), aes(x = band_level, y = temp, colour = factor(atlantis_level))) + 
  geom_point() + facet_wrap(~box)



# 
# 
# plot(romsmap(boxes[15, ], roms_ll), add = F, border = "red", lwd = 4)
# plot(sqrt(ru[[10]]^2 + rv[[10]]^2), add = TRUE, 
#      col = colorRampPalette(c("#132B43", high = "#56B1F7"))(10), zlim = c(0, 0.1))
# plot(romsmap(boxes[15, ], roms_ll), add = TRUE, border = "red", lwd = 4)
# ggplot(mutate(z_index %>% 
#                 filter(box == 15), 
#               x = xFromCell(ru, cell), y = yFromCell(ru, cell)), 
#        aes(x = x, y = y, fill = sqrt(u * u + v * v))) + geom_tile() + coord_equal() + 
#   facet_wrap(~roms_level)
# 
# 
# rm(ru, rv)
# 
# 
# ggplot(z_index %>% mutate(x = xFromCell(rdata, cell), 
#                           y = yFromCell(rdata, cell)), 
#        aes(x = x, y = y, fill = temp)) + geom_tile() + 
#   coord_equal()
# 
# ggplot(mutate(z_index %>% 
#                 filter(roms_level %% 3 == 0), 
#               x = xFromCell(ru, cell), y = yFromCell(ru, cell)), 
#        aes(x = x, y = y, fill = temp)) + geom_tile() + coord_equal() + 
#   facet_wrap(~roms_level)
