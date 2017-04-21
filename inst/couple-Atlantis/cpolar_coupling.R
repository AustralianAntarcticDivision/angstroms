## ROMS Atlantis coupling MDSumner, Jessica Melbourne-Thomas April 2017

## NOTES

#https://confluence.csiro.au/display/Atlantis/Hydro+FAQ

#https://github.com/cecilieha/NoBA/blob/master/create_flux_ncdf_atlantis.R

##corners <- read.delim("https://raw.githubusercontent.com/cecilieha/NoBA/master/corners_neighbours_nordic.txt")

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
## devtools::install_github(c("mdsumner/angstroms"))
library(angstroms)
library(rbgm) ## read BGM
library(bgmfiles) ## archive of BGM files
library(raadtools)
library(ncdump)
## get the Circumpolar files
cpolar <- raadtools:::cpolarfiles()
roms_file <- cpolar$fullname[2]

## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- bgmfile(bfile)


## get the longitude/latitude arrays
roms_ll <- romscoords(roms_file, transpose = TRUE)
## get a representative data array (slice is depth, then time, 1 is surface)
depth_time <- c(1, 1)
## note we don't bother matching u/v grids correctly here, yet
roms_grid <- sqrt(romsdata(roms_file, "u", transpose = TRUE, slice = depth_time) ^2 + 
                  romsdata(roms_file, "v", transpose = TRUE, slice = depth_time) ^2)


## map the BGM boxes and faces to the ROMS grid
roms_box <- romsmap(project_to(boxSpatial(bgm), "+init=epsg:4326"), roms_ll)
roms_face <- romsmap(project_to(faceSpatial(bgm), "+init=epsg:4326"), roms_ll)

## plot to check
plot(roms_grid, asp = "", col = viridis::viridis(100), addfun = function() plot(roms_box, add = TRUE, border = "white"), nr = 2)

## we need the unsullied boxes to identify points inside them
boxes <- boxSpatial(bgm)
## which box does each point fall in
ind_box <- over(project_to(coords_points(roms_ll), boxes) , as(boxes, "SpatialPolygons"))

## build the index for each box to the ROMS cells it contains
## and each face for the ROMS cells it traverses
library(dplyr)
box_roms_index <- tibble(box = ind_box, cell = seq_len(ncell(roms_ll))) %>% 
  filter(!is.na(box))

ind_face <- cellFromLine(roms_grid, roms_face)
face_roms_index <- tibble(face = rep(seq_len(nrow(roms_face)), lengths(ind_face)), 
                          cell = unlist(ind_face))


angstroms:::rawdata(roms_file, "Cs_r")
angstroms::romsdata(roms_file, "")

## test 
l <- vector("list", 31 * 31)
i <- 1
box_data <- fortify(roms_box) %>% 
  mutate(.id = as.integer(.id) - 1) %>%  
  inner_join(roms_box@data, c(".id"= "box_id"))
face_data <- fortify(roms_face) %>% 
  mutate(.id = as.integer(.id) - 1) %>%  
  inner_join(roms_face@data, c(".id"= ".fx0"))

system.time({
  for(ilevel in 1:31) {
    for (itime in 1:31) {
      i <- i + 1
      depth_time <- c(ilevel, itime)
      box_roms_index$u <- extract(romsdata(roms_file, "u", slice = depth_time), box_roms_index$cell)
      box_roms_index$v <- extract(romsdata(roms_file, "v", slice = depth_time), box_roms_index$cell)
      face_roms_index$u <- extract(romsdata(roms_file, "u", slice = depth_time), face_roms_index$cell)
      face_roms_index$v <- extract(romsdata(roms_file, "v", slice = depth_time), face_roms_index$cell)
      
      
      box_uv_summary <- box_roms_index %>% group_by(box) %>% summarize(u = mean(u), v = mean(v))
      roms_box$mag <- sqrt(box_uv_summary$u^2 + box_uv_summary$v^2)
      library(ggplot2)
      
      face_uv_summary <- face_roms_index %>% group_by(face) %>% summarize(u = mean(u), v = mean(v))
      roms_face$mag <- sqrt(face_uv_summary$u^2 + face_uv_summary$v^2)
      
     l[[i]] <- ggplot(face_data, aes(x = long, y = lat, group = group, colour = mag)) + 
        
        geom_path( 
                  lwd=1.05, lty = 1)
      
      
      
      
    }
    print(ilevel)
  }
})

box_uv_summary <- box_roms_index %>% group_by(box) %>% summarize(u = mean(u), v = mean(v))
roms_box$mag <- sqrt(box_uv_summary$u^2 + box_uv_summary$v^2)
library(ggplot2)

face_uv_summary <- face_roms_index %>% group_by(face) %>% summarize(u = mean(u), v = mean(v))
roms_face$mag <- sqrt(face_uv_summary$u^2 + face_uv_summary$v^2)

box_data <- fortify(roms_box) %>% 
  mutate(.id = as.integer(.id) - 1) %>%  
  inner_join(roms_box@data, c(".id"= "box_id"))
face_data <- fortify(roms_face) %>% 
  mutate(.id = as.integer(.id) - 1) %>%  
  inner_join(roms_face@data, c(".id"= ".fx0"))
ggplot(box_data, aes(x = long, y = lat, group = group, fill = mag)) + 
  geom_polygon()  +
  geom_path(data = face_data %>% rename(facemag = mag), 
            aes(x = long, y = lat, group = group, colour = facemag, fill = NULL), 
            lwd=1.05, lty = 2)


