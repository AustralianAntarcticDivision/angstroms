library(raadtools)
ice <- readice("2017-04-08")[[1]]
xy <- rgdal::project(coordinates(ice), projection(ice), inv = TRUE)
ice_coords <- setExtent(brick(setValues(ice, xy[, 1]), setValues(ice, xy[, 2])), extent(0, ncol(ice), 0, nrow(ice)))
ice_fake <- setExtent(ice, extent(0, ncol(ice), 0, nrow(ice)))
projection(ice_coords) <- NA
projection(ice_fake) <- NA

devtools::use_data(ice_coords)
devtools::use_data(ice_fake)