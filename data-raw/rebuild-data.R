## Rebuild bundled data as files on disk (not serialized R objects)
## terra SpatRaster objects can't be serialized, so we store as .tif/.gpkg
## and use active bindings in .onLoad to create live objects

library(terra)

dir.create("inst/raster", recursive = TRUE)
dir.create("inst/vector", recursive = TRUE)

## antarctica - store as GeoPackage
ant <- terra::vect(
  subset(rnaturalearth::countries110, SOVEREIGNT == "Antarctica", select = "SOV_A3")
)
terra::writeVector(ant, "inst/vector/antarctica.gpkg", overwrite = TRUE)

## ice_coords and ice_fake - need original polar stereo ice data
if (requireNamespace("raadtools", quietly = TRUE)) {
  ice_r <- raadtools::readice("2017-04-08")[[1]]
  ice_t <- terra::rast(ice_r)
  xy <- terra::xyFromCell(ice_t, seq_len(terra::ncell(ice_t)))
  ll <- reproj::reproj(xy, target = "+proj=longlat +datum=WGS84",
                       source = terra::crs(ice_t))

  lon_layer <- terra::setValues(ice_t[[1]], ll[, 1])
  lat_layer <- terra::setValues(ice_t[[1]], ll[, 2])
  ice_coords <- c(lon_layer, lat_layer)
  names(ice_coords) <- c("lon", "lat")
  terra::ext(ice_coords) <- terra::ext(0, ncol(ice_coords), 0, nrow(ice_coords))
  terra::crs(ice_coords) <- ""

  ice_fake <- ice_t
  terra::ext(ice_fake) <- terra::ext(0, ncol(ice_fake), 0, nrow(ice_fake))
  terra::crs(ice_fake) <- ""

  terra::writeRaster(ice_coords, "inst/raster/ice_coords.tif", overwrite = TRUE)
  terra::writeRaster(ice_fake, "inst/raster/ice_fake.tif", overwrite = TRUE)
} else {
  ## fallback: convert existing .rda objects
  message("raadtools not available, converting existing rda objects")
  load("data/ice_coords.rda")
  load("data/ice_fake.rda")
  ic <- terra::rast(ice_coords)
  terra::ext(ic) <- terra::ext(0, ncol(ic), 0, nrow(ic))
  terra::crs(ic) <- ""
  terra::writeRaster(ic, "inst/raster/ice_coords.tif", overwrite = TRUE)

  ifa <- terra::rast(ice_fake)
  terra::ext(ifa) <- terra::ext(0, ncol(ifa), 0, nrow(ifa))
  terra::crs(ifa) <- ""
  terra::writeRaster(ifa, "inst/raster/ice_fake.tif", overwrite = TRUE)
}

## Remove old .rda files once inst/ data is confirmed good
## file.remove("data/antarctica.rda", "data/ice_coords.rda", "data/ice_fake.rda")

message("Done. inst/raster/ and inst/vector/ populated.")
message("Remember to remove data/*.rda and LazyData from DESCRIPTION once confirmed.")
