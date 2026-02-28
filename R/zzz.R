# Active bindings for bundled data
# terra SpatRaster objects are C++ references that can't survive
# R serialization, so we store as .tif and bind on load.
# See https://stackoverflow.com/a/70722230

.onLoad <- function(libname, pkgname) {
  rasters <- c("ice_coords", "ice_fake")
  for (nm in rasters) {
    local({
      name <- nm
      f <- function() {
        r <- terra::rast(
          system.file(file.path("raster", paste0(name, ".tif")),
                      package = "angstroms", mustWork = TRUE)
        )
        terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
        terra::crs(r) <- ""
        r
      }
      makeActiveBinding(name, f, topenv())
    })
  }

  ## antarctica as SpatVector from GeoPackage
  makeActiveBinding("antarctica", function() {
    terra::vect(
      system.file("vector/antarctica.gpkg",
                   package = "angstroms", mustWork = TRUE)
    )
  }, topenv())
}
