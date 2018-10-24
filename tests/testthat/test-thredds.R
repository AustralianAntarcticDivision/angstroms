context("thredds")
skip_on_travis()
skip_on_appveyor()
ur <- "http://tds.marine.rutgers.edu/thredds/dodsC/roms/gom/g9/daily_avg"
library(raster)
h <- raster(ur, varname = "h", ncdf = TRUE)
hr <- romsdata(ur, varname = "h")
#h0 <- raster(ur, varname = "h")
test_that("Thredds server works via ncdf4", {
  ## should get raster/ncdf4 not rgdal::readGDAL
  ## https://github.com/hypertidy/angstroms/issues/6
  expect_that(dim(h), equals(c(100, 160, 1)))
  expect_that(dim(hr), equals(dim(h)) )
})

# test_that("Thredds server works via rgdal", {
#   skip_if_not(requireNamespace("rgdal", quietly = TRUE))
#   ## this is what readGDAL sees
#   expect_that(dim(h0), equals(c(4, 2, 1)))
# })
