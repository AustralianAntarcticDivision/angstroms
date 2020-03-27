test_that("sp build works", {
  expect_s4_class(
     romsmap(antarctica, ice_coords), "SpatialPolygonsDataFrame")
})
