test_that("coords_points returns matrix", {
  pts <- coords_points(ice_coords)
  expect_true(is.matrix(pts))
  expect_equal(ncol(pts), 2)
  expect_equal(colnames(pts), c("x", "y"))
  expect_equal(nrow(pts), terra::ncell(ice_coords))
})
