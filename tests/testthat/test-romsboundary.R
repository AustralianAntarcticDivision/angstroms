test_that("boundary returns wk geometry", {
  b <- romsboundary(ice_coords)
  expect_s3_class(b, "wk_wkt")
})
