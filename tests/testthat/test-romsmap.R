test_that("romsmap works with sf input", {
  skip_if_not_installed("sf")
  idx <- romsmap(antarctica, ice_coords)
  ## returns a wk geometry in index space
  expect_true(inherits(idx, "wk_wkb") || inherits(idx, "wk_wkt") || is.data.frame(idx))
})

test_that("romsmap works with matrix input", {
  pts <- matrix(c(147, -42, 150, -45), ncol = 2, byrow = TRUE)
  idx <- romsmap(pts, ice_coords)
  expect_true(is.data.frame(idx))
  expect_equal(ncol(idx), 2)
  expect_equal(nrow(idx), 2)
})
