library(testthat)
library(terra)

set.seed(235)
n_cells <- 20 * 20
r_clim <- rast(ncols = 20, nrows = 20, nlyrs = 3)
values(r_clim) <- cbind(
  1:n_cells * 0.1 + rnorm(n_cells, 0, 2),
  1:n_cells * 0.05 + rnorm(n_cells, 0, 1),
  rnorm(n_cells, 10, 3)
)
names(r_clim) <- c("varA", "varB", "varC")
terra::crs(r_clim) <- "EPSG:4326"


# Escribir el test
testthat::test_that("vif_filter works correctly", {
  r_filtered <- vif_filter(x = r_clim, th = 5)
  expect_s4_class(r_filtered, "SpatRaster")
  expect_true(terra::nlyr(r_filtered) <= terra::nlyr(r_clim))
  expect_equal(terra::crs(r_filtered), terra::crs(r_clim))
  expect_equal(dim(r_filtered)[1:2], dim(r_clim)[1:2])

})
