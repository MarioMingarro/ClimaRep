library(testthat)
library(terra)
library(sf)
library(stats)

set.seed(2458)
n_cells <- 100 * 100
r_clim <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
values(r_clim) <- c((rowFromCell(r_clim, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
                    (rowFromCell(r_clim, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
                    (colFromCell(r_clim, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
                    (colFromCell(r_clim, 1:n_cells) + (rowFromCell(r_clim, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                    (colFromCell(r_clim, 1:n_cells) / (rowFromCell(r_clim, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                    (colFromCell(r_clim, 1:n_cells) * (rowFromCell(r_clim, 1:n_cells))+ 0.1 + rnorm(n_cells, 0, 4)),
                    (colFromCell(r_clim, 1:n_cells) * (colFromCell(r_clim, 1:n_cells))+ 0.1 + rnorm(n_cells, 0, 4)))
names(r_clim) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
terra::crs(r_clim) <- "EPSG:4326"
testthat::test_that("vif_filter works correctly", {
  vif_result <- vif_filter(r_clim, th = 5)
  r_clim_filtered <- vif_result$filtered_raster
  expect_s4_class(r_clim_filtered, "SpatRaster")
  expect_true(terra::nlyr(r_clim_filtered) <= terra::nlyr(r_clim))
  expect_equal(terra::crs(r_clim_filtered), terra::crs(r_clim))
  expect_equal(dim(r_clim_filtered)[1:2], dim(r_clim)[1:2])
})
