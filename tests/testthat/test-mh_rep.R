library(testthat)
library(terra)
library(sf)
library(tidyterra)
library(stats)

testthat::test_that("mh_present runs and writes output", {
  set.seed(2458)
  n_cells <- 100 * 100
  r_clim_present <- terra::rast(ncols = 100,
                                nrows = 100,
                                nlyrs = 7)
  values(r_clim_present) <- c(
    (terra::rowFromCell(r_clim_present, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
    (terra::rowFromCell(r_clim_present, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
    (terra::colFromCell(r_clim_present, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
    (terra::colFromCell(r_clim_present, 1:n_cells) +
       (terra::rowFromCell(r_clim_present, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
    (terra::colFromCell(r_clim_present, 1:n_cells) /
       (terra::rowFromCell(r_clim_present, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
    (terra::colFromCell(r_clim_present, 1:n_cells) *
       (terra::rowFromCell(r_clim_present, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4))),
    (terra::colFromCell(r_clim_present, 1:n_cells) *
       (terra::colFromCell(r_clim_present, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4)))
  )
  names(r_clim_present) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
  terra::crs(r_clim_present) <- "EPSG:4326"
  r_clim_present_filtered <- vif_filter(r_clim_present, th = 5)
  hex_grid <- sf::st_sf(sf::st_make_grid(sf::st_as_sf(terra::as.polygons(
    terra::ext(r_clim_present)
  )), square = FALSE))
  sf::st_crs(hex_grid) <- "EPSG:4326"
  polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
  polygons$name <- c("Pol_1", "Pol_2")
  sf::st_crs(polygons) <- sf::st_crs(hex_grid)
  study_area_polygon <- sf::st_as_sf(as.polygons(terra::ext(r_clim_present)))
  sf::st_crs(study_area_polygon) <- "EPSG:4326"

  output_dir <- file.path(tempdir(), "mh_rep_test")
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
  dir.create(output_dir, showWarnings = FALSE)

  mh_rep(
    polygon = polygons,
    col_name = "name",
    climate_variables = r_clim_present_filtered,
    th = 0.95,
    dir_output = output_dir,
    save_raw = TRUE
  )

  rep_dir <- file.path(output_dir, "Representativeness")
  charts_dir <- file.path(output_dir, "Charts")
  mh_raw_dir <- file.path(output_dir, "Mh_Raw")

  expect_true(dir.exists(rep_dir), label = "Representativeness directory exists")
  expect_true(dir.exists(charts_dir), label = "Charts directory exists")
  expect_true(dir.exists(mh_raw_dir), label = "Mh_Raw directory exists")

  expect_gt(length(list.files(rep_dir)), 0, label = "Representativeness directory is not empty")
  expect_gt(length(list.files(charts_dir)), 0, label = "Charts directory is not empty")
  expect_gt(length(list.files(mh_raw_dir)), 0, label = "Mh_Raw directory is not empty")
  unlink(output_dir, recursive = TRUE)
})
