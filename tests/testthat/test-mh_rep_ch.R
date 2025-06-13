library(testthat)
library(terra)
library(sf)
library(stats)

testthat::test_that("mh_rep_ch works correctly", {
  set.seed(2458)
  n_cells <- 100 * 100
  r_clim_present <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
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
  r_clim_present_filtered <- r_clim_present[[1:5]]
  r_clim_future <- r_clim_present_filtered + 2
  names(r_clim_future) <- names(r_clim_present_filtered)
  hex_grid <- sf::st_sf(sf::st_make_grid(sf::st_as_sf(terra::as.polygons(
    terra::ext(r_clim_present))), square = FALSE))
  sf::st_crs(hex_grid) <- "EPSG:4326"
  polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
  polygons$name <- c("Pol_1", "Pol_2")
  study_area_polygon <- sf::st_as_sf(terra::as.polygons(terra::ext(r_clim_present)))
  sf::st_crs(study_area_polygon) <- "EPSG:4326"
  terra::plot(r_clim_present[[1]])
  terra::plot(polygons,
              add = TRUE,
              color = "transparent",
              lwd = 3)
  terra::plot(
    study_area_polygon,
    add = TRUE,
    col = "transparent",
    lwd = 3,
    border = "red")

  output_dir <- file.path(tempdir(), "ClimaRepOutput_mh_rep_ch_test")
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
  dir.create(output_dir, showWarnings = FALSE)

  output_mh_rep_ch <- mh_rep_ch(
    polygon = polygons,
    col_name = "name",
    present_climate_variables = r_clim_present_filtered,
    future_climate_variables = r_clim_future,
    study_area = study_area_polygon,
    th = 0.95,
    model = "ExampleModel",
    year = "2070",
    dir_output = output_dir,
    save_raw = TRUE)

  change_dir <- file.path(output_dir, "Change")
  charts_dir <- file.path(output_dir, "Charts")
  mh_raw_pre_dir <- file.path(output_dir, "Mh_Raw_Pre")
  mh_raw_fut_dir <- file.path(output_dir, "Mh_Raw_Fut")

  expect_true(dir.exists(change_dir), label = "Change directory exists")
  expect_true(dir.exists(charts_dir), label = "Charts directory exists")
  expect_true(dir.exists(mh_raw_pre_dir), label = "Mh_Raw_Pre directory exists")
  expect_true(dir.exists(mh_raw_fut_dir), label = "Mh_Raw_Fut directory exists")

  # Comprobar que las carpetas no están vacías
  expect_gt(length(list.files(change_dir)), 0, label = "Change directory is not empty")
  expect_gt(length(list.files(charts_dir)), 0, label = "Charts directory is not empty")
  expect_gt(length(list.files(mh_raw_pre_dir)), 0, label = "Mh_Raw_Pre directory is not empty")
  expect_gt(length(list.files(mh_raw_fut_dir)), 0, label = "Mh_Raw_Fut directory is not empty")
  unlink(output_dir, recursive = TRUE)
})
