library(testthat)
library(terra)
library(sf)


# --- Escribir el Test ---
testthat::test_that("mh_present_future runs and writes output", {

  set.seed(2458)
  n_cells <- 100 * 100
  r_clim_present <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
  values(r_clim_present) <- c((rowFromCell(r_clim_present, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
                              (rowFromCell(r_clim_present, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
                              (colFromCell(r_clim_present, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
                              (colFromCell(r_clim_present, 1:n_cells) + (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                              (colFromCell(r_clim_present, 1:n_cells) / (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                              (colFromCell(r_clim_present, 1:n_cells) * (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                              (colFromCell(r_clim_present, 1:n_cells) * (colFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)))
  names(r_clim_present) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
  terra::crs(r_clim_present) <- "EPSG:4326"
  terra::plot(r_clim_present)
  r_clim_present_filtered <- ClimaRep::vif_filter(r_clim_present, th = 5)
  hex_grid <- sf::st_sf(
    sf::st_make_grid(
      sf::st_as_sf(
        terra::as.polygons(
          terra::ext(r_clim_present))),
      square = FALSE))
  sf::st_crs(hex_grid) <- "EPSG:4326"
  polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
  polygons$name <- c("Pol_1", "Pol_2")
  sf::st_crs(polygons) <- sf::st_crs(hex_grid)
  study_area_polygon <- sf::st_as_sf(as.polygons(terra::ext(r_clim_present)))
  sf::st_crs(study_area_polygon) <- "EPSG:4326"
  terra::plot(r_clim_present[[1]])
  terra::plot(polygons, add = TRUE, color= "transparent", lwd = 3)
  terra::plot(study_area_polygon, add = TRUE, col = "transparent", lwd = 3, border = "red")
  r_clim_future <- r_clim_present_filtered + 2
  names(r_clim_future) <- names(r_clim_present_filtered)
  terra::crs(r_clim_future) <- terra::crs(r_clim_present_filtered)
  terra::plot(r_clim_future)
  ClimaRep::mh_rep_ch(
  polygon = polygons,
  col_name = "name",
  present_climatic_variables = r_clim_present_filtered,
  future_climatic_variables = r_clim_future,
  study_area = study_area_polygon,
  th = 0.95,
  model = "MODEL",
  year = "2070",
  dir_output = tempdir(),
  save_intermediate_raster = TRUE)


  expect_true(dir.exists(temp_dir_out_future))
  output_files <- list.files(temp_dir_out_future, recursive = TRUE)
  expect_true(length(output_files) > 0, info = "No files were written to the output directory.")
})
testthat::teardown({
  if (dir.exists(temp_dir_out_future)) {
    unlink(temp_dir_out_future, recursive = TRUE)
  }
})
