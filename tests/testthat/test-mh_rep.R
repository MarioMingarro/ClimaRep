library(testthat)
library(terra)
library(sf)
library(tidyterra)
library(stats)

set.seed(2458)
n_cells <- 100 * 100
r_clim_present <- terra::rast(ncols = 100,
                              nrows = 100,
                              nlyrs = 7)
values(r_clim_present) <- c((rowFromCell(r_clim_present, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
                            (rowFromCell(r_clim_present, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
                            (colFromCell(r_clim_present, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
                            (colFromCell(r_clim_present, 1:n_cells) + (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                            (colFromCell(r_clim_present, 1:n_cells) / (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                            (colFromCell(r_clim_present, 1:n_cells) * (rowFromCell(r_clim_present, 1:n_cells))+ 0.1 + rnorm(n_cells, 0, 4)),
                            (colFromCell(r_clim_present, 1:n_cells) * (colFromCell(r_clim_present, 1:n_cells))+ 0.1 + rnorm(n_cells, 0, 4)))
names(r_clim_present) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
terra::crs(r_clim_present) <- "EPSG:4326"
r_clim_present_filtered <- vif_filter(r_clim_present, th = 5)
hex_grid <- sf::st_sf(
  sf::st_make_grid(
    sf::st_as_sf(
      terra::as.polygons(
        terra::ext(r_clim_present)
      )), square = FALSE))
sf::st_crs(hex_grid) <- "EPSG:4326"
polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
polygons$name <- c("Pol_1", "Pol_2")
sf::st_crs(polygons) <- sf::st_crs(hex_grid)
study_area_polygon <- sf::st_as_sf(as.polygons(terra::ext(r_clim_present)))
sf::st_crs(study_area_polygon) <- "EPSG:4326"
testthat::test_that("mh_present runs and writes output", {
  mh_rep(
    polygon = polygons,
    col_name = "name",
    climatic_variables = r_clim_present_filtered,
    th = 0.95,
    dir_output = file.path(tempdir(), "ClimaRep"),
    save_raw = TRUE
  )
  dir_output = file.path(tempdir(), "ClimaRep")
  output_files <- list.files(dir_output, recursive = TRUE)
  expect_true(length(output_files) > 0, info = "No files were written to the output directory.")
})
