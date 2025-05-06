library(testthat)
library(terra)
library(sf)

# --- Configuración de Datos ---
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

hex_grid <- st_sf(st_make_grid(st_as_sf(as.polygons(terra::ext(r_clim))), square = FALSE))
st_crs(hex_grid) <- "EPSG:4326"
polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
polygons$name <- c("Polygon_1", "Polygon_2")
st_crs(polygons) <- st_crs(hex_grid)

temp_dir_out <- file.path(tempdir(), "test_mh_present_output")

unlink(temp_dir_out, recursive = TRUE)


testthat::test_that("mh_present runs and writes output", {
  mh_present(
    polygon = polygons,
    col_name = "name",
    climatic_variables = r_clim, # Usar r_clim
    th = 0.95,
    dir_output = temp_dir_out,
    save_intermediate_raster = TRUE
  )

  # 1. El directorio de salida debe haber sido creado
  expect_true(dir.exists(temp_dir_out))

  # 2. Debe haber archivos escritos en el directorio de salida
  output_files <- list.files(temp_dir_out, recursive = TRUE)
  expect_true(length(output_files) > 0, info = "No files were written to the output directory.")

})

# --- Limpieza (Teardown) ---
testthat::teardown({
  if (dir.exists(temp_dir_out)) {
    unlink(temp_dir_out, recursive = TRUE)
  }
})
