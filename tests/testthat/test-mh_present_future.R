library(testthat)
library(terra)
library(sf)


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

# --- Configuración de Datos (añadir escenario futuro) ---

r_clim_future <- r_clim + 2 # Añadir un incremento constante de 2 a cada variable
names(r_clim_future) <- names(r_clim) # Mantener los mismos nombres de variables
terra::crs(r_clim_future) <- terra::crs(r_clim) # Mantener el mismo CRS

# Crear un área de estudio (ej. la extensión del raster como un polígono)
study_area_poly <- st_as_sf(as.polygons(terra::ext(r_clim)))
st_crs(study_area_poly) <- "EPSG:4326"


# --- Definir Directorio de Salida Temporal para el Test ---
temp_dir_out_future <- file.path(tempdir(), "test_mh_present_future_output")

# Asegurarse de que el directorio no exista o esté vacío antes del test
unlink(temp_dir_out_future, recursive = TRUE)

# --- Escribir el Test ---
testthat::test_that("mh_present_future runs and writes output", {

  # --- Preparación ---
  # (Los datos r_clim, r_clim_future, polygons y study_area_poly ya están definidos arriba)
  # Definir parámetros adicionales requeridos por la función
  test_model <- "TestModel"
  test_year <- "2050"

  # --- Ejecutar la función ---
  mh_present_future(
    polygon = polygons,
    col_name = "name",
    present_climatic_variables = r_clim,
    future_climatic_variables = r_clim_future,
    study_area = study_area_poly, # Usar el polígono del área de estudio
    th = 0.95,
    model = test_model,
    year = test_year,
    dir_output = temp_dir_out_future,
    save_intermediate_raster = TRUE
  )

  # --- Expectativas (Assertions) ---

  # 1. El directorio de salida debe haber sido creado
  expect_true(dir.exists(temp_dir_out_future))

  # 2. Debe haber archivos escritos en el directorio de salida
  # La función probablemente creará subdirectorios o nombres de archivo basados en el modelo/año.
  # Verificamos que el directorio no esté vacío de forma recursiva.
  output_files <- list.files(temp_dir_out_future, recursive = TRUE)
  expect_true(length(output_files) > 0, info = "No files were written to the output directory.")

  # Opcional: Si sabes los nombres exactos de los archivos esperados, puedes añadir
  # expect_true(file.exists(file.path(temp_dir_out_future, "expected_file_name.tif")))

})

# --- Limpieza (Teardown) ---
testthat::teardown({
  if (dir.exists(temp_dir_out_future)) {
    unlink(temp_dir_out_future, recursive = TRUE)
  }
})
