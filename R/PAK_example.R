set.seed(2458)
n_cells <- 100 * 100
r_clim_present <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
values(r_clim_present) <- c((rowFromCell(r_clim_present, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
                            (rowFromCell(r_clim_present, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
                            (colFromCell(r_clim_present, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
                            (colFromCell(r_clim_present, 1:n_cells) + (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                            (colFromCell(r_clim_present, 1:n_cells) / (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                            (colFromCell(r_clim_present, 1:n_cells) * (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
                            (colFromCell(r_clim_present, 1:n_cells) * (colFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4))

)
names(r_clim_present) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
terra::crs(r_clim_present) <- "EPSG:4326"
terra::plot(r_clim_present)
r_clim_present_filtered <- ClimaRep::vif_filter(r_clim_present, th = 5)
hex_grid <- sf::st_sf(st_make_grid(st_as_sf(as.polygons(terra::ext(r_clim_present))), square = FALSE))
sf::st_crs(hex_grid) <- "EPSG:4326"
polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
polygons$name <- c("Pol_1", "Pol_2")
sf::st_crs(polygons) <- sf::st_crs(hex_grid)
study_area_polygon <- sf::st_as_sf(terra::as.polygons(terra::ext(r_clim_present)))
terra::st_crs(study_area_polygon) <- "EPSG:4326"


terra::plot(r_clim_present[[1]])
terra::plot(polygons, add = TRUE, color= "transparent", lwd = 3)
terra::plot(study_area_polygon, add = TRUE, col = "transparent", lwd = 3, border = "red")

ClimaRep::mh_rep(
polygon = polygons,
col_name = "name",
climatic_variables = r_clim_present_filtered,
th = 0.9,
dir_output = tempdir(),
save_intermediate_raster = TRUE
)

list.files(tempdir())
list.files(file.path(tempdir(), "Charts"))

mh_rep_raw <- rast(list.files(file.path(tempdir(), "MahalanobisRaw"),  pattern = "\\.tif$", full.names = TRUE))
terra::plot(mh_rep_raw[[1]])
terra::plot(polygons[1,], add = TRUE, color= "transparent", lwd = 3)

mh_rep_result <- rast(list.files(file.path(tempdir(), "Representativeness"),  pattern = "\\.tif$", full.names = TRUE))
terra::plot(mh_rep_result[[1]])
terra::plot(polygons[1,], add = TRUE, color= "transparent", lwd = 3)

list.files(file.path(tempdir(), "Representativeness"))

r_clim_future <- r_clim_present_filtered + 2
names(r_clim_future) <- names(r_clim_present_filtered)
terra::crs(r_clim_future) <- terra::crs(r_clim_present_filtered)
plot(r_clim_future)

mh_rep_ch(
polygon = polygons,
col_name = "name",
present_climatic_variables = r_clim_present_filtered,
future_climatic_variables = r_clim_future,
study_area = study_area_polygon,
th = 0.95,
model = "MIROC6",
year = "2070",
dir_output = tempdir(),
save_intermediate_raster = TRUE)
