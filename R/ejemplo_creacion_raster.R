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
study_area <- st_as_sf(as.polygons(terra::ext(r_clim)))
st_crs(study_area) <- "EPSG:4326"
hex_grid <- st_make_grid(study_area, square = FALSE) %>%
  st_sf()
polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
polygons$name <- c("Polygon 1", "Polygon 2")
st_crs(polygons) <- st_crs(study_area)


plot(r_clim[[1]])
plot(polygons, add = TRUE, col = "transparent")


getwd()
dir_out <- "results_present_analysis"
if (!dir.exists(dir_out)) dir.create(dir_out)


mh_representativeness(
polygon = polygons,
col_name = "name",
climatic_variables = r_clim,
th = 0.95, # Use a threshold, e.g., 90th percentile
dir_output = dir_out,
save_intermediate_raster = FALSE
)
