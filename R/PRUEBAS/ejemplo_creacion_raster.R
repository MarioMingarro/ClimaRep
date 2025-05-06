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
polygons$name <- c("Polygon 1", "Polygon 2")
st_crs(polygons) <- st_crs(hex_grid)


# plot(r_clim)
# plot(polygons, add = TRUE, col = "transparent")



getwd()
dir_out <- "results_present_analysis"
if (!dir.exists(dir_out)) dir.create(dir_out)


mh_present(
polygon = polygons,
col_name = "name",
climatic_variables = kk,
th = 0.95, # Use a threshold, e.g., 90th percentile
dir_output = dir_out,
save_intermediate_raster = T
)
