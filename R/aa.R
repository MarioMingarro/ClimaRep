library(terra)
library(sf)
r <- rast(matrix(rnorm(100), 10, 10))
s <- rast(list(r, 2*r + 0.1, r + 0.5, -r))
names(s) <- c("layer1", "layer2", "layer3", "layer4")


dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENT/"
present_climatic_variables <- terra::rast(list.files(dir_present_climate_data, "\\.tif$", full.names = T))

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/MURCIA.shp")
reference_system <-"EPSG:4326"
study_area <- st_transform(study_area, "EPSG:4326" )

# Crop raster to study area
present_climatic_variables <-  terra::mask(terra::crop(present_climatic_variables, study_area), study_area)

result_vif_filter <- vif_filter(present_climatic_variables, th = 10)


print(result_raster)
plot(present_climatic_variables[1,])
