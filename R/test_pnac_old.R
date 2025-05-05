library(terra)
library(sf)
library(tictoc)
dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA_OLD/PRESENTE/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA_OLD/RCP85_2070/"

present_climatic_variables <- terra::rast(list.files(dir_present_climate_data,  full.names = T))
future_climatic_variables <- terra::rast(list.files(dir_future_climate_data,  full.names = T))
names(future_climatic_variables) <- names(present_climatic_variables)

# Reference system ----
terra::crs(present_climatic_variables)
reference_system <- terra::crs(present_climatic_variables)
##############################################################################################################################
##############################################################################################################################


dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/OLD/TEST_PNAC/"
study_area <- read_sf("C:/A_TRABAJO/A_CLIMAREP_TEST/DATA/Peninsula_Iberica_89.shp")
#study_area <- read_sf("C:/A_TRABAJO/A_CLIMAREP_TEST/DATA/MURCIA.shp")
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/OLD/TEST_PNAC/national_parks.shp")

study_area <- st_transform(study_area, crs(reference_system))
study_area <- st_make_valid(study_area)
polygon <- st_transform(polygon, crs(reference_system))
polygon <- st_make_valid(polygon)

polygon<- st_intersection(st_crop(polygon, st_bbox(study_area)), study_area)





# Crop raster to study area
present_climatic_variables <-  terra::mask (crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(crop(future_climatic_variables,  study_area), study_area)
present_climatic_variables  <-  terra::mask(terra::crop(present_climatic_variables, future_climatic_variables), future_climatic_variables)
future_climatic_variables  <-  terra::mask(terra::crop(future_climatic_variables, present_climatic_variables), present_climatic_variables)

###########################################
tic()
present_climatic_variables <- vif_filter(present_climatic_variables, th = 10)
toc()
###########################################
future_climatic_variables <- terra::subset(future_climatic_variables, names(future_climatic_variables) %in% names(present_climatic_variables))


tic()
resultados <- pa_mh_present_future(
  polygon = polygon,
  col_name = "ORIG_NAME",
  present_climatic_variables = present_climatic_variables,
  future_climatic_variables = future_climatic_variables,
  study_area = study_area,
  th = 0.9,
  model = "GFDL",
  year = "2070",
  dir_output = dir_result,
  save_raster = F)
toc()

