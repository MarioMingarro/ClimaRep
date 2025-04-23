library(terra)
library(sf)
dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENT/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/FUTURE/GFDL/"
dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/PKG/"
#study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/Peninsula_Iberica_89.shp")
study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MURCIA/MURCIA.shp")
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MURCIA/WDPA_MURCIA_89N.shp")

## Load data ----
present_climatic_variables <- terra::rast(list.files(dir_present_climate_data, "\\.tif$", full.names = T))
# Crear un vector con los nombres de las variables a excluir
exclude_vars <- c("bio8", "bio9", "bio18", "bio19")

# Crear un patrón de expresión regular para excluir estas variables
exclude_pattern <- paste0("bio(", paste(gsub("bio", "", exclude_vars), collapse = "|"), ")")

# Seleccionar las variables deseadas
present_climatic_variables <- subset(present_climatic_variables, grep(exclude_pattern, names(present_climatic_variables), invert = TRUE, value = TRUE))

future_climatic_variables <- terra::rast(list.files(dir_future_climate_data, "\\.tif$", full.names = T))
future_climatic_variables <- subset(future_climatic_variables, grep(exclude_pattern, names(future_climatic_variables), invert = TRUE, value = TRUE))


names(present_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio2",
                                       "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7")

names(future_climatic_variables) <- names(present_climatic_variables)

# Reference system ----
terra::crs(present_climatic_variables)
reference_system <-"EPSG:4326"

study_area <- st_transform(study_area, crs(reference_system))
study_area <- st_make_valid(study_area)
polygon <- st_transform(polygon, crs(reference_system))
polygon <- st_make_valid(polygon)

polygon <- st_intersection(st_crop(polygon, st_bbox(study_area)), study_area)


# Crop raster to study area
present_climatic_variables <-  terra::mask (crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(crop(future_climatic_variables,  study_area), study_area)

###########################################
present_climatic_variables <- vif_filter(present_climatic_variables, th = 10)
###########################################
future_climatic_variables <- terra::subset(future_climatic_variables, names(future_climatic_variables) %in% names(present_climatic_variables$filtered_object))



library(ggplot2)
library(terra)
library(sf)
library(tidyterra)
kk <- polygon[5,]
pp <- pa_mh_present_future(
    present_climatic_variables,
    future_climatic_variables,
    kk,
    output_dir = "results",
    model = "kkdvk",
    year = "2070",
    th = 0.9,
    reference_system = "EPSG:4326" ,
    study_area = NULL)


mm <- calculate_mahalanobis_and_threshold(
    present_climatic_variables,
    future_climatic_variables,
    polygon,
    output_dir = "results",
    model = "NULL",
    year = "2070",
    th = 0.9,
    reference_system = "EPSG:4326",
    study_area = study_area
)
