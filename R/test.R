library(terra)
library(sf)
library(tictoc)
dir_present_climate_data <- "C:/A_TRABAJO/A_CLIMAREP_TEST/DATA/CLIMA/PRESENT/"
dir_future_climate_data <- "C:/A_TRABAJO/A_CLIMAREP_TEST/DATA/CLIMA/FUTURE/GFDL/"

## Load data ----
present_climatic_variables <- terra::rast(list.files(dir_present_climate_data, "\\.tif$", full.names = T))
# Crear un vector con los nombres de las variables a excluir
exclude_vars <- c("bio8", "bio9", "bio18", "bio19")
#
# # Crear un patrón de expresión regular para excluir estas variables
exclude_pattern <- paste0("bio(", paste(gsub("bio", "", exclude_vars), collapse = "|"), ")")
#
# # Seleccionar las variables deseadas
present_climatic_variables <- subset(present_climatic_variables, grep(exclude_pattern, names(present_climatic_variables), invert = TRUE, value = TRUE))

future_climatic_variables <- terra::rast(list.files(dir_future_climate_data, "\\.tif$", full.names = T))
future_climatic_variables <- subset(future_climatic_variables, grep(exclude_pattern, names(future_climatic_variables), invert = TRUE, value = TRUE))


names(present_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio2",
                                       "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7")

names(future_climatic_variables) <- names(present_climatic_variables)

# Reference system ----
terra::crs(present_climatic_variables)
reference_system <- terra::crs(present_climatic_variables)
##############################################################################################################################
##############################################################################################################################


dir_result <- "C:/A_TRABAJO/A_CLIMAREP_TEST/N2000_RESULTS/"
study_area <- read_sf("C:/A_TRABAJO/A_CLIMAREP_TEST/DATA/Peninsula_Iberica_89.shp")
#study_area <- read_sf("C:/A_TRABAJO/A_CLIMAREP_TEST/DATA/MURCIA.shp")
polygon <- read_sf("C:/A_TRABAJO/A_CLIMAREP_TEST/DATA/NATURA_2000.shp")
# polygon <- filter(polygon, DESIG_ENG %in% c( "Nacional Park"
# ))
polygon <- dplyr::filter(polygon, GIS_AREA >= 10)
study_area <- st_transform(study_area, crs(reference_system))
study_area <- st_make_valid(study_area)
polygon <- st_transform(polygon, crs(reference_system))
polygon <- st_make_valid(polygon)

polygon<- st_intersection(st_crop(polygon, st_bbox(study_area)), study_area)





# Crop raster to study area
present_climatic_variables <-  terra::mask (crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(crop(future_climatic_variables,  study_area), study_area)

###########################################
tic()
present_climatic_variables <- vif_filter(present_climatic_variables, th = 10)
toc()
###########################################
future_climatic_variables <- terra::subset(future_climatic_variables, names(future_climatic_variables) %in% names(present_climatic_variables))

polygon <- polygon[1:100,]
tic()
resultados <- for

    pa_mh_present_future(
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



process_rasters_by_category <- function(folder_path, output_filename = "combined_category_counts.tif", category_values = c(1, 2, 3)) {

    # 1. Encontrar todos los archivos GeoTIFF en la carpeta
    # Usamos ignore.case = TRUE para .tif y .TIFF
    raster_files <- list.files(folder_path, pattern = "\\.tif$|\\.tiff$", full.names = TRUE, ignore.case = TRUE)

    if (length(raster_files) == 0) {
        message("No se encontraron archivos GeoTIFF en la carpeta especificada: ", folder_path)
        return(NULL)
    }

    message("Encontrados ", length(raster_files), " archivos raster. Procesando...")

    # Inicializar listas para guardar los rásteres binarios para cada categoría
    # Usaremos una lista de listas, donde cada sublista es para una categoría
    binary_layers_by_category <- vector("list", length(category_values))
    names(binary_layers_by_category) <- paste0("Cat_", category_values)

    # Guardar el primer raster para comparar propiedades espaciales
    first_raster_props <- NULL

    # 2. Procesar cada archivo raster
    for (i in seq_along(raster_files)) {
        file <- raster_files[i]
        message("Procesando: ", basename(file), " (", i, "/", length(raster_files), ")")

        # Intentar leer el raster, manejando posibles errores
        r <- try(rast(file))
        if (inherits(r, "try-error")) {
            warning("No se pudo leer el archivo raster: ", basename(file), ". Saltando.")
            next
        }

        # Asegurarse de que es un raster de una sola capa (si tiene varias, tomar la primera)
        if (nlyr(r) > 1) {
            r <- r[[1]]
            warning("El raster '", basename(file), "' tiene múltiples capas, usando solo la primera.")
        }

        # Para el primer raster válido, guardar sus propiedades espaciales
        if (is.null(first_raster_props)) {
            first_raster_props <- r
        } else {
            # Opcional: Añadir verificaciones de alineación. Terra a menudo maneja esto
            # pero es bueno ser consciente si los rasters no están perfectamente alineados
            # if (!terra::compareCRS(r, first_raster_props)) { warning("Desajuste de CRS para ", basename(file), ". Los resultados pueden no ser fiables."); }
            # if (terra::ext(r) != terra::ext(first_raster_props)) { warning("Desajuste de extensión para ", basename(file), ". Los resultados pueden no ser fiables."); }
            # if (any(terra::res(r) != terra::res(first_raster_props))) { warning("Desajuste de resolución para ", basename(file), ". Los resultados pueden no ser fiables."); }
        }

        # Crear capas binarias (0/1) para cada categoría especificada
        for (k in seq_along(category_values)) {
            cat_value <- category_values[k]
            # ifel(condición, valor si TRUE, valor si FALSE)
            # Si el valor de la celda es igual a la categoría_valor, asignar 1, de lo contrario 0
            # Los NA en el raster original se convierten en 0 aquí
            binary_layer <- terra::ifel(r == cat_value, 1, 0)

            # Añadir la capa binaria a la lista correspondiente a la categoría
            binary_layers_by_category[[k]] <- c(binary_layers_by_category[[k]], binary_layer)
        }
    }

    # Verificar si se procesaron rasters válidos
    if (all(sapply(binary_layers_by_category, is.null))) {
        message("No se procesaron datos raster válidos.")
        return(NULL)
    }

    # 3. Apilar las capas binarias para cada categoría y calcular la suma
    message("Apilando y sumando las capas binarias por categoría...")

    output_layers <- list()

    for (k in seq_along(category_values)) {
        cat_value <- category_values[k]
        cat_name <- names(binary_layers_by_category)[k]
        list_of_layers <- binary_layers_by_category[[k]]

        if (length(list_of_layers) == 0 || is.null(list_of_layers[[1]])) {
            warning("No se generaron capas binarias para la categoría ", cat_value, ". Posiblemente la categoría no estaba presente o hubo errores de lectura.")
            # Crear una capa vacía con las mismas propiedades si no hay datos para esta categoría
            if (!is.null(first_raster_props)) {
                count_layer <- rast(first_raster_props) # Crea un raster con mismas ext/res/crs
                values(count_layer) <- 0 # Llenar con ceros
            } else {
                count_layer <- NULL # No hay propiedades de referencia si ningún raster fue válido
            }
        } else {
            # Apilar las capas binarias de esta categoría
            stack_cat <- if (length(list_of_layers) > 1) rast(list_of_layers) else list_of_layers[[1]]

            # Calcular la suma (conteo) para esta categoría en todos los rasters
            # na.rm = TRUE asegura que los NA resultantes de operaciones intermedias no afecten la suma final
            count_layer <- app(stack_cat, fun = sum, na.rm = TRUE)
        }


        # Asignar nombre a la capa de conteo
        if (!is.null(count_layer)) {
            names(count_layer) <- paste0("Count_Category", cat_value)
            output_layers[[k]] <- count_layer
        } else {
            warning("No se pudo crear la capa de conteo para la categoría ", cat_value)
        }
    }

    # Combinar las capas de conteo resultantes en un solo stack
    if (length(output_layers) == 0) {
        message("No se pudo generar ninguna capa de salida.")
        return(NULL)
    }

    output_stack <- rast(output_layers) # Combina la lista de capas en un stack

    # 4. Escribir el raster de salida
    output_path <- file.path(folder_path, output_filename)
    message("Escribiendo stack de rasters de salida a: ", output_path)
    # Usar datatype="INT2U" o "INT4U" si sabes que el conteo no excederá 255 o 65535 respectivamente
    # Si el conteo puede ser mayor, usar "INT4U" o dejar que writeRaster lo determine
    writeRaster(output_stack, output_path, overwrite = TRUE, datatype = "INT2U")

    message("Procesamiento completado.")

    # Devolver el stack de salida (invisiblemente)
    return(invisible(output_stack))
}

# --- Cómo usar la función ---

# 1. Especifica la ruta a la carpeta con tus archivos .tif
# Asegúrate de que esta carpeta exista y contenga los archivos
# Por ejemplo:
mi_carpeta_con_rasters <- "C:/A_TRABAJO/A_CLIMAREP_TEST/N2000_RESULTS/Shared/"

# 2. Llama a la función
# El archivo de salida "combined_category_counts.tif" se creará dentro de esa misma carpeta
combined_counts_raster <- process_rasters_by_category(mi_carpeta_con_rasters)

# Puedes especificar un nombre diferente para el archivo de salida si lo deseas
#combined_counts_raster <- process_rasters_by_category(mi_carpeta_con_rasters, output_filename = "conteo_categorias_final.tif")

# Si tus categorías tienen valores diferentes a 1, 2, 3, especifícalos:
# combined_counts_raster <- process_rasters_by_category(mi_carpeta_con_rasters, category_values = c(10, 20, 30))


# Después de ejecutarla, el archivo de salida estará en tu_carpeta/combined_category_counts.tif
# Puedes cargarlo y verlo así:
# resultado_final <- rast(file.path(mi_carpeta_con_rasters, "combined_category_counts.tif"))
# plot(resultado_final)
