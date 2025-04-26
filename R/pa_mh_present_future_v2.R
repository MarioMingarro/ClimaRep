pa_mh_present_future <- function(polygon,
                                 col_name,
                                 present_climatic_variables,
                                 future_climatic_variables,
                                 study_area,
                                 th = .9,
                                 model,
                                 year,
                                 dir_output = "output/",
                                 save_raster = TRUE
) {

  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))

  # Verificar y crear directorios si no existen
  dir_present <- file.path(dir_output, "Present")
  dir_future <- file.path(dir_output, "Future")
  dir_shared <- file.path(dir_output, "Shared")
  dir_charts <- file.path(dir_output, "Charts")

  dirs_to_create <- if(save_raster) {
    c(dir_present, dir_future, dir_charts, dir_shared)
  } else {
    c(dir_charts, dir_shared)
  }

  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  })


  # Verificar y ajustar sistemas de referencia coordenados (CRS)
  reference_system <- terra::crs(present_climatic_variables)

  if(terra::crs(future_climatic_variables) != reference_system){
    message("Ajustando CRS de future_climatic_variables para que coincida")
    future_climatic_variables <- terra::project(future_climatic_variables, reference_system)
  }

  if(sf::st_crs(polygon) != reference_system){
    message("Ajustando CRS de polygon para que coincida")
    polygon <- sf::st_transform(polygon, reference_system)
  }

  if(sf::st_crs(study_area) != reference_system){
    message("Ajustando CRS de study_area para que coincida")
    study_area <- sf::st_transform(study_area, reference_system)
  }

  # Verificar que polygon tenga la columna especificada
  if(!col_name %in% names(polygon)){
    stop("El objeto polygon debe contener una columna llamada '", col_name, "' con los nombres de los polígonos")
  }

  # Iterar sobre cada polígono
  for(j in 1:nrow(polygon)) {
    pol <- polygon[j, ]
    pol_name <- as.character(pol[[col_name]])  # Acceso correcto a la columna dinámica

    message("\nProcesando polígono: ", pol_name, " (", j, " de ", nrow(polygon), ")")

    # Función interna para procesar datos climáticos
    process_climate_data <- function(climatic_variables, poly) {
      # Extraer y procesar datos para el polígono
      raster_polygon <- terra::mask(terra::crop(climatic_variables, poly), poly)
      data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
      data_polygon <- na.omit(data_polygon)
      return(data_polygon)
    }

    # Procesar datos presentes para el polígono
    data_polygon <- process_climate_data(present_climatic_variables, pol)

    # Verificar que hay datos en el polígono
    if(nrow(data_polygon) == 0) {
      warning("No hay datos disponibles para el polígono: ", pol_name, ". Saltando...")
      next
    }

    # Convertir datos presentes a dataframe
    data_p <- terra::as.data.frame(present_climatic_variables, xy = TRUE)
    data_p <- na.omit(data_p)
    data_p$Period <- "Present"

    # Convertir datos futuros a dataframe
    data_f <- terra::as.data.frame(future_climatic_variables, xy = TRUE)
    data_f <- na.omit(data_f)
    data_f$Period <- "Future"

    # Combinar datos presentes y futuros
    data_p_f <- rbind(data_p, data_f)

    # Calcular matriz de covarianza una sola vez (mejor eficiencia)
    cov_matrix <- cov(data_p_f[, 3:(ncol(data_p_f)-1)], use = "complete.obs")

    # Función para calcular Mahalanobis y crear rasters
    calculate_mh_raster <- function(data, period_name, output_dir, plot_title) {
      # Calcular distancia Mahalanobis
      mh_values <- mahalanobis(data[, 3:(ncol(data)-1)],
                               colMeans(data_polygon[, 3:ncol(data_polygon)]),
                               cov_matrix)

      # Crear dataframe con resultados
      mh_df <- cbind(data[, c(1:2, ncol(data))], mh = mh_values)

      # Convertir a raster
      mh_raster <- terra::rast(mh_df[, c(1, 2, 4)], crs = reference_system)
      names(mh_raster) <- "mh_distance"

      # Guardar raster
      if(save_raster) {
        if(period_name == "Present") {
          terra::writeRaster(mh_raster,
                             paste0(output_dir, "/MH_PRESENT_", pol_name, ".tif"),
                             overwrite = TRUE)
        } else if(period_name == "Future") {
          terra::writeRaster(mh_raster,
                             paste0(output_dir, "/MH_FUTURE_", model, "_", year, "_", pol_name, ".tif"),
                             overwrite = TRUE)
        }
      }

      return(mh_raster)
    }

    # Procesar periodo presente
    mh_raster_p <- calculate_mh_raster(data_p, "Present", dir_present,
                                       paste0("Present MH distance of ", pol_name))

    # Procesar periodo futuro
    mh_raster_f <- calculate_mh_raster(data_f, "Future", dir_future,
                                       paste0(year, " ", model, " Future MH distance of ", pol_name))

    # Función para calcular umbral y crear raster binario
    calculate_threshold_raster <- function(mh_raster, poly, threshold_value, output_dir, file_prefix) {
      # Convertir a puntos y calcular umbral
      puntos_todos <- terra::as.points(mh_raster)
      puntos_todos <- sf::st_as_sf(puntos_todos)
      colnames(puntos_todos) <- c("mh", "geometry")
      puntos_dentro <- sf::st_intersection(puntos_todos, poly)

      # Aplicar umbral
      puntos_todos$th <- as.numeric(puntos_todos$mh <= threshold_value)

      # Crear raster binario
      res <- terra::res(mh_raster)
      bbox <- terra::ext(mh_raster)
      raster_template <- terra::rast(ext = bbox,
                                     resolution = res,
                                     crs = terra::crs(mh_raster))

      puntos_vect <- terra::vect(puntos_todos)
      raster_th <- terra::rasterize(puntos_vect, raster_template, field = "th")

      # Guardar raster
      if(save_raster) {
        terra::writeRaster(raster_th,
                           paste0(output_dir, "/", file_prefix, "_", pol_name, ".tif"),
                           overwrite = TRUE)
      }

      return(raster_th)
    }

    # Calcular umbral usando solo datos presentes dentro del polígono
    puntos_presente <- terra::as.points(mh_raster_p)
    puntos_presente <- sf::st_as_sf(puntos_presente)
    puntos_dentro_presente <- sf::st_intersection(puntos_presente, pol)

    if(nrow(puntos_dentro_presente) == 0) {
      warning("No hay puntos presentes dentro del polígono para calcular umbral: ", pol_name, ". Saltando...")
      next
    }

    th_value <- quantile(puntos_dentro_presente$mh, probs = th, na.rm = TRUE)

    # Crear rasters binarios para presente y futuro
    raster_th_p <- calculate_threshold_raster(mh_raster_p, pol, th_value, dir_present, "TH_MH_PRESENT")
    raster_th_f <- calculate_threshold_raster(mh_raster_f, pol, th_value, dir_future, paste0("TH_MH_FUTURE_", model, "_", year))

    # Calcular áreas de representatividad compartida
    compartida <- raster_th_p * raster_th_f
    solo_presente <- (raster_th_p == 1) * (raster_th_f == 0)
    solo_futura <- (raster_th_p == 0) * (raster_th_f == 1)

    # Crear raster categorico final
    raster_th_p_f <- compartida + (solo_presente * 2) + (solo_futura * 3)
    resultado_factor <- as.factor(raster_th_p_f)

    # Crear y guardar gráfico
    l <- ggplot2::ggplot() +
      geom_spatraster(data = resultado_factor) +
      geom_sf(data = study_area, color = "gray50", fill = NA, linewidth = 1.5) +
      geom_sf(data = pol, color = "black", fill = NA) +
      scale_fill_manual(
        values = c("0" = "grey90",
                   "1" = "gold",
                   "2" = "aquamarine3",
                   "3" = "coral3"),
        labels = c("0" = "Non-representativeness areas",
                   "1" = "Stable representativeness",
                   "2" = "Present representativeness",
                   "3" = "Future representativeness"),
        na.value = "transparent",
        na.translate = FALSE
      ) +
      coord_sf() +
      theme_minimal() +
      labs(title = pol_name, fill = " ")

    ggsave(filename = paste0(dir_charts, "/", pol_name, "_rep_shared.jpeg"),
           plot = l, width = 10, height = 8, dpi = 300)

    # Guardar raster final
    terra::writeRaster(raster_th_p_f,
                       paste0(dir_shared, "/TH_SHARED_", model, "_", year, "_", pol_name, ".tif"),
                       overwrite = TRUE)

    message("Procesamiento completado para: ", pol_name)
  }

  message("\nProcesamiento completado para todos los polígonos")
}
