pa_mh_present_future <- function(polygon,
                                 col_name,
                                 present_climatic_variables,
                                 future_climatic_variables,
                                 study_area,
                                 th = 0.9,
                                 model,
                                 year,
                                 dir_output = "output/",
                                 save_raster = TRUE) {

  # Configuración inicial de warnings y directorios
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))

  # Verificar y crear directorios
  dir_present <- file.path(dir_output, "Present")
  dir_future <- file.path(dir_output, "Future")
  dir_shared <- file.path(dir_output, "Shared")
  dir_charts <- file.path(dir_output, "Charts")
  dirs_to_create <- c(dir_present, dir_future, dir_charts, dir_shared)

  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  })

  # --- Transformaciones CRS Únicas ---
  reference_system <- terra::crs(present_climatic_variables)

  if (terra::crs(future_climatic_variables) != reference_system) {
    message("Ajustando CRS de future_climatic_variables")
    future_climatic_variables <- terra::project(future_climatic_variables, reference_system)
  }
  if (sf::st_crs(polygon) != reference_system) {
    message("Ajustando CRS de polygon")
    polygon <- sf::st_transform(polygon, reference_system)
  }
  if (sf::st_crs(study_area) != reference_system) {
    message("Ajustando CRS de study_area")
    study_area <- sf::st_transform(study_area, reference_system)
  }



  # Preprocesamiento global
  message("Preprocesando datos globales...")

  # Convertir a dataframes una sola vez
  data_p <- na.omit(terra::as.data.frame(present_climatic_variables, xy = TRUE))
  data_p$Period <- "Present"
  data_f <- na.omit(terra::as.data.frame(future_climatic_variables, xy = TRUE))
  data_f$Period <- "Future"
  data_p_f <- rbind(data_p, data_f)
  cov_matrix <- cov(data_p_f[, 3:(ncol(data_p_f)-1)], use = "complete.obs")

  # Procesamiento por polígono
  for(j in 1:nrow(polygon)) {
    pol <- polygon[j, ]
    pol_name <- as.character(pol[[col_name]])

    message("\nProcesando polígono: ", pol_name, " (", j, " de ", nrow(polygon), ")")

    # Extraer datos del polígono
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    if(all(is.na(terra::values(raster_polygon)))) {
      warning("No hay datos disponibles para el polígono: ", pol_name, ". Saltando...")
      next
    }

    # Calcular vector de medias
    mu <- terra::global(raster_polygon, "mean", na.rm = TRUE)$mean

    # Función optimizada para cálculo de Mahalanobis
    calculate_mh <- function(data) {
      coords <- data[, 1:2]
      climatic_data <- as.matrix(data[, 3:(ncol(data)-1)])
      mh_values <- mahalanobis(climatic_data, mu, cov_matrix)
      terra::rast(cbind(coords, mh_values), type = "xyz", crs = reference_system)
    }

    # Calcular distancias Mahalanobis
    mh_present <- calculate_mh(data_p)
    mh_future <- calculate_mh(data_f)

    if(save_raster) {
        terra::writeRaster(mh_present,
                           paste0(dir_present, "/MH_PRESENT_", pol_name, ".tif"),
                           overwrite = TRUE)
        terra::writeRaster(mh_future,
                           paste0(dir_future, "/MH_FUTURE_", model, "_", year, "_", pol_name, ".tif"),
                           overwrite = TRUE)
      }

    # Calcular umbral
    mh_poly <- terra::mask(mh_present, pol)
    th_value <- quantile(terra::values(mh_poly), probs = th, na.rm = TRUE)

    if(anyNA(th_value)) {
      warning("No se pudo calcular umbral para ", pol_name, ". Saltando...")
      next
    }

    # Función de clasificación optimizada
    classify_mh <- function(mh_raster, threshold) {
      terra::ifel(mh_raster <= threshold, 1, 0)
    }

    # Crear rasters binarios
    th_present <- classify_mh(mh_present, th_value)
    th_future <- classify_mh(mh_future, th_value)

    # Calcular áreas de representatividad
    shared <- th_present * th_future
    solo_presente <- th_present - shared
    solo_futura <- th_future - shared
    raster_final <- shared + (solo_presente * 2) + (solo_futura * 3)

    terra::writeRaster(
      raster_final,
      file.path(dir_output, "Shared", paste0("TH_SHARED_", model, "_", year, "_", pol_name, ".tif")),
      overwrite = TRUE
    )


    # Convertir a factor y manejar NAs explícitamente
    raster_final <- terra::as.factor(raster_final)

    # Crear gráfico con manejo explícito de niveles
    p <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = raster_final) +
      ggplot2::geom_sf(data = study_area, color = "gray50", fill = NA, linewidth = 1.5) +
      ggplot2::geom_sf(data = pol, color = "black", fill = NA) +
      ggplot2::scale_fill_manual(
        name = "Categorías",
        values = c(
          "0" = "grey90",
          "1" = "gold",
          "2" = "aquamarine3",
          "3" = "coral3"
        ),
        labels = c(
          "0" = "No representativo",
          "1" = "Representatividad estable",
          "2" = "Representatividad presente",
          "3" = "Representatividad futura"
        ),
        na.value = "transparent",
        drop = FALSE  # Fuerza a mostrar todos los niveles aunque no estén presentes
      ) +
      ggplot2::ggtitle(pol_name) +
      ggplot2::theme_minimal()

    ggplot2::ggsave(
      filename = file.path(dir_output, "Charts", paste0(pol_name, "_rep_shared.jpeg")),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }

  message("\nProcesamiento completado para todos los polígonos")
  return(invisible(NULL))
}
