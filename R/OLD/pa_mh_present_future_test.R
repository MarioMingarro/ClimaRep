pa_mh_present_future <- function(
    present_climatic_variables,
    future_climatic_variables,
    polygon,
    output_dir = "results",
    model = "NULL",
    year = "2070",
    th = 0.9,
    reference_system = "EPSG:4326" ,
    study_area = NULL) {

  # Crear directorios de salida (usando R base)
  dirs_to_create <- file.path(output_dir, c("present", "future", "charts"))
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }

  # Función para preparar datos (sin dplyr)
  prepare_data <- function(raster_data, period) {
    df <- as.data.frame(raster_data, xy = TRUE)
    df <- df[complete.cases(df), ]
    df <- cbind(df[, 1:2], Period = period, df[, 3:ncol(df)])
    return(df)
  }

  data_present <- prepare_data(present_climatic_variables, "Present")
  data_future <- prepare_data(future_climatic_variables, "Future")

  # Función para calcular distancia Mahalanobis (base R)
  calculate_mh <- function(data, polygon_data) {
    polygon_mask <- mask(crop(present_climatic_variables, polygon_data), polygon_data)
    polygon_df <- as.data.frame(polygon_mask, xy = TRUE)
    polygon_df <- polygon_df[complete.cases(polygon_df), ]

    mahalanobis(
      as.matrix(data[, 4:ncol(data)]),
      colMeans(polygon_df[, 3:ncol(polygon_df)]),
      cov(data[, 4:ncol(data)]),
      inverted = FALSE
    )
  }

  # Procesar cada polígono
  process_polygon <- function(j) {
    pol <- polygon[j, ]
    name <- polygon$NAME[j]

    # Presente
    mh_present <- calculate_mh(data_present, pol)
    mh_df_present <- cbind(data_present[, 1:3], mh = mh_present)
    mh_raster_present <- rast(mh_df_present[, c(1, 2, 4)], crs = reference_system)

    # Calcular umbral
    points_present <- as.points(mh_raster_present)
    points_present_sf <- sf::st_as_sf(points_present)
    points_inside <- sf::st_intersection(points_present_sf, pol)
    th_value <- quantile(points_inside$mh, probs = th, na.rm = TRUE)

    # Aplicar umbral
    points_present_sf$th <- as.numeric(points_present_sf$mh <= th_value)
    raster_th_present <- create_threshold_raster(points_present_sf, mh_raster_present)

    # Futuro
    combined_data <- rbind(data_present, data_future)
    mh_future <- calculate_mh(combined_data, pol)
    mh_df_future <- cbind(combined_data[, 1:3], mh = mh_future)

    # Filtrar solo futuro (base R en lugar de dplyr)
    mh_df_future <- mh_df_future[mh_df_future$Period == "Future", ]
    mh_raster_future <- rast(mh_df_future[, c(1, 2, 4)], crs = reference_system)

    # Aplicar umbral a futuro
    points_future <- as.points(mh_raster_future)
    points_future_sf <- sf::st_as_sf(points_future)
    points_future_sf$th <- as.numeric(points_future_sf$mh <= th_value)
    raster_th_future <- create_threshold_raster(points_future_sf, mh_raster_future)

    # Guardar resultados futuros
    writeRaster(
      mh_raster_future,
      file.path(output_dir, "future", paste0("MH_", model, "_", year, "_", name, ".tif")),
      overwrite = TRUE
    )

    # Análisis comparativo
    compare_scenarios(raster_th_present, raster_th_future, pol, name, output_dir)

    # Análisis de distancia
    distance_analysis(points_present_sf, points_future_sf, pol, name, output_dir, model, year)
  }

  # Función auxiliar para crear raster de umbral
  create_threshold_raster <- function(points_sf, template_raster) {
    res <- res(template_raster)
    bbox <- ext(template_raster)

    raster_template <- rast(
      ext = bbox,
      nrows = round((bbox[4] - bbox[3]) / res[2]),
      ncols = round((bbox[2] - bbox[1]) / res[1])
    )

    rasterized <- rasterize(vect(points_sf), raster_template, field = "th")
    crs(rasterized) <- crs(template_raster)
    return(rasterized)
  }

  # Función para comparar escenarios (requiere ggplot2)
  compare_scenarios <- function(raster_present, raster_future, pol, name, output_dir) {
    shared <- raster_present * raster_future
    only_present <- (raster_present == 1) * (raster_future == 0)
    only_future <- (raster_present == 0) * (raster_future == 1)

    combined <- shared + (only_present * 2) + (only_future * 3)
    combined_factor <- as.factor(combined)

    # Gráfico mínimo (se mantiene ggplot2 para calidad de visualización)
    plot <- ggplot() +
      geom_spatraster(data = combined_factor) +
      geom_sf(data = study_area, color = "gray50", fill = NA, linewidth = 1.5) +
      geom_sf(data = pol, color = "black", fill = NA) +
      scale_fill_manual(
        values = c("0" = "grey90", "1" = "gold", "2" = "aquamarine3", "3" = "coral3"),
        labels = c("Non-representativeness", "Stable", "Present only", "Future only"),
        na.value = "transparent"
      ) +
      coord_sf() +
      theme_minimal() +
      labs(title = name, fill = "")

    png(file.path(output_dir, "charts", paste0(name, "_rep_shared.png")),
        width = 1000, height = 800)
    print(plot)
    dev.off()

    writeRaster(
      combined,
      file.path(output_dir, "future", paste0("TH_MH_PRESENT_", model, "_", year, "_", name, ".tif")),
      overwrite = TRUE
    )
  }

  # Función para análisis de distancia (base R)
  distance_analysis <- function(points_present, points_future, pol, name, output_dir, model, year) {
    # Calcular distancias mínimas
    calc_min_distance <- function(points) {
      dist_matrix <- sf::st_distance(points, pol)
      apply(dist_matrix, 1, min) / 1000  # Convertir a km
    }

    points_present$dist <- round(calc_min_distance(points_present), 0)
    points_future$dist <- round(calc_min_distance(points_future), 0)

    # Agregar coordenadas
    coords <- sf::st_coordinates(points_present)
    points_present <- cbind(as.data.frame(points_present), coords)
    coords <- sf::st_coordinates(points_future)
    points_future <- cbind(as.data.frame(points_future), coords)

    # Calcular estadísticas por distancia (base R)
    calculate_stats <- function(points, suffix) {
      points <- points[!is.na(points$dist), ]
      stats <- aggregate(th ~ dist, data = points,
                         FUN = function(x) c(
                           n_total = length(x),
                           n_zeros = sum(x == 0),
                           n_ones = sum(x == 1),
                           pct_zeros = mean(x == 0) * 100,
                           pct_ones = mean(x == 1) * 100
                         ))
      stats <- cbind(stats[1], as.data.frame(stats[[2]]))
      names(stats) <- paste0(names(stats), "_", suffix)
      return(stats)
    }

    stats_present <- calculate_stats(points_present, "p")
    stats_future <- calculate_stats(points_future, "f")

    # Combinar resultados
    combined_stats <- merge(stats_present, stats_future, by.x = "dist_p", by.y = "dist_f")
    names(combined_stats)[1] <- "dist"

    # Calcular acumulados
    combined_stats$ones_p_accum <- cumsum(combined_stats$n_ones_p)
    combined_stats$ones_f_accum <- cumsum(combined_stats$n_ones_f)

    # Gráfico mínimo (base R)
    png(file.path(output_dir, "charts", paste0(name, "_scenarios_difference.png")),
        width = 1000, height = 800)
    plot(combined_stats$dist, combined_stats$ones_p_accum, type = "l", col = "aquamarine3",
         lwd = 2, xlab = "Distance (km)", ylab = "Cumulative representative cells",
         main = paste(name, year, model, sep = " - "))
    lines(combined_stats$dist, combined_stats$ones_f_accum, col = "coral3", lwd = 2)
    abline(h = max(combined_stats$ones_p_accum), lty = 2, lwd = 1.5)
    legend("bottomright", legend = c("Present", "Future"),
           col = c("aquamarine3", "coral3"), lwd = 2)
    dev.off()
  }

  # Procesar todos los polígonos
  for (j in seq_len(nrow(polygon))) {
    process_polygon(j)
  }
}
