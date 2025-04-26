#' Multivariate Climate Representativeness Analysis (Present-Future)
#'
#' @description Calculates Mahalanobis-based climatic representativeness for current and future
#' climate scenarios across polygon regions. Identifies stable, present-only, future-only,
#' and non-representative areas.
#'
#' @param polygon An sf object containing analysis regions
#' @param col_name Character. Name of column in polygon with region identifiers
#' @param present_climatic_variables SpatRaster. Current climate variables raster stack
#' @param future_climatic_variables SpatRaster. Future climate projections raster stack
#' @param study_area sf object. Study area boundary polygon
#' @param th Numeric (0-1). Quantile threshold for classification (default: 0.9)
#' @param model Character. Climate model name (used in filenames)
#' @param year Character. Projection year (used in filenames)
#' @param dir_output Character. Output directory path (default: "output/")
#' @param save_raster Logical. Save intermediate rasters? (default: TRUE)
#'
#' @return Invisibly returns NULL. Writes to disk:
#' - GeoTIFF rasters in Present/, Future/, Shared/ subdirectories
#' - JPEG maps in Charts/ subdirectory
#'
#' @details
#' Key workflow steps:
#' 1. CRS harmonization using present_climatic_variables as reference
#' 2. Global covariance matrix calculation from combined present-future data
#' 3. Per-polygon processing: Mahalanobis distance, thresholding, classification
#' 4. Output generation of rasters and visualizations
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' # Example data paths
#' poly_path <- system.file("extdata", "regions.gpkg", package = "yourpackage")
#' present_path <- system.file("extdata", "present_climate.tif", package = "yourpackage")
#' future_path <- system.file("extdata", "future_climate.tif", package = "yourpackage")
#' study_path <- system.file("extdata", "study_area.gpkg", package = "yourpackage")
#'
#' pa_mh_present_future(
#'   polygon = st_read(poly_path),
#'   col_name = "region_id",
#'   present_climatic_variables = rast(present_path),
#'   future_climatic_variables = rast(future_path),
#'   study_area = st_read(study_path),
#'   th = 0.85,
#'   model = "MIROC6",
#'   year = "2070",
#'   dir_output = "results"
#' )
#' }
#'
#' @export
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor
#' @importFrom sf st_crs st_transform st_geometry
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile

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


  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))

  dir_present <- file.path(dir_output, "Present")
  dir_future <- file.path(dir_output, "Future")
  dir_shared <- file.path(dir_output, "Shared")
  dir_charts <- file.path(dir_output, "Charts")
  dirs_to_create <- c(dir_present, dir_future, dir_charts, dir_shared)

  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  })

  reference_system <- terra::crs(present_climatic_variables)

  if (terra::crs(future_climatic_variables) != reference_system) {
    message("Adjusting CRS of future_climatic_variables to match reference system.")
    future_climatic_variables <- terra::project(future_climatic_variables, reference_system)
  }
  if (sf::st_crs(polygon) != reference_system) {
    message("Adjusting CRS of polygon to match reference system.")
    polygon <- sf::st_transform(polygon, reference_system)
  }
  if (sf::st_crs(study_area) != reference_system) {
    message("Adjusting CRS of study_area to match reference system.")
    study_area <- sf::st_transform(study_area, reference_system)
  }

  message("Starting process")

  data_p <- na.omit(terra::as.data.frame(present_climatic_variables, xy = TRUE))
  data_p$Period <- "Present"
  data_f <- na.omit(terra::as.data.frame(future_climatic_variables, xy = TRUE))
  data_f$Period <- "Future"
  data_p_f <- rbind(data_p, data_f)
  cov_matrix <- cov(data_p_f[, 3:(ncol(data_p_f)-1)], use = "complete.obs")


  for(j in 1:nrow(polygon)) {
    pol <- polygon[j, ]
    pol_name <- as.character(pol[[col_name]])

    message("\nProcessing polygon: ", pol_name, " (", j, " of ", nrow(polygon), ")")

    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    if(all(is.na(terra::values(raster_polygon)))) {
      warning("No available data for: ", pol_name, ". Skipping...")
      next
    }

    mu <- terra::global(raster_polygon, "mean", na.rm = TRUE)$mean

    calculate_mh <- function(data) {
      coords <- data[, 1:2]
      climatic_data <- as.matrix(data[, 3:(ncol(data)-1)])
      mh_values <- mahalanobis(climatic_data, mu, cov_matrix)
      terra::rast(cbind(coords, mh_values), type = "xyz", crs = reference_system)
    }

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

    mh_poly <- terra::mask(mh_present, pol)
    th_value <- quantile(terra::values(mh_poly), probs = th, na.rm = TRUE)

    if(anyNA(th_value)) {
      warning("No threshold was obtained for: ", pol_name, ". Skipping...")
      next
    }

    classify_mh <- function(mh_raster, threshold) {
      terra::ifel(mh_raster <= threshold, 1, 0)
    }

    th_present <- classify_mh(mh_present, th_value)
    th_future <- classify_mh(mh_future, th_value)

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
    p <- suppressMessages(ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = raster_final) +
      ggplot2::geom_sf(data = study_area, color = "gray50", fill = NA, linewidth = 1) +
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
          "0" = "Non-representative",
          "1" = "Stable representativeness",
          "2" = "Present representativeness",
          "3" = "Future representativeness"
        ),
        na.value = "transparent",
        drop = FALSE  # Fuerza a mostrar todos los niveles aunque no estén presentes
      ) +
      ggplot2::ggtitle(pol_name) +
      ggplot2::theme_minimal())

    ggplot2::ggsave(
      filename = file.path(dir_output, "Charts", paste0(pol_name, "_rep_shared.jpeg")),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }

  message("\nAll processes were completed")
  return(invisible(NULL))
}
