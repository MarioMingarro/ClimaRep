#' Multivariate temporal climate representativeness analysis
#'
#' @description Calculates Mahalanobis-based climatic representativeness for different periods of time
#' for polygons within a defined study area. Representativeness is assessed by comparing
#' climate conditions in present and future periods against a reference climate space.
#' The function identifies areas that remain stable, lose, or gain representativeness over time.
#'
#' @param polygon An sf object containing the analysis regions (polygons).
#' @param col_name Character. Name of the column in the `polygon` object that contains unique identifiers for each polygon.
#' @param present_climatic_variables SpatRaster. A non correlated raster stack of climate variables representing present conditions.
#' @param future_climatic_variables SpatRaster. A raster stack containing the *same* climate variables as `present_climatic_variables` but representing future projected conditions.
#' @param study_area sf object. A polygon defining the overall boundary of the study area. Used for clipping/masking operations.
#' @param th Numeric (0-1). Quantile threshold used to define representativeness. Pixels with a Mahalanobis distance below the distance corresponding to this quantile are considered representative (default: 0.95).
#' @param model Character. Name of the climate model used for future projections (used in output filenames).
#' @param year Character. Projection year for future climate data (used in output filenames).
#' @param dir_output Character. Path to the directory where output files will be saved. The function will create subdirectories within this path (default: "output/").
#' @param save_raster Logical. If TRUE, saves intermediate rasters generated during the process (e.g., raw Mahalanobis distance). The final classification rasters listed in `@return` are saved regardless of this parameter setting (default: FALSE).
#' Note on parameters `model` and `year`: These parameters are mandatory. They are used in the output filenames and directory structure to facilitate better data management and organization of results.


#' @return Invisibly returns NULL. Writes to disk:
#' - Classification GeoTIFF rasters (representing stable, loss, and new areas) for Present/, Future/, and Shared/ classifications in the `dir_output` subdirectories (saved regardless of `save_raster` setting).
#' - JPEG maps visualizing the classification results for each polygon in the Charts/ subdirectory within `dir_output`.
#' - *Optionally*, intermediate rasters (if `save_raster = TRUE`).
#'
#' @details
#' This function performs a multivariate analysis using Mahalanobis distance to assess and compare
#' the climatic representativeness of areas within specified polygons for present and future periods.
#'
#' Key workflow steps:
#' 1.  **CRS Harmonization:** Ensures all spatial inputs (`polygon`, `present_climatic_variables`, `future_climatic_variables`, `study_area`) share the same Coordinate Reference System (CRS), using the CRS of `present_climatic_variables` as the reference.
#' 2.  **Reference Climate Space:** Calculates the multivariate mean and covariance matrix based on the `present_climatic_variables` and `future_climatic_variables` within study area*. This defines the reference climate space against which representativeness is measured for each polygon.
#' 3.  **Per-Polygon Processing:** For each polygon:
#'     * Clips and masks the climate rasters (`present_climatic_variables` and `future_climatic_variables`) to the polygon's extent and shape.
#'     * Calculates the Mahalanobis distance for *each pixel* within the polygon for *both* the present and future periods, using the mean and covariance matrix calculated in step 2.
#'     * Applies the quantile threshold (`th`) to the Mahalanobis distances to classify each pixel as "representative" or "not representative" for both the present and future periods within that polygon.
#' 4.  **Classification of Change:** Based on the representativeness status in the present and future, each pixel within the polygon is classified into one of the following categories:
#'     * **Stable:** Representative in both present and future periods.
#'     * **Loss:** Representative in the present period, but not in the future.
#'     * **New:** Not representative in the present period, but representative in the future.
#'     * (Optionally mention others like 'Not representative in either period' if they are explicitly handled/mapped).
#' 5.  **Output Generation:** Generates classification rasters (merged or per-polygon) and visualization maps (JPEG) of the results for each polygon and the overall study area. Rasters are saved in GeoTIFF format and maps as JPEG files within the specified output directory structure.
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
#'pa_mh_present_future(
#'polygon = st_read(poly_path),
#'col_name = "region_id",
#'present_climatic_variables = rast(present_path),
#'future_climatic_variables = rast(future_path),
#'study_area = st_read(study_path),
#'th = 0.85,
#'model = "MIROC6",
#'year = "2070",
#'dir_output = "results"
#' )
#' }
#'
#' @export
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor
#' @importFrom sf st_crs st_transform st_geometry
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile
#'
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

  if (save_raster) {
    dirs_to_create <- c(dir_present, dir_future, dir_shared, dir_charts)
  } else {
    dirs_to_create <- c(dir_shared, dir_charts)
  }
  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  })

  reference_system_check <- terra::crs(present_climatic_variables, describe = TRUE)$code
  reference_system <- terra::crs(present_climatic_variables)

  if (terra::crs(future_climatic_variables, describe = TRUE)$code != reference_system_check) {
    message("Adjusting CRS of future_climatic_variables to match reference system.")
    future_climatic_variables <- terra::project(future_climatic_variables, reference_system)
  }
  if (sf::st_crs(polygon)$epsg != reference_system_check) {
    message("Adjusting CRS of polygon to match reference system.")
    polygon <- sf::st_transform(polygon, reference_system)
  }
  if (sf::st_crs(study_area)$epsg != reference_system_check) {
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
    pol_name <- as.character(tolower(pol_name))
    pol_name <- iconv(pol_name, to = 'ASCII//TRANSLIT')
    pol_name <- gsub("[^a-z0-9_]+", "_", pol_name)
    pol_name <- gsub("__+", "_", pol_name)
    pol_name <- gsub("^_|_$", "", pol_name)

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

    raster_final <- terra::as.factor(raster_final)

    p <- suppressMessages(ggplot2::ggplot() +
                            tidyterra::geom_spatraster(data = raster_final) +
                            ggplot2::geom_sf(data = study_area, color = "gray50", fill = NA, linewidth = 1) +
                            ggplot2::geom_sf(data = pol, color = "black", fill = NA) +
                            ggplot2::scale_fill_manual(
                              name = " ",
                              values = c(
                                "0" = "grey90",
                                "1" = "aquamarine4",
                                "2" = "coral1",
                                "3" = "aquamarine2"
                              ),
                              labels = c(
                                "0" = "Non-representative",
                                "1" = "Stable representativeness",
                                "2" = "Loss representativeness",
                                "3" = "New representativeness"
                              ),
                              na.value = "transparent",
                              na.translate = FALSE,
                              drop = FALSE
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
