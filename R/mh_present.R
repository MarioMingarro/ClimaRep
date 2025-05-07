#' @title Multivariate Climatic Representativeness Analysis for a Single Period
#'
#' @description Calculates Mahalanobis-based climatic representativeness for polygons
#' within a defined area, using data from a single time period. Representativeness
#' is assessed by comparing the climate conditions of each pixel within a polygon
#' against the reference climate space defined by that same polygon.
#'
#' @param polygon An sf object containing the analysis regions (polygons).
#' @param col_name Character. Name of the column in the `polygon` object that contains unique identifiers for each polygon.
#' @param climatic_variables SpatRaster. A raster stack of climate variables (preferably not strongly correlated) representing the conditions of the analysis period (e.g., the present).
#' @param th Numeric (0-1). Quantile threshold used to define representativeness. Pixels with a Mahalanobis distance below the distance corresponding to this quantile are considered representative (default: 0.95).
#' @param dir_output Character. Path to the directory where output files will be saved. The function will create subdirectories within this path (default: "output_representativeness/").
#' @param save_intermediate_raster Logical. If TRUE, saves intermediate Mahalanobis distance rasters calculated per polygon before classification. The final classification rasters are always saved (default: FALSE).
#'
#' @return Invisibly returns NULL. Writes to disk:
#' - Classification GeoTIFF rasters (representing Representative and Non-representative areas) for each polygon in the `Representativeness/` subdirectory within `dir_output`.
#' - JPEG maps visualizing the classification results for each polygon in the `Charts/` subdirectory within `dir_output`.
#' - *Optionally*, intermediate Mahalanobis distance rasters (if `save_intermediate_raster = TRUE`).
#'
#' @details
#' This function performs a multivariate analysis using Mahalanobis distance to assess
#' the climatic representativeness of areas within specified polygons for a single time period.
#'
#' Key workflow steps:
#' 1.  **CRS Harmonization:** Ensures all spatial inputs (`polygon`, `climatic_variables`) share the same Coordinate Reference System (CRS), using the CRS of `climatic_variables` as the reference.
#' 2.  **Per-Polygon Processing:** For each polygon in the `polygon` object:
#'     * Clips and masks the climate variables raster (`climatic_variables`) to the polygon's area.
#'     * Calculates the multivariate mean and covariance matrix using the data from the clipped and masked raster (ignoring NAs). This defines the reference climate space *for that specific polygon*.
#'     * Calculates the Mahalanobis distance for *each pixel with valid data* within the polygon, using the mean and covariance matrix calculated for that polygon.
#'     * Calculates the threshold (`th_value`) as the `th` quantile of the Mahalanobis distances calculated within the polygon.
#'     * Classifies each pixel within the polygon as "Representative" (distance <= threshold, value 1) or "Non-Representative" (distance > threshold, value 0).
#' 3.  **Output Generation:** Generates a classification raster (GeoTIFF) for each polygon and a corresponding visualization map (JPEG). Files are saved within the directory structure specified in `dir_output`.
#'
#' It is important that the `climatic_variables` are not strongly correlated, as Mahalanobis distance assumes independence or uses the covariance matrix to account for correlation. Consider performing a collinearity analysis (e.g., using VIF) beforehand.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' # Create raster data with 3 layers
#' r_clim <- rast(ncols=20, nrows=20)
#' values(r_clim) <- 1:ncell(r_clim)
#' r_clim$var1 <- values(r_clim) * 0.1 + rnorm(ncell(r_clim), 0, 2)
#' r_clim$var2 <- values(r_clim) * 0.05 + rnorm(ncell(r_clim), 0, 1)
#' r_clim$var3 <- rnorm(ncell(r_clim), 10, 3)
#' names(r_clim) <- c("id", "temp", "prec", "elev") # Example variable names
#'
#' # Create two analysis polygons
#' polygons_sf <- st_make_grid(st_as_sf(as.polygons(ext(r_clim))), n = 2) %>%
#'   st_sf(name = c("AreaA", "AreaB"))
#'
#' # Define an output directory
#' dir_out <- "results_present_analysis"
#' if (!dir.exists(dir_out)) dir.create(dir_out)
#'
#' # Run the representativeness analysis for the 'present' climate data
#' pa_mh_representativeness(
#'   polygon = polygons_sf,
#'   col_name = "name",
#'   climatic_variables = r_clim,
#'   th = 0.9, # Use a threshold, e.g., 90th percentile
#'   dir_output = dir_out,
#'   save_intermediate_raster = FALSE
#' )
#'
#' # Check the output files
#' cat("Output files in Representativeness/:\n")
#' print(list.files(file.path(dir_out, "Representativeness")))
#' cat("\nOutput files in Charts/:\n")
#' print(list.files(file.path(dir_out, "Charts")))
#' }
#'
#' @export
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor values ncell
#' @importFrom sf st_crs st_transform st_geometry st_as_sf st_make_grid
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave element_text
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile na.omit
#'
mh_representativeness <- function(polygon,
                                     col_name,
                                     climatic_variables,
                                     th = 0.95,
                                     dir_output = "output_representativeness/",
                                     save_intermediate_raster = FALSE) {
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))

  if (!inherits(polygon, "sf"))
    stop("Parameter 'polygon' must be an sf object.")
  if (!inherits(climatic_variables, "SpatRaster"))
    stop("Parameter 'climatic_variables' must be a SpatRaster object.")
  if (!is.character(col_name) || length(col_name) != 1 || !(col_name %in% names(polygon))) {
    stop("Parameter 'col_name' must be a single character string naming a column in 'polygon'.")
  }
  if (!is.numeric(th) || length(th) != 1 || th < 0 || th > 1) {
    stop("Parameter 'th' must be a single numeric value between 0 and 1.")
  }
  if (!is.character(dir_output) || length(dir_output) != 1) {
    stop("Parameter 'dir_output' must be a single character string.")
  }
  if (terra::nlyr(climatic_variables) < 2) {
    warning("climatic_variables has fewer than 2 layers. Mahalanobis distance is typically for multiple variables. Proceeding with single variable analysis if applicable.")
  }

  dir_rep <- file.path(dir_output, "Representativeness")
  dir_charts <- file.path(dir_output, "Charts")
  dirs_to_create <- c(dir_rep, dir_charts)

  if (save_intermediate_raster) {
    dir_mh_raw <- file.path(dir_output, "MahalanobisRaw")
    dirs_to_create <- c(dirs_to_create, dir_mh_raw)
  }

  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  })


  message("Validating and adjusting Coordinate Reference Systems (CRS)...")
  reference_system_check <- terra::crs(climatic_variables, describe = TRUE)$code
  reference_system <- terra::crs(climatic_variables)

  if (sf::st_crs(polygon)$epsg != reference_system_check) {
    message("Adjusting CRS of polygon to match reference system.")
    polygon <- sf::st_transform(polygon, reference_system)
  }


  message("Starting per-polygon processing...")

  data_p <- na.omit(terra::as.data.frame(climatic_variables, xy = TRUE))
  data_p$Period <- "Present"
  cov_matrix <- cov(data_p[, 3:(ncol(data_p)-1)], use = "complete.obs")

  for(j in 1:nrow(polygon)) {
    pol <- polygon[j, ]
    pol_name <- as.character(pol[[col_name]])
    pol_name <- as.character(tolower(pol_name))
    pol_name <- iconv(pol_name, to = 'ASCII//TRANSLIT')
    pol_name <- gsub("[^a-z0-9_]+", "_", pol_name)
    pol_name <- gsub("__+", "_", pol_name)
    pol_name <- gsub("^_|_$", "", pol_name)

    message("\nProcessing polygon: ", pol_name, " (", j, " of ", nrow(polygon), ")")

    raster_polygon <- terra::mask(terra::crop(climatic_variables, pol), pol)
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

    if(save_intermediate_raster) {
      terra::writeRaster(mh_present,
                         file.path(dir_mh_raw, paste0("MH_PRESENT_", pol_name, ".tif")),
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

    raster_final <- th_present

    terra::writeRaster(
      raster_final,
      file.path(dir_output, "Representativeness", paste0("TH_REP_", pol_name, ".tif")),
      overwrite = TRUE
    )

    raster_final <- terra::as.factor(raster_final)

    p <- suppressMessages(ggplot2::ggplot() +
                            tidyterra::geom_spatraster(data = raster_final) +
                            ggplot2::geom_sf(data = pol, color = "black", fill = NA) +
                            ggplot2::scale_fill_manual(
                              name = " ",
                              values = c(
                                "0" = "grey90",
                                "1" = "aquamarine4"),
                              labels = c(
                                "0" = "Non-representative",
                                "1" = "Representative"
                              ),
                              na.value = "transparent",
                              na.translate = FALSE,
                              drop = FALSE
                            ) +
                            ggplot2::ggtitle(pol_name) +
                            ggplot2::theme_minimal())

    ggplot2::ggsave(
      filename = file.path(dir_output, "Charts", paste0(pol_name, "_rep.jpeg")),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }

  message("\nAll processes were completed")
  cat(paste("\nOutput files in: ",getwd(),"/", dir_output))
  return(invisible(NULL))
}
