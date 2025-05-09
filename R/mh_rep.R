#' @title Multivariate Climatic Representativeness Analysis
#'
#' @description Calculates Mahalanobis-based climatic representativeness for input polygons within a defined area, using climate data from a single time period (e.g., present).
#' Representativeness is assessed by comparing the multivariate climate conditions of each cell to the reference climate space defined by the climate conditions *within* that specific polygon.
#'
#' @param polygon An sf object containing the analysis regions (polygons).
#' @param col_name Character. Name of the column in the `polygon` object that contains unique identifiers for each polygon.
#' @param climatic_variables SpatRaster. A raster stack of climate variables representing the conditions of the analysis period (preferably not strongly correlated).
#' @param th Numeric (0-1). Quantile threshold used to define representativeness. Cells with a Mahalanobis distance below or equal to the distance corresponding to this quantile are classified as representative (default: 0.95).
#' @param dir_output Character. Path to the directory where output files will be saved. The function will create subdirectories within this path (default: "output_representativeness/").
#' @param save_intermediate_raster Logical. If TRUE, saves the intermediate continuous Mahalanobis distance rasters calculated for each polygon before binary classification. The final binary classification rasters are always saved (default: FALSE).
#'
#' @return Invisibly returns NULL. Writes the following outputs to disk within subdirectories of `dir_output`:
#' \itemize{
#'   \item Classification GeoTIFF rasters: Binary rasters (typically coded 1 for Representative, 0 for Non-representative) for each input polygon are saved in the `Representativeness/` subdirectory.
#'   \item Visualization Maps: JPEG (or PNG) image files visualizing the classification results for each polygon are saved in the `Charts/` subdirectory.
#'   \item Intermediate Mahalanobis Distance Rasters: *Optionally* saved as GeoTIFF files in the `Mh_Raw/` subdirectory if `save_intermediate_raster = TRUE`.
#' }
#'
#' @details
#' This function performs a multivariate analysis using Mahalanobis distance to assess
#' the climatic representativeness of input polygons based on climate data from a single time period.
#' It evaluates how well the climate conditions of locations within and outside each polygon match the range of multivariate conditions found *within* that particular polygon.
#'
#' Key workflow steps include:
#' \enumerate{
#'   \item **CRS Harmonization:** Ensures all spatial inputs (`polygon`, `climatic_variables`) share the same Coordinate Reference System (CRS), using the CRS of `climatic_variables` as the reference.
#'   \item **Per-Polygon Processing:** For each polygon in the `polygon` object:
#'   \itemize{
#'     \item Clips and masks the climate variables raster (`climatic_variables`) to the polygon's extent and boundary.
#'     \item Calculates the multivariate mean and covariance matrix using the climate data from the clipped and masked raster (handling NA values). This defines the reference climate space specific to that polygon.
#'     \item Calculates the Mahalanobis distance for each cell with valid data within the polygon's extent, relative to the polygon's calculated mean and covariance matrix.
#'     \item Determines a threshold (`th_value`) based on the `th` quantile of the Mahalanobis distances calculated *within the polygon itself*.
#'     \item Classifies each cell within the polygon's extent as "Representative" (Mahalanobis distance $\le$ threshold, typically assigned value 1) or "Non-Representative" (distance $>$ threshold, typically assigned value 0).
#'   }
#'   \item **Output Generation:** Saves the resulting binary classification raster for each polygon as a GeoTIFF file and generates a corresponding visualization map (JPEG/PNG). These are saved within the specified output directory structure.
#' }
#'
#' It is important to note that Mahalanobis distance assumes multivariate normality and is sensitive to collinearity among variables. While the covariance matrix accounts for correlations, it is strongly recommended that the `climatic_variables` are not strongly correlated. Consider performing a collinearity analysis (e.g., using Variance Inflation Factor - VIF) beforehand, perhaps using the `vif_filter` function from this package.
#'
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor values ncell
#' @importFrom sf st_crs st_transform st_geometry st_as_sf st_make_grid
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave element_text
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile na.omit
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' set.seed(2458)
#' n_cells <- 100 * 100
#' r_clim_present <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
#' values(r_clim_present) <- c((rowFromCell(r_clim_present, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
#'                             (rowFromCell(r_clim_present, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
#'                             (colFromCell(r_clim_present, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
#'                             (colFromCell(r_clim_present, 1:n_cells) + (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
#'                             (colFromCell(r_clim_present, 1:n_cells) / (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
#'                             (colFromCell(r_clim_present, 1:n_cells) * (rowFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)),
#'                             (colFromCell(r_clim_present, 1:n_cells) * (colFromCell(r_clim_present, 1:n_cells))* 0.1 + rnorm(n_cells, 0, 4)))
#' names(r_clim_present) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
#' terra::crs(r_clim_present) <- "EPSG:4326"
#' terra::plot(r_clim_present)
#' r_clim_present_filtered <- ClimaRep::vif_filter(r_clim_present, th = 5)
#' hex_grid <- sf::st_sf(
#' sf::st_make_grid(
#'   sf::st_as_sf(
#'     terra::as.polygons(
#'       terra::ext(r_clim_present))),
#'   square = FALSE))
#' sf::st_crs(hex_grid) <- "EPSG:4326"
#' polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
#' polygons$name <- c("Pol_1", "Pol_2")
#' sf::st_crs(polygons) <- sf::st_crs(hex_grid)
#' study_area_polygon <- sf::st_as_sf(as.polygons(terra::ext(r_clim_present)))
#' sf::st_crs(study_area_polygon) <- "EPSG:4326"
#' terra::plot(r_clim_present[[1]])
#' terra::plot(polygons, add = TRUE, color= "transparent", lwd = 3)
#' terra::plot(study_area_polygon, add = TRUE, col = "transparent", lwd = 3, border = "red")
#' ClimaRep::mh_rep(
#' polygon = polygons,
#' col_name = "name",
#' climatic_variables = r_clim_present_filtered,
#' th = 0.9, # Use a threshold, e.g., 90th percentile
#' dir_output = tempdir(),
#' save_intermediate_raster = TRUE)
#'
#' @export
#'
mh_rep <- function(polygon,
                   col_name,
                   climatic_variables,
                   th = 0.95,
                   dir_output = tempdir(),
                   save_intermediate_raster = FALSE) {
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))
  if (!inherits(polygon, "sf"))
    stop("Parameter 'polygon' must be an sf object.")
  if (!inherits(climatic_variables, "SpatRaster"))
    stop("Parameter 'climatic_variables' must be a SpatRaster object.")
  if (!is.character(col_name) ||
      length(col_name) != 1 || !(col_name %in% names(polygon))) {
    stop("Parameter 'col_name' must be a single character string naming a column in 'polygon'.")
  }
  if (!is.numeric(th) || length(th) != 1 || th < 0 || th > 1) {
    stop("Parameter 'th' must be a single numeric value between 0 and 1.")
  }
  if (!is.character(dir_output) || length(dir_output) != 1) {
    stop("Parameter 'dir_output' must be a single character string.")
  }
  if (terra::nlyr(climatic_variables) < 2) {
    warning(
      "climatic_variables has fewer than 2 layers. Mahalanobis distance is typically for multiple variables. Proceeding with single variable analysis if applicable."
    )
  }
  dir_rep <- file.path(dir_output, "Representativeness")
  dir_charts <- file.path(dir_output, "Charts")
  dirs_to_create <- c(dir_rep, dir_charts)
  if (save_intermediate_raster) {
    dir_mh_raw <- file.path(dir_output, "Mh_Raw")
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
  cov_matrix <- cov(data_p[, 3:(ncol(data_p) - 1)], use = "complete.obs")
  for (j in 1:nrow(polygon)) {
    pol <- polygon[j, ]
    pol_name <- as.character(pol[[col_name]])
    pol_name <- as.character(tolower(pol_name))
    pol_name <- iconv(pol_name, to = 'ASCII//TRANSLIT')
    pol_name <- gsub("[^a-z0-9_]+", "_", pol_name)
    pol_name <- gsub("__+", "_", pol_name)
    pol_name <- gsub("^_|_$", "", pol_name)
    message("\nProcessing polygon: ",
            pol_name,
            " (",
            j,
            " of ",
            nrow(polygon),
            ")")
    raster_polygon <- terra::mask(terra::crop(climatic_variables, pol), pol)
    if (all(is.na(terra::values(raster_polygon)))) {
      warning("No available data for: ", pol_name, ". Skipping...")
      next
    }
    mu <- terra::global(raster_polygon, "mean", na.rm = TRUE)$mean
    calculate_mh <- function(data) {
      coords <- data[, 1:2]
      climatic_data <- as.matrix(data[, 3:(ncol(data) - 1)])
      mh_values <- mahalanobis(climatic_data, mu, cov_matrix)
      terra::rast(cbind(coords, mh_values),
                  type = "xyz",
                  crs = reference_system)
    }
    mh_present <- calculate_mh(data_p)
    if (save_intermediate_raster) {
      terra::writeRaster(mh_present, file.path(dir_mh_raw, paste0("MH_PRESENT_", pol_name, ".tif")), overwrite = TRUE)
    }
    mh_poly <- terra::mask(mh_present, pol)
    th_value <- quantile(terra::values(mh_poly),
                         probs = th,
                         na.rm = TRUE)
    if (anyNA(th_value)) {
      warning("No threshold was obtained for: ", pol_name, ". Skipping...")
      next
    }
    classify_mh <- function(mh_raster, threshold) {
      terra::ifel(mh_raster <= threshold, 1, 0)
    }
    th_present <- classify_mh(mh_present, th_value)
    raster_final <- th_present
    terra::writeRaster(raster_final,
                       file.path(
                         dir_output,
                         "Representativeness",
                         paste0("TH_REP_", pol_name, ".tif")
                       ),
                       overwrite = TRUE)
    raster_final <- terra::as.factor(raster_final)
    p <- suppressMessages(
      ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = raster_final) +
        ggplot2::geom_sf(
          data = pol,
          color = "black",
          fill = NA
        ) +
        ggplot2::scale_fill_manual(
          name = " ",
          values = c("0" = "grey90", "1" = "aquamarine4"),
          labels = c("0" = "Non-representative", "1" = "Representative"),
          na.value = "transparent",
          na.translate = FALSE,
          drop = FALSE
        ) +
        ggplot2::ggtitle(pol_name) +
        ggplot2::theme_minimal()
    )
    ggplot2::ggsave(
      filename = file.path(
        dir_output,
        "Charts",
        paste0(pol_name, "_representativeness.jpeg")
      ),
      plot = p,
      width = 10,
      height = 10,
      dpi = 300
    )
  }
  message("\nAll processes were completed")
  cat(paste("\nOutput files in:", dir_output))
  return(invisible(NULL))
}
