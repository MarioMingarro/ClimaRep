#' @title Multivariate Climate Representativeness Analysis
#'
#' @description Calculates Mahalanobis-based climate representativeness for input `polygon` within a defined area.
#' Representativeness is assessed by comparing the multivariate climate conditions of each cell to the reference climate space defined by the climate conditions within `study_area`.
#'
#' @param polygon An `sf` object containing the definied areas.
#' @param col_name `Character`. Name of the column in the `polygon` object that contains unique identifiers for each polygon.
#' @param climate_variables `SpatRaster`. A raster stack of climate variables representing the conditions of the analysis period.
#' @param th `Numeric` (0-1). Percentile threshold used to define representativeness. Cells with a Mahalanobis distance below or equal to the `th` are classified as representative (default: 0.95).
#' @param dir_output `Character`. Path to the directory where output files will be saved. The function will create subdirectories within this path.
#' @param save_raw `Logical.` If `TRUE`, saves the intermediate continuous Mahalanobis distance rasters calculated for each polygon before binary classification. The final binary classification rasters are always saved (default: FALSE).
#'
#' @return Invisibly returns NULL. Writes the following outputs to disk within subdirectories of `dir_output`:
#' \itemize{
#'   \item Classification (`.tif` ) rasters: Binary rasters (`1` for Representative and `0` for Non-representative) for each input `polygon` are saved in the `Representativeness/` subdirectory.
#'   \item Visualization (`.jpeg`) maps: Image files visualizing the classification results for each `polygon` are saved in the `Charts/` subdirectory.
#'   \item Raw Mahalanobis distance rasters: *Optionally* saved as `.tif` files in the `Mh_Raw/` subdirectory if `save_raw = TRUE`.
#' }
#'
#' @details
#' This function performs a multivariate analysis using Mahalanobis distance to assess
#' the climate representativeness of input polygons based on climate data from a single time period.
#'
#' Here are the key steps:
#' \enumerate{
#'   \item Ensure all spatial inputs (`polygon`, `climate_variables`) share the same Coordinate Reference System (CRS), using the CRS of `climate_variables` as the reference.
#'   \item For each polygon in the `polygon` object:
#'   \itemize{
#'     \item Crop and mask the climate variables raster (`climate_variables`) to the boundary of the current polygon.
#'     \item Calculate the multivariate mean and covariance matrix using the climate data from the clipped and masked raster (handling NA values). This defines the reference climate conditions for the current polygon.
#'     \item Calculate the Mahalanobis distance for each cell within the (`study_area`)'s extent relative to the multivariate centroid and covariance matrix calculated for the current polygon.
#'     \item Apply the specified threshold (`th`) to the calculated Mahalanobis distances to determine which cells are considered representative. This threshold is typically a percentile of the Mahalanobis distances calculated for the cells originally within the current polygon.
#'     \item Classify each cell within the (`study_area`)'s extent as `Representative = 1` (mh distance is below the threshold ) or `Non-Representative = 0` (mh distance is above the threshold).
#'   }
#'   \item Output Generation: Saves the binary classification raster (`.tif`) for each polygon and generates a corresponding visualization map (`jpeg`). These are saved within the specified output directory (`dir_output`).
#' }
#'
#' It is important to note that Mahalanobis distance assumes multivariate normality and is sensitive to collinearity among variables.
#' While the covariance matrix accounts for correlations, it is strongly recommended that the `climate_variables` are not strongly correlated.
#' Consider performing a collinearity analysis beforehand, perhaps using the `vif_filter` function from this package.
#'
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor values ncell
#' @importFrom sf st_crs st_transform st_geometry st_as_sf st_make_grid
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave element_text
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile na.omit
#' @importFrom utils packageVersion
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#' set.seed(2458)
#' n_cells <- 100 * 100
#' r_clim_present <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
#' values(r_clim_present) <- c(
#'   (rowFromCell(r_clim_present, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
#'   (rowFromCell(r_clim_present, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
#'   (colFromCell(r_clim_present, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
#'   (colFromCell(r_clim_present, 1:n_cells) +
#'     (rowFromCell(r_clim_present, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'   (colFromCell(r_clim_present, 1:n_cells) /
#'     (rowFromCell(r_clim_present, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'   (colFromCell(r_clim_present, 1:n_cells) *
#'     (rowFromCell(r_clim_present, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4)),
#'   (colFromCell(r_clim_present, 1:n_cells) *
#'     (colFromCell(r_clim_present, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4)))
#' names(r_clim_present) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
#' terra::crs(r_clim_present) <- "EPSG:4326"
#' terra::plot(r_clim_present)
#' r_clim_present_filtered <- vif_filter(r_clim_present, th = 5)
#' hex_grid <- sf::st_sf(
#'   sf::st_make_grid(
#'     sf::st_as_sf(
#'       terra::as.polygons(
#'         terra::ext(r_clim_present)
#'       ),
#'     square = FALSE
#'   )
#' )
#' sf::st_crs(hex_grid) <- "EPSG:4326"
#' polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
#' polygons$name <- c("Pol_1", "Pol_2")
#' sf::st_crs(polygons) <- sf::st_crs(hex_grid)
#' study_area_polygon <- sf::st_as_sf(as.polygons(terra::ext(r_clim_present)))
#' sf::st_crs(study_area_polygon) <- "EPSG:4326"
#' terra::plot(r_clim_present[[1]])
#' terra::plot(polygons, add = TRUE, color = "transparent", lwd = 3)
#' terra::plot(study_area_polygon, add = TRUE, col = "transparent", lwd = 3, border = "red")
#' mh_rep(
#'   polygon = polygons,
#'   col_name = "name",
#'   climate_variables = r_clim_present_filtered,
#'   th = 0.95,
#'   dir_output = file.path(tempdir(), "ClimaRep"),
#'   save_raw = TRUE
#' )
#' }
#' @export
#'
mh_rep <- function(polygon,
                   col_name,
                   climate_variables,
                   th = 0.95,
                   dir_output = file.path(tempdir(), "ClimaRep"),
                   save_raw = FALSE) {
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))
  if (!inherits(polygon, "sf"))
    stop("Parameter 'polygon' must be an sf object.")
  if (!is.character(col_name) ||
      length(col_name) != 1 || !(col_name %in% names(polygon))) {
    stop("Parameter 'col_name' must be a single character string naming a column in 'polygon'.")
  }
  if (!inherits(climate_variables, "SpatRaster"))
    stop("Parameter 'climate_variables' must be a SpatRaster object.")
  if (!is.numeric(th) || length(th) != 1 || th < 0 || th > 1) {
    stop("Parameter 'th' must be a single numeric value between 0 and 1.")
  }
  if (!is.character(dir_output) || length(dir_output) != 1) {
    stop("Parameter 'dir_output' must be a single character string.")
  }
  if (terra::nlyr(climate_variables) < 2) {
    warning(
      "climate_variables has fewer than 2 layers. Mahalanobis distance is typically for multiple variables. Proceeding with single variable analysis if applicable."
    )
  }
  dir_rep <- file.path(dir_output, "Representativeness")
  dir_charts <- file.path(dir_output, "Charts")
  dirs_to_create <- c(dir_rep, dir_charts)
  if (save_raw) {
    dir_mh_raw <- file.path(dir_output, "Mh_Raw")
    dirs_to_create <- c(dirs_to_create, dir_mh_raw)
  }
  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  })
  message("Validating and adjusting Coordinate Reference Systems (CRS)...")
  reference_system_check <- terra::crs(climate_variables, describe = TRUE)$code
  reference_system <- terra::crs(climate_variables)
  if (sf::st_crs(polygon)$epsg != reference_system_check) {
    message("Adjusting CRS of polygon to match reference system.")
    polygon <- sf::st_transform(polygon, reference_system)
  }
  message("Starting per-polygon processing...")
  data_p <- na.omit(terra::as.data.frame(climate_variables, xy = TRUE))
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
    raster_polygon <- terra::mask(terra::crop(climate_variables, pol), pol)
    if (all(is.na(terra::values(raster_polygon)))) {
      warning("No available data for: ", pol_name, ". Skipping...")
      next
    }
    mu <- terra::global(raster_polygon, "mean", na.rm = TRUE)$mean
    calculate_mh <- function(data) {
      coords <- data[, 1:2]
      climate_data <- as.matrix(data[, 3:(ncol(data) - 1)])
      mh_values <- mahalanobis(climate_data, mu, cov_matrix)
      terra::rast(cbind(coords, mh_values),
                  type = "xyz",
                  crs = reference_system)
    }
    mh_present <- calculate_mh(data_p)
    if (save_raw) {
      terra::writeRaster(mh_present, file.path(dir_mh_raw, paste0("MH_rep_", pol_name, ".tif")), overwrite = TRUE)
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
                         paste0("TH_Rep_", pol_name, ".tif")
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
        paste0(pol_name, "_rep.jpeg")
      ),
      plot = p,
      width = 10,
      height = 10,
      dpi = 300
    )
  }
  message("All processes were completed")
  cat(paste("Output files in:", dir_output))
  return(invisible(NULL))
}
