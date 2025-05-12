#' @title Multivariate Temporal Climate Representativeness Change Analysis
#'
#' @description Calculates Mahalanobis-based climatic representativeness for input polygons across two time periods (present and future) within a defined study area. The function identifies areas of climate representativeness persistence, loss, or gain.
#' Representativeness is assessed by comparing the multivariate climate conditions of each cell to the reference climate space defined by the climate conditions *within* each specific input polygon in the *present* period, scaled by the overall climate variability across both time periods within the study area.
#'
#' @param polygon An sf object containing the analysis regions (input polygons).
#' @param col_name Character. Name of the column in the `polygon` object that contains unique identifiers for each input polygon.
#' @param present_climatic_variables SpatRaster. A raster stack of climate variables representing present conditions (preferably not strongly correlated).
#' @param future_climatic_variables SpatRaster. A raster stack containing the *same* climate variables as `present_climatic_variables` but representing future projected conditions. Must have the same extent and resolution as `present_climatic_variables`.
#' @param study_area sf object. A single polygon defining the overall boundary of the study area. Used to define the extent for overall climate variability calculation and output rasters.
#' @param th Numeric (0-1). Quantile threshold used to define representativeness based on Mahalanobis distance. Cells with a Mahalanobis distance below or equal to the distance corresponding to this quantile (calculated from distances *within* the present polygon) are classified as representative (default: 0.9).
#' @param model Character. Name or identifier of the climate model used for future projections (e.g., "MIROC6"). This parameter is mandatory and used in output filenames and subdirectory names.
#' @param year Character. Projection year or period for future climate data (e.g., "2070"). This parameter is mandatory and used in output filenames and subdirectory names.
#' @param dir_output Character. Path to the directory where output files will be saved. The function will create subdirectories within this path (default: "output/").
#' @param save_raw Logical. If TRUE, saves intermediate continuous Mahalanobis distance rasters for both present and future periods within the study area extent. The final classification rasters listed in `@return` are saved regardless of this parameter setting (default: TRUE).
#'
#' @return Invisibly returns NULL. Writes the following outputs to disk within subdirectories of `dir_output`:
#' \itemize{
#'   \item Classification Change Raster: A GeoTIFF raster showing the change category for each input polygon (e.g., codes 0, 1, 2, 3) is saved in the `Change/` subdirectory.
#'   The categories typically represent: 0 - Not Representative in either period; 1 - Persistence (Representative in both); 2 - Loss (Representative Present, Not Future); 3 - Gain (Not Representative Present, Representative Future).
#'   \item Visualization Maps: JPEG (or PNG) image files visualizing the change classification results for each input polygon are saved in the `Charts/` subdirectory.
#'   \item Intermediate Raw Mahalanobis Distance Rasters: *Optionally* saved as GeoTIFF files if `save_raw = TRUE`. The raw present distances (relative to the present polygon mean, scaled by overall variance) are saved in `Mh_Raw_Pre/`.
#'   The raw future distances (relative to the present polygon mean, scaled by overall variance) are saved in `Mh_Raw_Fut/`.
#' }
#'
#' @details
#' This function performs a multivariate analysis using Mahalanobis distance to assess and compare the climatic representativeness of input polygons across present and future periods within a defined study area.
#' Representativeness for both periods is determined relative to the multivariate climate space defined by the conditions within each input polygon in the present time period.
#' The calculation uses the mean climate of the present polygon, but accounts for the overall climate variability (covariance) observed across the entire study area for both present and future time periods combined.
#'
#' Key workflow steps include:
#' \enumerate{
#'   \item **CRS Harmonization:** Ensures all spatial inputs (`polygon`, `present_climatic_variables`, `future_climatic_variables`, `study_area`) share the same Coordinate Reference System (CRS), using the CRS of `present_climatic_variables` as the reference. Inputs are reprojected if necessary.
#'   \item **Reference Variability Calculation:** Calculates the multivariate covariance matrix using climate data from *all cells within the study area* for *both* present and future time periods combined. This captures the overall climate variability against which distances are scaled.
#'   \item **Per-Polygon Processing:** For each input polygon:
#'   \itemize{
#'     \item Calculates the multivariate mean climate using climate data *only from within that polygon* for the **present** time period. This defines the center of the reference climate space for the polygon.
#'     \item Calculates the Mahalanobis distance for *all cells within the study area*, relative to the polygon's **present** mean (calculated in the previous step) and the **overall** present+future covariance matrix (calculated in step 2).
#'     This results in a Mahalanobis distance raster for the present period and a Mahalanobis distance raster for the future period.
#'     \item Determines a threshold (`th_value`) based on the `th` quantile of the Mahalanobis distances calculated *only from within the present polygon*.
#'     \item Classifies all cells within the study area for both present and future periods as "Representative" (Mahalanobis distance < threshold) or "Non-Representative" (distance > threshold).
#'   }
#'   \item **Classification of Change:** Compares the binary representativeness status of each cell between the present and future periods to classify it into one of four categories:
#'   \itemize{
#'     \item Not Representative in either period (e.g., code 0)
#'     \item Persistence (Representative in both present and future periods; e.g., code 1)
#'     \item Loss (Representative in the present period, but not in the future; e.g., code 2)
#'     \item Gain (Not representative in the present period, but representative in the future; e.g., code 3)
#'   \item **Output Generation:** Saves the resulting categorical change raster for each input polygon as a GeoTIFF file and generates a corresponding visualization map (JPEG/PNG). Intermediate raw distance rasters are optionally saved.
#'   All files are saved within the specified output directory structure, using the `model` and `year` parameters in filenames for identification.
#' }}
#'
#' It is important to note that Mahalanobis distance assumes multivariate normality and is sensitive to collinearity among variables.
#' While the covariance matrix accounts for correlations, it is strongly recommended that the climate variables (`present_climatic_variables` and `future_climatic_variables`) are not strongly correlated.
#' Consider performing a collinearity analysis (e.g., using Variance Inflation Factor - VIF) beforehand, perhaps using the `vif_filter` function from this package.
#'
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor
#' @importFrom sf st_crs st_transform st_geometry
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile cor
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
#'     (colFromCell(r_clim_present, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4))
#' )
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
#'   climatic_variables = r_clim_present_filtered,
#'   th = 0.95,
#'   dir_output = file.path(tempdir(), "ClimaRep"),
#'   save_raw = TRUE
#' )
#' }
#' @export
#'
mh_rep_ch <- function(polygon,
                      col_name,
                      present_climatic_variables,
                      future_climatic_variables,
                      study_area,
                      th = 0.95,
                      model,
                      year,
                      dir_output = file.path(tempdir(), "ClimaRep"),
                      save_raw = TRUE) {
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))
  if (!inherits(polygon, "sf"))
    stop("Parameter 'polygon' must be an sf object.")
  if (!is.character(col_name) ||
      length(col_name) != 1 || !(col_name %in% names(polygon))) {
    stop("Parameter 'col_name' must be a single character string naming a column in 'polygon'.")
  }
  if (!inherits(present_climatic_variables, "SpatRaster"))
    stop("Parameter 'present_climatic_variables' must be a SpatRaster object.")
  if (!inherits(future_climatic_variables, "SpatRaster"))
    stop("Parameter 'future_climatic_variables' must be a SpatRaster object.")
  if (!inherits(study_area, "sf"))
    stop("Parameter 'study_area' must be an sf object.")
  if (!is.numeric(th) || length(th) != 1 || th < 0 || th > 1) {
    stop("Parameter 'th' must be a single numeric value between 0 and 1.")
  }
  if (!is.character(model) || length(model) != 1) {
    stop("Parameter 'model' must be a single character string.")
  }
  if (!is.character(year) || length(year) != 1) {
    stop("Parameter 'year' must be a single character string.")
  }
  if (!is.character(dir_output) || length(dir_output) != 1) {
    stop("Parameter 'dir_output' must be a single character string.")
  }
  if (terra::nlyr(present_climatic_variables) < 2) {
    warning(
      "climatic_variables has fewer than 2 layers. Mahalanobis distance is typically for multiple variables. Proceeding with single variable analysis if applicable."
    )
  }
  dir_present <- file.path(dir_output, "Mh_Raw_Pre")
  dir_future <- file.path(dir_output, "Mh_Raw_Fut")
  dir_change <- file.path(dir_output, "Change")
  dir_charts <- file.path(dir_output, "Charts")
  dirs_to_create <- c(dir_present, dir_future, dir_charts, dir_change)
  if (save_raw) {
    dirs_to_create <- c(dir_present, dir_future, dir_change, dir_charts)
  } else {
    dirs_to_create <- c(dir_change, dir_charts)
  }
  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  })
  message("Validating and adjusting Coordinate Reference Systems (CRS)...")
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
  cov_matrix <- cov(data_p_f[, 3:(ncol(data_p_f) - 1)], use = "complete.obs")
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
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
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
    mh_future <- calculate_mh(data_f)
    if (save_raw) {
      terra::writeRaster(
        mh_present,
        paste0(dir_present, "/MahalanobisRaw_", pol_name, ".tif"),
        overwrite = TRUE
      )
      terra::writeRaster(
        mh_future,
        paste0(
          dir_future,
          "/MahalanobisRaw_post_",
          model,
          "_",
          year,
          "_",
          pol_name,
          ".tif"
        ),
        overwrite = TRUE
      )
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
    th_future <- classify_mh(mh_future, th_value)
    change <- th_present * th_future
    solo_presente <- th_present - change
    solo_futura <- th_future - change
    raster_final <- change + (solo_presente * 2) + (solo_futura * 3)
    terra::writeRaster(raster_final,
                       file.path(
                         dir_output,
                         "change",
                         paste0("TH_change_", model, "_", year, "_", pol_name, ".tif")
                       ),
                       overwrite = TRUE)

    raster_final <- terra::as.factor(raster_final)
    p <- suppressMessages(
      ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = raster_final) +
        ggplot2::geom_sf(
          data = study_area,
          color = "gray50",
          fill = NA,
          linewidth = 1
        ) +
        ggplot2::geom_sf(
          data = pol,
          color = "black",
          fill = NA
        ) +
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
        ggplot2::theme_minimal()
    )
    ggplot2::ggsave(
      filename = file.path(
        dir_output,
        "Charts",
        paste0(pol_name, "_rep_change.jpeg")
      ),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  message("All processes were completed")
  cat(paste("Output files in: ", dir_output))
  return(invisible(NULL))
}
