#' @title Multivariate temporal climate tepresentativeness change analysis
#'
#' @description This function calculates Mahalanobis-based climate representativeness (or forward climate analogs) for input polygon across two time periods (present and future) within a defined area.
#'
#' The function categorizes cells based on how their climate representativeness changes, labeling them as Retained, Lost, or Novel.
#'
#' Representativeness is assessed by comparing the multivariate climate conditions of each cell, of the reference climate space (`present_climate_variables` and `future_climate_variables`), with the climate conditions within each specific input `polygon`.
#'
#' @param polygon An `sf` object containing the defined areas. **Its CRS will be used as the reference system.**
#' @param col_name `character`. Name of the column in the `polygon` object that contains unique identifiers for each polygon.
#' @param present_climate_variables A `SpatRaster` stack of climate variables representing current conditions.
#' @param future_climate_variables A `SpatRaster` stack containing the same climate variables as `present_climate_variables` but representing future projected conditions.
#' @param study_area A single `sf` polygon.
#' @param th `numeric` (0-1). Percentile threshold used to define representativeness. Cells with a Mahalanobis distance below or equal to the `th` are classified as representative (default: 0.95).
#' @param model `character`. Name or identifier of the climate model used (e.g., "MIROC6"). This parameter is used in output filenames and subdirectory names, allowing for better file management.
#' @param year `character`. Year or period of future climate data (e.g., "2070"). This parameter is used in output filenames and subdirectory names, allowing for better file management.
#' @param dir_output `character`. Path to the directory where output files will be saved. The function will create subdirectories within this path.
#' @param save_raw `logical`. If `TRUE`, saves the intermediate continuous Mahalanobis distance rasters calculated for each polygon before binary classification. The final binary classification rasters are always saved (default: `FALSE`).
#'
#' @return Writes the following outputs to disk within subdirectories of `dir_output`:
#' \itemize{
#'  \item Classification (`.tif` ) change rasters: Change category rasters (`0` for **Non-representative**, `1` for **Retained**, `2` for **Lost** and `3` for **Novel**) for each input polygon are saved in the `Change/` subdirectory.
#'  \item Visualization (`.jpeg`) maps: Image files visualizing the change classification results for each `polygon` are saved in the `Charts/` subdirectory.
#'  \item Raw Mahalanobis distance rasters: Optionally, they are saved as `.tif` files in the `Mh_Raw_Pre/` and `Mh_Raw_Fut/` subdirectories if `save_raw = TRUE`.
#' }
#'
#' @details
#' This function extends the approach used in `mh_rep` to assess Changes in Climate Representativeness (or forward climate analogs) over time.
#' While `mh_rep()` calculates representativeness in a single scenario, `mh_rep_ch()` adapts this by using the mean from the present polygon but a covariance matrix derived from the overall climate space across both present and future periods combined.
#'
#' Here are the key steps:
#' \enumerate{
#'  \item Checking CRS: Ensures that `future_climate_variables`, `climate_variables`, and `study_area` have matching CRSs with `polygon` by automatically transforming them if needed. The CRS of `polygon` will be used as the reference system.
#'  \item Crops and masks the `climate_variables` and `future_climate_variables` raster to the `study_area` to limit all subsequent calculations to the area of interest.
#'  \item Calculate the multivariate covariance matrix using climate data from all cells for both present and present-future time periods combined.
#'  \item For each polygon in the `polygon` object:
#'  \itemize{
#'    \item Crop and mask the present climate variables raster (`present_climate_variables`) to the boundary of the current polygon.
#'    \item Calculate the multivariate mean using the climate data from the previous step. This defines the climate centroid for the current polygon.
#'    Calculate the Mahalanobis distance for each cell relative to the centroid and the overall present and present-future covariance matrix.
#'    This results in a Mahalanobis distance raster for the present period and another for the future period.
#'    \item Apply the specified threshold (`th`) to Mahalanobis distances to determine which cells are considered representative. This threshold is a percentile of the Mahalanobis distances within the current polygon.
#'    \item Classify each cells, for both present and future periods, as Representative = `1` (Mahalanobis distance \eqn{\le} `th`) or Non-Representative = `0` (Mahalanobis distance $>$ `th`).
#'  }
#'  \item Compares the binary representativeness of each cell between the present and future periods and determines cells where conditions are:
#'  \itemize{
#'    \item `0`: **Non-representative**: Cells that are outside the defined Mahalanobis threshold in both present and future periods.
#'    \item `1`: **Retained**: Cells that are within the defined Mahalanobis threshold in both present and future periods. **Representative** if `Climarep::mh_rep()` is used
#'    \item `2`: **Lost**: Cells that are within the defined Mahalanobis threshold in the present period but outside it in the future period.
#'    \item `3`: **Novel**: Cells that are outside the defined Mahalanobis threshold in the present period but within it in the future period.
#'   }
#'  \item Saves the classification raster (`.tif`) and generates a corresponding visualization map (`.jpeg`) for each polygon. These are saved within the specified output directory (`dir_output`).
#'  All files are saved using the `model` and `year` parameters for better file management.
#' }
#'
#' It is important to note that Mahalanobis distance assumes is sensitive to collinearity among variables.
#' While the covariance matrix accounts for correlations, it is strongly recommended that the climate variables (`present_climate_variables`) are not strongly correlated.
#' Consider performing a collinearity analysis beforehand, perhaps using the `vif_filter()` function from this package.
#'
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor compareGeom resample
#' @importFrom sf st_crs st_transform st_geometry st_as_sf
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave element_text
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile
#' @importFrom utils packageVersion
#'
#' @examples
#' library(terra)
#' library(sf)
#' set.seed(2458)
#' n_cells <- 100 * 100
#' r_clim_present <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
#' values(r_clim_present) <- c(
#'    (terra::rowFromCell(r_clim_present, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
#'    (terra::rowFromCell(r_clim_present, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
#'    (terra::colFromCell(r_clim_present, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
#'    (terra::colFromCell(r_clim_present, 1:n_cells) +
#'      (terra::rowFromCell(r_clim_present, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'    (terra::colFromCell(r_clim_present, 1:n_cells) /
#'      (terra::rowFromCell(r_clim_present, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'    (terra::colFromCell(r_clim_present, 1:n_cells) *
#'      (terra::rowFromCell(r_clim_present, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4))),
#'    (terra::colFromCell(r_clim_present, 1:n_cells) *
#'      (terra::colFromCell(r_clim_present, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4)))
#' )
#' names(r_clim_present) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
#' terra::crs(r_clim_present) <- "EPSG:4326"
#'
#' vif_result <- ClimaRep::vif_filter(r_clim_present, th = 5)
#' print(vif_result$summary)
#' r_clim_present_filtered <- vif_result$filtered_raster
#' r_clim_future <- r_clim_present_filtered + 2
#' names(r_clim_future) <- names(r_clim_present_filtered)
#' hex_grid <- sf::st_sf(
#'    sf::st_make_grid(
#'      sf::st_as_sf(
#'        terra::as.polygons(
#'          terra::ext(r_clim_present_filtered))),
#'      square = FALSE))
#' sf::st_crs(hex_grid) <- "EPSG:4326"
#' polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
#' polygons$name <- c("Pol_1", "Pol_2")
#' study_area_polygon <- sf::st_as_sf(terra::as.polygons(terra::ext(r_clim_present_filtered)))
#' sf::st_crs(study_area_polygon) <- "EPSG:4326"
#' terra::plot(r_clim_present_filtered[[1]])
#' terra::plot(polygons, add = TRUE, color = "transparent", lwd = 3)
#' terra::plot(study_area_polygon, add = TRUE, col = "transparent", lwd = 3, border = "red")
#'
#' ClimaRep::mh_rep_ch(
#'    polygon = polygons,
#'    col_name = "name",
#'    present_climate_variables = r_clim_present_filtered,
#'    future_climate_variables = r_clim_future,
#'    study_area = study_area_polygon,
#'    th = 0.95,
#'    model = "ExampleModel",
#'    year = "2070",
#'    dir_output = file.path(tempdir(), "ClimaRepChange"),
#'    save_raw = TRUE)
#' @export
mh_rep_ch <- function(polygon,
                      col_name,
                      present_climate_variables,
                      future_climate_variables,
                      study_area,
                      th = 0.95,
                      model,
                      year,
                      dir_output = file.path(tempdir(), "ClimaRep"),
                      save_raw = FALSE) {
  if (!inherits(polygon, "sf"))
    stop("Parameter 'polygon' must be an sf object")
  if (!is.character(col_name) ||
      length(col_name) != 1 || !(col_name %in% names(polygon))) {
    stop("Parameter 'col_name' must be a single character string naming a column in 'polygon'")
  }
  if (!inherits(present_climate_variables, "SpatRaster"))
    stop("Parameter 'present_climate_variables' must be a SpatRaster object")
  if (!inherits(future_climate_variables, "SpatRaster"))
    stop("Parameter 'future_climate_variables' must be a SpatRaster object")
  if (!inherits(study_area, "sf"))
    stop("Parameter 'study_area' must be an sf object.")
  if (!is.numeric(th) || length(th) != 1 || th < 0 || th > 1) {
    stop("Parameter 'th' must be a single numeric value between 0 and 1")
  }
  if (!is.character(model) || length(model) != 1) {
    stop("Parameter 'model' must be a single character string")
  }
  if (!is.character(year) || length(year) != 1) {
    stop("Parameter 'year' must be a single character string")
  }
  if (!is.character(dir_output) || length(dir_output) != 1) {
    stop("Parameter 'dir_output' must be a single character string")
  }
  if (terra::nlyr(present_climate_variables) < 2) {
    warning("present_climate_variables has fewer than 2 layers. Mahalanobis distance is typically for multiple variables")
  }
  if (terra::nlyr(present_climate_variables) != terra::nlyr(future_climate_variables)) {
    stop("Number of layers in 'present_climate_variables' and 'future_climate_variables' must be the same")
  }

  message("Establishing output file structure")
  dir_present <- file.path(dir_output, "Mh_Raw_Pre")
  dir_future <- file.path(dir_output, "Mh_Raw_Fut")
  dir_change <- file.path(dir_output, "Change")
  dir_charts <- file.path(dir_output, "Charts")
  dirs_to_create <- c(dir_change, dir_charts)
  if (save_raw) {
    dirs_to_create <- c(dirs_to_create, dir_present, dir_future)
  }
  sapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  })
  message("Validating and adjusting Coordinate Reference Systems (CRS)")

  # Get the reference CRS from the polygon
  reference_system <- terra::crs(sf::st_crs(polygon)$wkt)
  reference_system_check <- sf::st_crs(polygon)$epsg

  if (is.na(reference_system_check) || reference_system_check == "") {
    stop("CRS for 'polygon' is undefined. Please set a valid CRS for the polygon")
  }
  if (terra::crs(present_climate_variables, describe = TRUE)$code != reference_system_check) {
    message("Adjusting CRS of present_climate_variables to match the polygon's system")
    present_climate_variables <- terra::project(present_climate_variables, reference_system)
  }
  if (terra::crs(future_climate_variables, describe = TRUE)$code != reference_system_check) {
    message("Adjusting CRS of future_climate_variables to match the polygon's system")
    future_climate_variables <- terra::project(future_climate_variables, reference_system)
  }
  if (sf::st_crs(study_area)$epsg != reference_system_check) {
    message("Adjusting CRS of study_area to match the polygon's system")
    study_area <- sf::st_transform(study_area, reference_system)
  }
  message("Defining the climate reference space")
  present_masked <- terra::mask(terra::crop(present_climate_variables, study_area),
                                study_area)
  future_masked <- terra::mask(terra::crop(future_climate_variables, study_area),
                               study_area)
  data_p_study <- na.omit(terra::as.data.frame(present_masked, xy = TRUE))
  data_f_study <- na.omit(terra::as.data.frame(future_masked, xy = TRUE))
  if (nrow(data_p_study) == 0 ||
      nrow(data_f_study) == 0) {
    stop("No valid climate data found within 'study_area' for one or both periods. Cannot calculate combined covariance matrix")
  }
  climate_data_cols <- 3:(ncol(data_p_study))
  data_combined_clim <- rbind(data_p_study[, climate_data_cols], data_f_study[, climate_data_cols])
  cov_matrix_pre <- suppressWarnings(cov(data_p_study[, climate_data_cols], use = "complete.obs"))
  cov_matrix_prefut <- suppressWarnings(cov(data_combined_clim, use = "complete.obs"))
  if (inherits(try(solve(cov_matrix_prefut), silent = TRUE)
               , "try-error")) {
    stop("Covariance matrix (combined present/future data) is singular (e.g., perfectly correlated or insufficient data). Consider filtering variables 'vif_filter()'")
  }
  message("Starting per-polygon processing:")
  classify_mh <- function(mh_raster, threshold) {
    terra::ifel(mh_raster <= threshold, 1, 0)
  }
  for (j in 1:nrow(polygon)) {
    pol <- polygon[j, ]
    pol_name <- as.character(pol[[col_name]])
    pol_name <- as.character(tolower(pol_name))
    pol_name <- iconv(pol_name, to = 'ASCII//TRANSLIT')
    pol_name <- gsub("[^a-z0-9_]+", "_", pol_name)
    pol_name <- gsub("__+", "_", pol_name)
    pol_name <- gsub("^_|_$", "", pol_name)
    message("\nProcessing polygon: ",pol_name, " (",j," of ",nrow(polygon),")")
    raster_polygon_present <- terra::mask(terra::crop(present_masked, pol), pol)
    if (all(is.na(terra::values(raster_polygon_present)))) {
      warning("No available data for: ",
              pol_name,
              " in the present period. Skipping.")
      next
    }
    mu_present_polygon <- terra::global(raster_polygon_present, "mean", na.rm = TRUE)$mean
    mh_values_present <- mahalanobis(as.matrix(data_p_study[, climate_data_cols]),
                                     mu_present_polygon,
                                     cov_matrix_pre)
    mh_values_future <- mahalanobis(as.matrix(data_f_study[, climate_data_cols]),
                                    mu_present_polygon,
                                    cov_matrix_prefut)
    mh_present <- terra::rast(cbind(data_p_study$x, data_p_study$y, mh_values_present),
                              type = "xyz",
                              crs = reference_system)
    mh_future <- terra::rast(cbind(data_f_study$x, data_f_study$y, mh_values_future),
                             type = "xyz",
                             crs = reference_system)
    if (save_raw) {
      terra::writeRaster(mh_present, file.path(dir_present, paste0("Mh_raw_Pre_", pol_name, ".tif")),
                         overwrite = TRUE)
      terra::writeRaster(mh_future, file.path(dir_future, paste0("Mh_raw_Fut_", model, "_", year, "_", pol_name, ".tif")),
                         overwrite = TRUE)
    }
    mh_poly <- terra::mask(mh_present, pol)
    th_value <- suppressWarnings(quantile(terra::values(mh_poly),
                                          probs = th,
                                          na.rm = TRUE))
    if (anyNA(th_value) || is.infinite(th_value)) {
      warning("No valid threshold was obtained for: ",
              pol_name,
              ". Skipping.")
      next
    }
    th_present <- classify_mh(mh_present, th_value)
    th_future <- classify_mh(mh_future, th_value)
    raster_final <- (th_present * th_future) +
      ((th_present - (th_present * th_future)) * 2) +
      ((th_future - (th_present * th_future)) * 3)
    terra::writeRaster(raster_final,
                       file.path(
                         dir_output,
                         "Change",
                         paste0("Th_change_", model, "_", year, "_", pol_name, ".tif")),
                       overwrite = TRUE)
    raster_final_factor <- terra::as.factor(raster_final)
    p <- suppressMessages(
      ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = raster_final_factor) +
        ggplot2::geom_sf(
          data = study_area,
          color = "gray50",
          fill = NA,
          linewidth = 1) +
        ggplot2::geom_sf(
          data = pol,
          color = "black",
          fill = NA) +
        ggplot2::scale_fill_manual(
          name = " ",
          values = c(
            "0" = "grey90", # Non-representative
            "1" = "aquamarine4", # Retained
            "2" = "coral1", # Lost
            "3" = "steelblue2"), # Novel
          labels = c(
            "0" = "Non-representative",
            "1" = "Retained",
            "2" = "Lost",
            "3" = "Novel"),
          na.value = "transparent",
          na.translate = FALSE,
          drop = FALSE) +
        ggplot2::ggtitle(pol_name) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = element_text(hjust = 0.5)))
    ggplot2::ggsave(
      filename = file.path(
        dir_output,
        "Charts",
        paste0(pol_name, "_rep_change.jpeg")
      ),
      plot = p,
      width = 10,
      height = 10,
      dpi = 300)
  }
  message("All processes were completed")
  message(paste("Output files in: ", dir_output))
  return(invisible(NULL))
}
