#' @title Multivariate Temporal Climate Representativeness Change Analysis
#'
#' @description Calculates Mahalanobis-based climate representativeness or forward climate analogs for input `polygon` across two time periods (present and future) within a defined study area. The function identifies areas of climate representativeness **Retained**, **Lost**, or **Novel**.
#' Representativeness is assessed by comparing the multivariate climate conditions of each cell to the reference climate space defined by the climate conditions within each specific input polygon in the **present** period, scaled by the overall climate variability across both time periods within the `study_area`.
#'
#' @param polygon An `sf` object containing the defined areas. **Must have the same CRS as `present_climate_variables`.**
#' @param col_name `Character`. Name of the column in the `polygon` object that contains unique identifiers for each polygon.
#' @param present_climate_variables `SpatRaster`. A raster stack of climate variables representing **present** conditions (preferably not strongly correlated). Its CRS will be used as the reference for all other spatial inputs.
#' @param future_climate_variables `SpatRaster`. A raster stack containing the *same* climate variables as `present_climate_variables` but representing **future** projected conditions. **Must have the same CRS, extent, and resolution as `present_climate_variables`.**
#' @param study_area `sf`. A single polygon defining the overall boundary of the climate space. **Must have the same CRS as `present_climate_variables`.**
#' @param th `Numeric` (0-1). Percentile threshold used to define representativeness. Cells with a Mahalanobis distance below or equal to the `th` are classified as representative (default: 0.95).
#' @param model `Character.` Name or identifier of the climate model used for future projections (e.g., "MIROC6"). This parameter is mandatory and used in output filenames and subdirectory names.
#' @param year `Character.` Projection year or period for future climate data (e.g., "2070"). This parameter is mandatory and used in output filenames and subdirectory names.
#' @param dir_output `Character`. Path to the directory where output files will be saved. The function will create subdirectories within this path.
#' @param save_raw `Logical.` If `TRUE`, saves the intermediate continuous Mahalanobis distance rasters calculated for each polygon before binary classification. The final binary classification rasters are always saved (default: FALSE).
#'
#' @return Invisibly returns NULL. Writes the following outputs to disk within subdirectories of `dir_output`:
#' \itemize{
#'  \item Classification (`.tif` ) change raster: Change category rasters (`0` for **Non-representative**, `1` for **Retained**, `2` for **Lost** and `3` for **Novel**) for each input `polygon` are saved in the `Change/` subdirectory.
#'  \item Visualization (`.jpeg`) maps: Image files visualizing the change classification results for each `polygon` are saved in the `Charts/` subdirectory.
#'  \item Raw Mahalanobis Distance Rasters: *Opcionalmente* se guardan como archivos `.tif` en los subdirectorios `Mh_Raw_Pre/` y `Mh_Raw_Fut/` si `save_raw = TRUE`.
#' }
#'
#' @details
#' This function extends the multivariate analysis approach used in `mh_rep` to assess changes in climate representativeness (or **forward climate analogs**) over time.
#' While `mh_rep()` calculates representativeness relative to the mean and covariance of climate conditions *without change*, `mh_rep_ch()` adapts this by using the mean from the present polygon but a covariance matrix derived from the overall climate space across both present and future periods combined.
#'
#' **Crucially, this function assumes that all spatial inputs (`polygon`, `present_climate_variables`, `future_climate_variables`, `study_area`) are already correctly aligned and share the same Coordinate Reference System (CRS) and consistent spatial properties (extent, resolution). If inputs do not meet these criteria, the function will stop with an informative error.**
#'
#' Key workflow steps include:
#' \enumerate{
#'  \item **Pre-check of spatial inputs**: Ensures that `present_climate_variables`, `future_climate_variables`, `polygon`, and `study_area` all have matching CRSs, and that `present_climate_variables` and `future_climate_variables` share identical extents and resolutions.
#'  \item Calculate the multivariate covariance matrix using climate data from *all cells within the study area* for *both* present and future time periods combined. This captures the overall climate variability against which distances are scaled.
#'  \item For each polygon in the `polygon` object:
#'  \itemize{
#'    \item Calculate the multivariate mean climate using climate data *only from within that polygon* for the **present** time period. This defines the center of the reference climate space for the polygon.
#'    \item Calculate the Mahalanobis distance for all cells within the study area, relative to the polygon's **present** conditions (calculated in the previous step) and the overall **present** and **future** covariance matrix (calculated in step 2).
#'    This results in a Mahalanobis distance raster for the present period and a Mahalanobis distance raster for the future period. The generated Mahalanobis rasters are resampled to ensure perfect alignment with the input climate rasters.
#'    \item Determine a threshold (`th`) based on the `th` percentile of the Mahalanobis distances calculated *only from within the present polygon*.
#'    \item Classify, based on the threshold (`th`), all cells within the `study_area` for both present and future periods as Representative or forward climate analogs (Mahalanobis distance \eqn{\le} `th`) or "Non-Representative" (Mahalanobis distance $>$ `th`).
#'  }
#'  \item Compares the binary representativeness status of each cell between the present and future periods and determines cells where conditions are:
#'  \itemize{
#'    \item `0`: **Non-represented**: Cells that are outside the defined Mahalanobis threshold in both present and future periods.
#'    \item `1`: **Retained**: Cells that are within the defined Mahalanobis threshold in both present and future periods.
#'    \item `2`: **Lost**: Cells that are within the defined Mahalanobis threshold in the present period but outside it in the future period.
#'    \item `3`: **Novel**: Cells that are outside the defined Mahalanobis threshold in the present period but within it in the future period.
#'  }
#'  \item Output Generation: Saves the classification raster (`.tif`) for each polygon and generates a corresponding visualization map (`.jpeg`). These are saved within the specified output directory (`dir_output`).
#'  All files are saved using the `model` and `year` parameters in filenames for identification.
#' }
#'
#' It is important to note that Mahalanobis distance assumes multivariate normality and is sensitive to collinearity among variables.
#' While the covariance matrix accounts for correlations, it is strongly recommended that the climate variables (`present_climate_variables` and `future_climate_variables`) are not strongly correlated.
#' Consider performing a collinearity analysis beforehand, perhaps using the `vif_filter` function from this package.
#'
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor compareGeom resample
#' @importFrom sf st_crs st_transform st_geometry st_as_sf
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave element_text
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile
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
#' # r_clim_present_filtered <- vif_filter(r_clim_present, th = 5) # Assuming vif_filter is available
#' r_clim_present_filtered <- r_clim_present[[1:5]] # Example if vif_filter not available
#' r_clim_future <- r_clim_present_filtered + 2
#' names(r_clim_future) <- names(r_clim_present_filtered)
#' hex_grid <- sf::st_sf(
#'    sf::st_make_grid(
#'      sf::st_as_sf(
#'        terra::as.polygons(
#'          terra::ext(r_clim_present))),
#'      square = FALSE))
#' sf::st_crs(hex_grid) <- "EPSG:4326"
#' polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
#' polygons$name <- c("Pol_1", "Pol_2")
#' study_area_polygon <- sf::st_as_sf(terra::as.polygons(terra::ext(r_clim_present)))
#' sf::st_crs(study_area_polygon) <- "EPSG:4326"
#' terra::plot(r_clim_present[[1]])
#' terra::plot(polygons, add = TRUE, color = "transparent", lwd = 3)
#' terra::plot(study_area_polygon, add = TRUE, col = "transparent", lwd = 3, border = "red")
#'
#' mh_rep_ch(
#'    polygon = polygons,
#'    col_name = "name",
#'    present_climate_variables = r_clim_present_filtered,
#'    future_climate_variables = r_clim_future,
#'    study_area = study_area_polygon,
#'    th = 0.95,
#'    model = "ExampleModel",
#'    year = "2070",
#'    dir_output = file.path(tempdir(), "ClimaRepOutput"),
#'    save_raw = TRUE)
#' }
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
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn))
  if (!inherits(polygon, "sf"))
    stop("Parameter 'polygon' must be an sf object.")
  if (!is.character(col_name) ||
      length(col_name) != 1 || !(col_name %in% names(polygon))) {
    stop("Parameter 'col_name' must be a single character string naming a column in 'polygon'.")
  }
  if (!inherits(present_climate_variables, "SpatRaster"))
    stop("Parameter 'present_climate_variables' must be a SpatRaster object.")
  if (!inherits(future_climate_variables, "SpatRaster"))
    stop("Parameter 'future_climate_variables' must be a SpatRaster object.")
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
  if (terra::nlyr(present_climate_variables) < 2) {
    warning(
      "present_climate_variables has fewer than 2 layers. Mahalanobis distance is typically for multiple variables. Proceeding with single variable analysis if applicable."
    )
  }
  if (terra::nlyr(present_climate_variables) != terra::nlyr(future_climate_variables)) {
    stop(
      "Number of layers in 'present_climate_variables' and 'future_climate_variables' must be the same."
    )
  }
  message("Validating CRS and spatial alignment of inputs...")
  ref_crs <- terra::crs(present_climate_variables, describe = TRUE)$code
  if (is.na(ref_crs) || ref_crs == "") {
    stop("CRS for 'present_climate_variables' is undefined. Please set a valid CRS.")
  }
  if (terra::crs(future_climate_variables, describe = TRUE)$code != ref_crs) {
    stop(
      "CRS mismatch: 'future_climate_variables' must have the same CRS as 'present_climate_variables'."
    )
  }
  if (sf::st_crs(polygon)$epsg != ref_crs) {
    stop("CRS mismatch: 'polygon' must have the same CRS as 'present_climate_variables'.")
  }
  if (sf::st_crs(study_area)$epsg != ref_crs) {
    stop("CRS mismatch: 'study_area' must have the same CRS as 'present_climate_variables'.")
  }

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
  message("Starting process")
  present_study_area_masked <- terra::mask(terra::crop(present_climate_variables, study_area),
                                           study_area)
  future_study_area_masked <- terra::mask(terra::crop(future_climate_variables, study_area),
                                          study_area)
  data_p_study <- na.omit(terra::as.data.frame(present_study_area_masked, xy = TRUE))
  data_f_study <- na.omit(terra::as.data.frame(future_study_area_masked, xy = TRUE))
  if (nrow(data_p_study) == 0 || nrow(data_f_study) == 0) {
    stop(
      "No valid climate data found within 'study_area' for one or both periods. Cannot calculate combined covariance matrix."
    )
  }
  climate_var_names <- names(present_climate_variables)
  data_p_study_clim <- data_p_study[, climate_var_names, drop = FALSE]
  data_f_study_clim <- data_f_study[, climate_var_names, drop = FALSE]
  if (ncol(data_p_study_clim) < 2) {
    stop(
      "Not enough climate variables (layers) to calculate a multivariate Mahalanobis distance. Need at least 2 layers."
    )
  }
  data_combined_clim <- rbind(data_p_study_clim, data_f_study_clim)
  cov_matrix <- cov(data_combined_clim, use = "complete.obs")
  if (inherits(try(solve(cov_matrix), silent = TRUE)
               , "try-error")) {
    stop(
      "Covariance matrix (calculated from combined present/future study area data) is singular. This can happen if variables are perfectly correlated or there's insufficient unique data points. Consider filtering collinear variables or checking data quality within the study area."
    )
  }
  full_data_present <- na.omit(terra::as.data.frame(present_climate_variables, xy = TRUE))
  full_data_future <- na.omit(terra::as.data.frame(future_climate_variables, xy = TRUE))
  full_data_present_clim <- full_data_present[, climate_var_names, drop = FALSE]
  full_data_future_clim <- full_data_future[, climate_var_names, drop = FALSE]
  full_data_present_coords <- full_data_present[, c("x", "y")]
  full_data_future_coords <- full_data_future[, c("x", "y")]
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
    raster_polygon_present <- terra::mask(terra::crop(present_climate_variables, pol), pol)
    if (all(is.na(terra::values(raster_polygon_present)))) {
      warning("No available data for: ",
              pol_name,
              " in the present period. Skipping...")
      next
    }
    mu_present_polygon <- terra::global(raster_polygon_present, "mean", na.rm = TRUE)$mean
    mh_values_present <- mahalanobis(as.matrix(full_data_present_clim),
                                     mu_present_polygon,
                                     cov_matrix)
    mh_values_future <- mahalanobis(as.matrix(full_data_future_clim),
                                    mu_present_polygon,
                                    cov_matrix)
    mh_present_full <- terra::rast(
      cbind(full_data_present_coords, mh_values_present),
      type = "xyz",
      crs = terra::crs(present_climate_variables)
    )
    mh_future_full <- terra::rast(
      cbind(full_data_future_coords, mh_values_future),
      type = "xyz",
      crs = terra::crs(future_climate_variables)
    )
    mh_present <- terra::mask(mh_present_full, study_area)
    mh_future <- terra::mask(mh_future_full, study_area)
    if (save_raw) {
      terra::writeRaster(mh_present, file.path(
        dir_present,
        paste0("MahalanobisRaw_Pre_", pol_name, ".tif")
      ), overwrite = TRUE)
      terra::writeRaster(mh_future, file.path(
        dir_future,
        paste0(
          "MahalanobisRaw_Fut_",
          model,
          "_",
          year,
          "_",
          pol_name,
          ".tif"
        )
      ), overwrite = TRUE)
    }
    mh_poly_for_th <- terra::mask(mh_present_full, pol)
    th_value <- quantile(terra::values(mh_poly_for_th),
                         probs = th,
                         na.rm = TRUE)
    if (anyNA(th_value) || is.infinite(th_value)) {
      warning("No valid threshold was obtained for: ",
              pol_name,
              ". Skipping...")
      next
    }
    classify_mh <- function(mh_raster, threshold) {
      terra::ifel(mh_raster <= threshold, 1, 0)
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
                         paste0("TH_change_", model, "_", year, "_", pol_name, ".tif")
                       ),
                       overwrite = TRUE)
    raster_final_factor <- terra::as.factor(raster_final)
    p <- suppressMessages(
      ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = raster_final_factor) +
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
            # Non-represented
            "1" = "aquamarine4",
            # Retained
            "2" = "coral1",
            # Lost
            "3" = "aquamarine2"     # Novel
          ),
          labels = c(
            "0" = "Non-represented",
            "1" = "Retained",
            "2" = "Lost",
            "3" = "Novel"
          ),
          na.value = "transparent",
          na.translate = FALSE,
          drop = FALSE
        ) +
        ggplot2::ggtitle(pol_name) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = element_text(hjust = 0.5))
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
  cat(paste("Output files in: ", dir_output, "\n"))
  return(invisible(NULL))
}
