#' @title Multivariate Climate Representativeness Analysis
#'
#' @description This function calculates Mahalanobis-based Climate Representativeness for input polygon within a defined area.
#'
#' Representativeness is assessed by comparing the multivariate climate conditions of each cell, of the reference climate space (`climate_variables`), with the climate conditions within each specific input `polygon`.
#'
#' @param polygon An `sf` object containing the defined areas. **Must have the same CRS as** `climate_variables`.
#' @param col_name `character`. Name of the column in the `polygon` object that contains unique identifiers for each polygon.
#' @param climate_variables A `SpatRaster` stack of climate variables representing the conditions. Its CRS will be used as the reference system.
#' @param th `numeric` (0-1). Percentile threshold used to define representativeness. Cells with a Mahalanobis distance below or equal to the `th` are classified as representative (default: 0.95).
#' @param dir_output `character`. Path to the directory where output files will be saved. The function will create subdirectories within this path.
#' @param save_raw `logical`. If `TRUE`, saves the intermediate continuous Mahalanobis distance rasters calculated for each polygon before binary classification. The final binary classification rasters are always saved (default: `FALSE`).
#'
#' @return Writes the following outputs to disk within subdirectories of `dir_output`:
#' \itemize{
#'  \item Classification (`.tif` ) rasters: Binary rasters (`0` for **Non-representative** and`1` for **Representative**) for each input polygon are saved in the `Representativeness/` subdirectory.
#'  \item Visualization (`.jpeg`) maps: Image files visualizing the classification results for each `polygon` are saved in the `Charts/` subdirectory.
#'  \item Raw Mahalanobis distance rasters: Optionally saved as `.tif` files in the `Mh_Raw/` subdirectory if `save_raw = TRUE`.
#' }
#'
#' @details
#' This function performs a multivariate analysis using Mahalanobis distance to assess the Climate Representativeness of input `polygons` for a single time period.
#'
#' Crucially, this function assumes that all spatial inputs (`polygon`, `climate_variables`) are already correctly aligned and share the same Coordinate Reference System (CRS). If inputs do not meet these criteria, the function will stop with an informative error.
#'
#' Here are the key steps:
#' \enumerate{
#'  \item Checking of spatial inputs: Ensures that `polygon` and `climate_variables` have matching CRSs.
#'  \item Calculate the multivariate covariance matrix using climate data from all cells.
#'  \item For each polygon in the `polygon` object:
#'  \itemize{
#'    \item Crop and mask the climate variables raster (`climate_variables`) to the boundary of the current polygon.
#'    \item Calculate the multivariate mean using the climate data from the previous step. This defines the climate centroid for the current polygon.
#'    \item Calculate the Mahalanobis distance for each cell relative to the centroid and covariance matrix.
#'    \item Apply the specified threshold (`th`) to Mahalanobis distances to determine which cells are considered representative. This threshold is a percentile of the Mahalanobis distances within the current polygon.
#'    \item Classify each cell as Representative = `1` (Mahalanobis distance \eqn{\le} `th`) or Non-Representative = `0` (Mahalanobis distance $>$ `th`).
#'  }
#'  \item Saves the binary classification raster (`.tif`) and generates a corresponding visualization map (`.jpeg`) for each polygon. These are saved within the specified output directory (`dir_output`).
#' }
#'
#' It is important to note that Mahalanobis distance is sensitive to collinearity among variables.
#' While the covariance matrix accounts for correlations, it is strongly recommended that the `climate_variables` are not strongly correlated.
#' Consider performing a collinearity analysis beforehand, perhaps using the `vif_filter` function from this package.
#'
#' @importFrom terra crs project crop mask global as.data.frame rast writeRaster as.factor values ncell compareGeom
#' @importFrom sf st_crs st_transform st_geometry st_as_sf st_make_grid
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual ggtitle theme_minimal ggsave element_text
#' @importFrom tidyterra geom_spatraster
#' @importFrom stats mahalanobis cov quantile na.omit
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
#' hex_grid <- sf::st_sf(
#'    sf::st_make_grid(
#'      sf::st_as_sf(
#'        terra::as.polygons(
#'          terra::ext(r_clim_present_filtered))),
#'      square = FALSE)
#'      )
#' sf::st_crs(hex_grid) <- "EPSG:4326"
#' polygons <- hex_grid[sample(nrow(hex_grid), 2), ]
#' polygons$name <- c("Pol_A", "Pol_B")
#' terra::plot(r_clim_present_filtered[[1]])
#' terra::plot(polygons, add = TRUE, color = "transparent", lwd = 3)
#'
#' ClimaRep::mh_rep(
#'    polygon = polygons,
#'    col_name = "name",
#'    climate_variables = r_clim_present_filtered,
#'    th = 0.95,
#'    dir_output = file.path(tempdir(), "ClimaRep"),
#'    save_raw = TRUE
#'    )
#' @export
mh_rep <- function(polygon,
                   col_name,
                   climate_variables,
                   th = 0.95,
                   dir_output = file.path(tempdir(), "ClimaRep"),
                   save_raw = FALSE) {
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
      "Climate_variables has fewer than 2 layers. Mahalanobis distance is typically for multiple variables."
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
  message("Validating and adjusting Coordinate Reference Systems (CRS)")
  reference_system_check <- terra::crs(climate_variables, describe = TRUE)$code
  reference_system <- terra::crs(climate_variables)
  if (sf::st_crs(polygon)$epsg != reference_system_check) {
    message("Adjusting CRS of polygon to match reference system.")
    polygon <- sf::st_transform(polygon, reference_system)
  }
  message("Starting per-polygon processing:")
  data_p <- na.omit(terra::as.data.frame(climate_variables, xy = TRUE))
  if (ncol(data_p) < 3) {
    stop("Not enough variables in 'climate_variables' to calculate Mahalanobis distance.")
  }
  climate_data_cols <- 3:(ncol(data_p) - 0)
  cov_matrix <- suppressWarnings(cov(data_p[, climate_data_cols], use = "complete.obs"))
  if (inherits(try(solve(cov_matrix), silent = TRUE)
               , "try-error")) {
    stop(
      "Covariance matrix is singular (e.g., perfectly correlated or insufficient data). Consider filtering variables 'vif_filter()'."
    )
  }
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
      warning("No available data for: ", pol_name, ". Skipping.")
      next
    }
    mu <- terra::global(raster_polygon, "mean", na.rm = TRUE)$mean
    mh_values <- mahalanobis(as.matrix(data_p[, climate_data_cols]), mu, cov_matrix)
    mh_present <- terra::rast(cbind(data_p$x, data_p$y, mh_values),
                              type = "xyz",
                              crs = reference_system)
    if (save_raw) {
      terra::writeRaster(mh_present, file.path(dir_mh_raw, paste0("Mh_raw_", pol_name, ".tif")), overwrite = TRUE)
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
    classify_mh <- function(mh_raster, threshold) {
      terra::ifel(mh_raster <= threshold, 1, 0)
    }
    th_present <- classify_mh(mh_present, th_value)
    raster_final <- th_present
    terra::writeRaster(raster_final,
                       file.path(
                         dir_output,
                         "Representativeness",
                         paste0("Th_rep_", pol_name, ".tif")
                       ),
                       overwrite = TRUE)
    raster_final_factor <- terra::as.factor(raster_final)
    p <- suppressMessages(
      ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = raster_final_factor) +
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
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = element_text(hjust = 0.5))
    )
    ggplot2::ggsave(
      filename = file.path(dir_output, "Charts", paste0(pol_name, "_rep.jpeg")),
      plot = p,
      width = 10,
      height = 10,
      dpi = 300
    )
  }
  message("All processes were completed")
  message(paste("Output files in: ", dir_output))
  return(invisible(NULL))
}
