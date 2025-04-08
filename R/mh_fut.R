## FUTURE----
#' Function to process the data of the future series
#'
#' This function processes the data for a future time series (defined by a
#' specific model and year) to calculate Mahalanobis distance and a
#' thresholded raster.
#'
#' @param j An index or identifier for the variable being processed.
#' @param th Threshold probability for defining the presence/absence based on
#'   the Mahalanobis distance. Default is 0.95.
#'
#' @return This function does not explicitly return an object but creates and
#'   assigns several raster objects to the global environment and writes them
#'   to disk. These include:
#'   \itemize{
#'     \item \code{MH\_\{model\}\_\{year\}\_\{name[j]\}.tif}: A raster layer of
#'       Mahalanobis distances for the future scenario.
#'     \item \code{TH\_MH\_\{model\}\_\{year\}\_\{name[j]\}.tif}: A thresholded
#'       raster layer based on the Mahalanobis distance for the future scenario.
#'   }
#'   It also prints plots of these raster layers.
#'
#' @examples
#' \dontrun{
#'   # Assuming you have the necessary data objects like
#'   # data_future_climatic_variables, polygon, future_climatic_variables,
#'   # names, reference_system, dir_future, model, and year defined in your environment.
#'   # Example usage might look like:
#'   # pa_mh_future(j = 1)
#' }
pa_mh_future <- function(j, th = .95) {
  data_f <- data_future_climatic_variables
  mh_f <- data.frame(matrix(1,
                            nrow = nrow(data_f),
                            ncol = length(names)))
  names(mh_f) <- names

  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(future_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)

    mh <- mahalanobis(data_f[,4:length(data_f)],
                      colMeans(data_polygon[,3:length(data_polygon)]),
                      cov(data_f[,4:length(data_f)]),
                      inverted = FALSE)

    mh_f[,i] <- mh
  }

  mh_f <- cbind(data_f[,c(1:3)], mh_f)
  mh_raster_f <- terra::rast(mh_f[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_f) <- colnames(mh_f[j+3])
  plot(mh_raster_f)
  terra::writeRaster(mh_raster_f, paste0(dir_future, "MH_", model, "_", year, "_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_", model, "_", year, "_", names[j]), mh_raster_f, envir = .GlobalEnv)

  puntos_todos_f <- terra::as.points(mh_raster_f)
  puntos_todos_f <- sf::st_as_sf(puntos_todos_f)
  colnames(puntos_todos_f) <- c("mh", "geometry")
  puntos_dentro_f <- sf::st_intersection(puntos_todos_f, pol)

  th_mh_f <- quantile(na.omit(puntos_dentro_f$mh), probs = th)
  puntos_todos_f$th <- dplyr::case_when(
    puntos_todos_f$mh > th_mh_f ~ 0,
    puntos_todos_f$mh <= th_mh_f  ~ 1
  )
  puntos_todos_f$th <- as.numeric(puntos_todos_f$th)

  res <- terra::res(mh_raster_f)
  bbox <- terra::ext(mh_raster_f)

  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- terra::rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_f <- terra::vect(puntos_todos_f)

  raster <- terra::rasterize(puntos_vect_f, raster_template, field = "th")
  terra::crs(raster) <- terra::crs(mh_raster_f)
  plot(raster)
  terra::writeRaster(raster,
                     paste0(dir_future, "TH_MH_", model, "_", year, "_", names[j], ".tif"),
                     overwrite = TRUE)

  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), mh_raster_p, envir = .GlobalEnv)
}
