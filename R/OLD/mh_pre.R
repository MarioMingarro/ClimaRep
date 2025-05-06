# Function to process the data of the present series
#' Function to process the data of the present series
#'
#' This function processes the data for the present time series to calculate
#' Mahalanobis distance and a thresholded raster.
#'
#' @param j An index or identifier for the variable being processed.
#' @param th Threshold probability for defining the presence/absence based on
#'   the Mahalanobis distance. Default is 0.95.
#'
#' @return This function does not explicitly return an object but creates and
#'   assigns several raster objects to the global environment and writes them
#'   to disk. These include:
#'   \itemize{
#'     \item \code{MH\_PRESENT\_\{name[j]\}.tif}: A raster layer of Mahalanobis
#'       distances.
#'     \item \code{TH\_MH\_PRESENT\_\{name[j]\}.tif}: A thresholded raster layer
#'       based on the Mahalanobis distance.
#'   }
#'   It also prints plots of these raster layers.
#'
#' @examples
#' \dontrun{
#'   # Assuming you have the necessary data objects like
#'   # data_present_climatic_variables, polygon, present_climatic_variables,
#'   # names, reference_system, and dir_present defined in your environment.
#'   # Example usage might look like:
#'   # pa_mh_present(j = 1)
#' }
pa_mh_present <- function(j, th = .95) {
  data_p <- data_present_climatic_variables
  mh_p <- data.frame(matrix(1,
                            nrow = nrow(data_p),
                            ncol = length(names)))

  names(mh_p) <- names

  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)

    mh <- mahalanobis(data_p[,4:length(data_p)],
                      colMeans(data_polygon[,3:length(data_polygon)]),
                      cov(data_p[,4:length(data_p)]),
                      inverted = FALSE)

    mh_p[,i] <- mh
  }

  mh_p <- cbind(data_p[,c(1:3)], mh_p)
  mh_raster_p <- terra::rast(mh_p[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[j+3])
  plot(mh_raster_p)
  terra::writeRaster(mh_raster_p, paste0(dir_present, "MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)


  puntos_todos_p <- terra::as.points(mh_raster_p)
  puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
  colnames(puntos_todos_p) <- c("mh", "geometry")
  puntos_dentro_p <- sf::st_intersection(puntos_todos_p, pol)

  th_mh_p <- quantile(na.omit(puntos_dentro_p$mh), probs = th)
  puntos_todos_p$th <- dplyr::case_when(
    puntos_todos_p$mh > th_mh_p ~ 0,
    puntos_todos_p$mh <= th_mh_p  ~ 1
  )
  puntos_todos_p$th <- as.numeric(puntos_todos_p$th)

  res <- terra::res(mh_raster_p)
  bbox <- terra::ext(mh_raster_p)

  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- terra::rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_p <- terra::vect(puntos_todos_p)

  raster <- terra::rasterize(puntos_vect_p, raster_template, field = "th")
  terra::crs(raster) <- terra::crs(mh_raster_p)
  plot(raster)
  terra::writeRaster(raster, paste0(dir_present, "TH_MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)

}
