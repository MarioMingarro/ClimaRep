#' @title Overlay Mahalanobis-based Climate Representativeness Change Classifications
#'
#' @description Combines multiple single-layer GeoTIFF rasters, typically outputs from `mh_rep_ch`
#' for different input polygons, into a single multi-layered raster. Each layer
#' in the output stack represents the count of how many input rasters had a
#' specific category value (e.g., Persistence, Loss, Gain) at each grid cell.
#' This allows visualizing spatial agreement or accumulation of change types across regions.
#'
#' @param folder_path Character. Path to the directory containing the input single-layer GeoTIFF
#' classification rasters (e.g., outputs from `mh_rep_ch`).
#' @param output_filename Character. The name for the output combined multi-layered GeoTIFF file.
#' The file will be saved within the `folder_path`. Defaults to "combined_category_counts_ap.tif".
#' @param category_values Numeric vector. A vector specifying the specific pixel
#'   values (`c(0, 1, 2, 3)`) to count in the input categories:
#'   0 = non-representative,
#'   1 = representative (mh_rep) or stable (mh_rep_ch),
#'   2 = loss (mh_rep_ch),
#'   3 = new (mh_rep_ch)). Defaults to `c(1, 2, 3)`.
#' @param add_to_environment Logical. If TRUE, the resulting multi-layered `SpatRaster`
#' object is assigned to a variable named `climarep_img` in the calling R environment.
#' Defaults to `FALSE`.
#'
#' @return Invisibly returns the resulting multi-layered `SpatRaster` object.
#' Writes the output stack to a GeoTIFF file in the specified `folder_path`.
#' If `add_to_environment` is TRUE, also assigns the stack to `climarep_img` in the parent environment.
#'
#' @details
#' This function is intended to be used after running `mh_rep` or `mh_rep_ch` for multiple
#' input polygons. Each run of `mh_rep_ch` produces a classification raster
#' (e.g., codes 0, 1, 2, 3) for a single polygon. `mh_overlay` takes these
#' individual polygon outputs and aggregates them by category.
#'
#' For each cell in the study area and for each specified `category_value`, the
#' function counts how many of the input classification rasters assigned that
#' specific `category_value` to that cell. The output is a stack where each
#' layer represents the counts for one category.
#'
#'
#' @importFrom terra rast ifel app writeRaster nlyr values
#'
#' @examples
#' \dontrun{
#' folder_path <- tempdir()
#' mh_overlay(folder_path,
#'   output_filename = "combined_category_counts.tif",
#'   category_values = c(1, 2, 3),
#'   add_to_environment = TRUE)
#' terra::plot(climarep_img)
#' }
#' @export
#'
mh_overlay <- function(folder_path, output_filename = "combined_category_counts.tif", category_values = c(1, 2, 3), add_to_environment = FALSE) {
  if (!is.character(folder_path) || length(folder_path) != 1) {
    stop("Parameter 'folder_path' must be a single character string.")
  }
  if (!is.character(output_filename) || length(output_filename) != 1) {
    stop("Parameter 'output_filename' must be a single character string.")
  }
  if (!is.numeric(category_values) || length(category_values) < 1) {
    stop("Parameter 'category_values' must be a numeric vector with at least one value (e.g. c(1, 2, 3)")
  }
  if (!is.logical(add_to_environment) || length(add_to_environment) != 1) {
    stop("Parameter 'add_to_environment' must be a single logical value (TRUE or FALSE).")
  }

  raster_files <- list.files(folder_path, pattern = "\\.tif$|\\.tiff$", full.names = TRUE, ignore.case = TRUE)
  if (length(raster_files) == 0) {
    message("No GeoTIFF files found in the specified folder: ", folder_path)
    return(NULL)
  }
  message("Found ", length(raster_files), " raster files. Processing...")
  binary_layers_by_category <- vector("list", length(category_values))
  names(binary_layers_by_category) <- paste0("Cat_", category_values)
  first_raster_props <- NULL
  for (i in seq_along(raster_files)) {
    file <- raster_files[i]
    message("Processing: ", basename(file), " (", i, "/", length(raster_files), ")")

    r <- try(terra::rast(file), silent = TRUE)
    if (inherits(r, "try-error")) {
      warning("Could not read raster file: ", basename(file), ". Skipping.")
      next
    }
    first_raster_props <- r
    if (terra::nlyr(r) > 1) {
      r <- r[[1]]
    }
    for (k in seq_along(category_values)) {
      cat_value <- category_values[k]
      binary_layer <- terra::ifel(r == cat_value, 1, 0)
      binary_layers_by_category[[k]] <- c(binary_layers_by_category[[k]], binary_layer)
    }
  }
  if (all(sapply(binary_layers_by_category, function(x) length(x) == 0))) {
    message("No valid raster data was processed from the found files.")
    return(NULL)
  }
  message("Stacking and summing binary layers by category...")
  output_layers <- list()
  for (k in seq_along(category_values)) {
    cat_value <- category_values[k]
    list_of_layers <- binary_layers_by_category[[k]]
    if (length(list_of_layers) == 0) {
      warning("No binary layers were generated for category ", cat_value, ". This category may not be present in any valid raster, or there were reading errors.")
      if (!is.null(first_raster_props)) {
        count_layer <- terra::rast(first_raster_props)
        terra::values(count_layer) <- 0
      } else {
        count_layer <- NULL
      }
    } else {
      stack_cat <- terra::rast(list_of_layers)
      count_layer <- terra::app(stack_cat, fun = sum, na.rm = TRUE)
    }
    if (!is.null(count_layer)) {
      names(count_layer) <- paste0("Count_Category_", cat_value)
      output_layers[[length(output_layers) + 1]] <- count_layer
    }
  }
  if (length(output_layers) == 0) {
    message("No output layers could be generated for the specified categories.")
    return(NULL)
  }
  output_stack <- terra::rast(output_layers)
  dir_output <- file.path(folder_path, output_filename)
  terra::writeRaster(output_stack, dir_output, overwrite = TRUE, datatype = "INT4U")
  if (add_to_environment) {
    message("Assigning final raster stack to variable 'climarep_img' in the calling environment.")
    assign("climarep_img", output_stack, envir = parent.frame())
  }
  message("All processes were completed")
  cat(paste("Output files in:", dir_output))
  return(invisible(output_stack))
}
