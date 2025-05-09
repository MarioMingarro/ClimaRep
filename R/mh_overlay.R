#' @title Count Category Overlap Across Raster Layers for Protected Areas Analysis
#'
#' @description
#' This function processes multiple GeoTIFF files from a specified folder,
#' treating them as layers containing categorical information relevant to
#' Protected Area (AP) analysis (e.g., representative, non-representative,
#' stable, future areas). It counts, for each pixel, how many input rasters
#' classify that pixel into each of the specified categories. The output is a
#' `SpatRaster` stack where each layer shows the total count of occurrences
#' for a specific category across all input layers at every location. This is
#' useful for assessing the spatial consistency and frequency of different
#' classifications in representativeness studies.
#'
#' @param folder_path A character string specifying the path to the folder
#'   containing the input GeoTIFF (.tif or .tiff) files. These files should
#'   represent different spatial datasets or scenarios classifying areas into
#'   the categories specified in `category_values`.
#' @param output_filename A character string specifying the name for the output
#'   GeoTIFF file. This file will contain a stack of layers, one for each
#'   category value specified in `category_values`, showing the count of
#'   occurrences. Defaults to "combined_category_counts_ap.tif".
#' @param category_values A numeric vector specifying the specific category
#'   values (pixel values) to count in the input rasters (e.g., `c(1, 2, 3, 4)`
#'   where 1=representative, 2=non-representative, etc.). Defaults to `c(1, 2, 3)`.
#'
#' @details
#' This function is designed to help analyze spatial datasets where multiple
#' raster layers provide potentially different classifications for the same
#' geographic area, common in ecological modeling or conservation planning
#' studies (e.g., predicting species distribution, habitat suitability, or
#' conservation status under different scenarios).
#'
#' The process involves:
#' 1. Listing all GeoTIFF files in the `folder_path`.
#' 2. Reading each file. If a raster has multiple layers, only the first is used.
#' 3. For each specified category value, a binary layer (1 or 0) is created
#'    from each input raster, indicating where that category is present.
#' 4. The binary layers are grouped by category, and for each category, the
#'    layers are summed across all input rasters. The resulting sum layer
#'    shows, for each pixel, the *number of input rasters* where that pixel
#'    belonged to that specific category.
#' 5. The sum layers for all categories are combined into a single output
#'    `SpatRaster` stack.
#' 6. The final stack is saved as a GeoTIFF file. The data type is set to
#'    "INT2U", suitable for counts up to 65535. Adjust if higher counts are
#'    expected.
#'
#' The function assumes the input rasters have compatible spatial properties
#' (extent, resolution, CRS), relying on `terra`'s ability to handle minor
#' differences during stacking and processing. It includes messaging for
#' progress and warnings for issues like unreadable files or categories not
#' found in the data.
#'
#' The output stack provides valuable input for subsequent analysis of spatial
#' agreement or frequency of classifications, aiding in identifying core areas
#' for each category or areas of high uncertainty.
#'
#' @return
#' An invisible `SpatRaster` object containing the stacked count layers. Each
#' layer corresponds to a category value from `category_values` and shows
#' the number of input rasters where that category was present at each location.
#' Returns `NULL` if no valid raster files are found or processed.
#'
#' @importFrom terra rast list.files ifel app writeRaster nlyr values file.path
#' @importFrom methods inherits
#' @importFrom base basename file.path seq_along sapply is.null length warning message
#'
#' @examples
#' \dontrun{
#' # requires terra package
#' library(terra)
#'
#' # Create a temporary directory and some dummy raster files
#' # simulating different scenarios or datasets for AP categories
#' temp_dir <- file.path(tempdir(), "ap_category_overlap_example")
#' dir.create(temp_dir, showWarnings = FALSE)
#'
#' # Simulate 3 raster layers. Values represent categories:
#' # 1 = Representative, 2 = Non-representative, 3 = Stable, 4 = Future
#' r_scenario1 <- rast(ncols=15, nrows=15, xmin=0, xmax=15, ymin=0, ymax=15)
#' values(r_scenario1) <- sample(c(1, 2, 3, 4, 99, NA), size=ncell(r_scenario1), replace=TRUE, prob=c(0.2, 0.3, 0.1, 0.1, 0.2, 0.1))
#' writeRaster(r_scenario1, file.path(temp_dir, "ap_cats_scenario_a.tif"), overwrite=TRUE)
#'
#' r_scenario2 <- rast(ncols=15, nrows=15, xmin=0, xmax=15, ymin=0, ymax=15)
#' values(r_scenario2) <- sample(c(1, 2, 3, 4, 50, NA), size=ncell(r_scenario2), replace=TRUE, prob=c(0.3, 0.2, 0.1, 0.1, 0.2, 0.1))
#' writeRaster(r_scenario2, file.path(temp_dir, "ap_cats_scenario_b.tif"), overwrite=TRUE)
#'
#' r_scenario3 <- rast(ncols=15, nrows=15, xmin=0, xmax=15, ymin=0, ymax=15)
#' values(r_scenario3) <- sample(c(1, 2, 3, 4, 99, NA), size=ncell(r_scenario3), replace=TRUE, prob=c(0.1, 0.1, 0.3, 0.3, 0.1, 0.1))
#' writeRaster(r_scenario3, file.path(temp_dir, "ap_cats_scenario_c.tif"), overwrite=TRUE)
#'
#' # Define the categories relevant to your AP analysis
#' ap_categories <- c(1, 2, 3, 4) # Representative, Non-repr, Stable, Future
#'
#' # Process the rasters, counting occurrences of these categories
#' output_stack <- conteo_superposicion_categorias_ap(
#'   folder_path = temp_dir,
#'   output_filename = "ap_category_overlap_counts.tif",
#'   category_values = ap_categories
#' )
#'
#' # Print the resulting stack (invisible by default, but can be captured)
#' if (!is.null(output_stack)) {
#'   print(output_stack)
#'   # You can inspect the counts for each category layer
#'   # For example, plot the count for the 'Representative' category (assuming it was category 1)
#'   # plot(output_stack[["Count_Category1"]], main="Count of 'Representative' Category Overlap")
#' }
#'
#' # Clean up the temporary files and directory
#' unlink(temp_dir, recursive = TRUE)
#'}
#' @export
mh_overlay <- function(folder_path, output_filename = "combined_category_counts_ap.tif", category_values = c(1, 2, 3)) {

  # 1. Find all GeoTIFF files in the folder
  raster_files <- list.files(folder_path, pattern = "\\.tif$|\\.tiff$", full.names = TRUE, ignore.case = TRUE)

  if (length(raster_files) == 0) {
    message("No se encontraron archivos GeoTIFF en la carpeta especificada: ", folder_path)
    return(NULL)
  }

  message("Found ", length(raster_files), " raster files. Processing...")

  # Initialize lists to store binary rasters for each category
  binary_layers_by_category <- vector("list", length(category_values))
  # Name the list elements by category value for easier tracking
  names(binary_layers_by_category) <- paste0("Cat_", category_values)

  # Store properties of the first valid raster as a reference
  first_raster_props <- NULL

  # 2. Process each raster file
  for (i in seq_along(raster_files)) {
    file <- raster_files[i]
    message("Processing: ", basename(file), " (", i, "/", length(raster_files), ")")

    # Try to read the raster, handling potential errors
    r <- try(terra::rast(file), silent = TRUE) # Use silent=TRUE to avoid printing default try-error messages
    if (inherits(r, "try-error")) {
      warning("Could not read raster file: ", basename(file), ". Skipping.")
      next
    }

    # Ensure it's a single-layer raster (take the first layer if multi-layer)
    if (terra::nlyr(r) > 1) {
      r <- r[[1]]
      warning("Raster '", basename(file), "' has multiple layers, using only the first.")
    }

    # For the first valid raster, store its properties
    if (is.null(first_raster_props)) {
      first_raster_props <- r
    } else {
      # Optional: Add strict alignment checks here if needed, e.g.:
      # if (!terra::compareCRS(r, first_raster_props)) warning("CRS mismatch for ", basename(file))
      # if (terra::ext(r) != terra::ext(first_raster_props)) warning("Extent mismatch for ", basename(file))
      # if (any(terra::res(r) != terra::res(first_raster_props))) warning("Resolution mismatch for ", basename(file))
    }

    # Create binary layers (0/1) for each specified category value
    for (k in seq_along(category_values)) {
      cat_value <- category_values[k]
      # ifel(condition, value_if_TRUE, value_if_FALSE)
      # Assign 1 if cell value equals category_value, otherwise 0.
      # Original NAs that don't match the category become 0.
      binary_layer <- terra::ifel(r == cat_value, 1, 0)

      # Add the binary layer to the list corresponding to the category
      binary_layers_by_category[[k]] <- c(binary_layers_by_category[[k]], binary_layer)
    }
  }

  # Check if any valid rasters were processed and binary layers created
  # We check if *any* category list received *any* layer
  if (all(sapply(binary_layers_by_category, function(x) length(x) == 0))) {
    message("No valid raster data was processed from the found files.")
    return(NULL)
  }

  # 3. Stack binary layers per category and calculate the sum
  message("Stacking and summing binary layers by category...")

  output_layers <- list()

  for (k in seq_along(category_values)) {
    cat_value <- category_values[k]
    list_of_layers <- binary_layers_by_category[[k]]

    if (length(list_of_layers) == 0) {
      warning("No binary layers were generated for category ", cat_value, ". This category may not be present in any valid raster, or there were reading errors.")
      # Create an empty layer with the correct properties if no data for this category but some rasters were valid
      if (!is.null(first_raster_props)) {
        count_layer <- terra::rast(first_raster_props) # Create a raster with same ext/res/crs
        terra::values(count_layer) <- 0 # Fill with zeros
      } else {
        # If no valid rasters were processed at all, skip creating this layer
        count_layer <- NULL
      }
    } else {
      # Stack the binary layers for this category
      # If only one layer for a category, rast() handles it correctly
      stack_cat <- terra::rast(list_of_layers)

      # Calculate the sum (count) for this category across all rasters
      # na.rm = TRUE is important if rasters have different extents or NAs
      count_layer <- terra::app(stack_cat, fun = sum, na.rm = TRUE)
    }

    # Assign name to the count layer and add to output list if created
    if (!is.null(count_layer)) {
      names(count_layer) <- paste0("Count_Category", cat_value)
      output_layers[[length(output_layers) + 1]] <- count_layer # Append to output_layers
    } else {
      # Warning already given above if count_layer is NULL
    }
  }

  # Verify if any output layers were successfully created
  if (length(output_layers) == 0) {
    message("No output layers could be generated for the specified categories.")
    return(NULL)
  }

  # Combine the resulting count layers into a single stack
  output_stack <- terra::rast(output_layers)

  # 4. Write the output raster
  output_path <- file.path(folder_path, output_filename)
  message("Writing output raster stack to: ", output_path)
  # Use datatype suitable for counts. INT2U (unsigned 16-bit) max value 65535.
  # Change to "INT4U" (unsigned 32-bit, max 4,294,967,295) if more layers are expected.
  terra::writeRaster(output_stack, output_path, overwrite = TRUE, datatype = "INT2U")

  message("Processing complete.")

  # Return the output stack (invisibly)
  return(invisible(output_stack))
}
