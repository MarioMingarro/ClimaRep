library(testthat)
library(terra)


testthat::test_that("mh_overlay counts categories correctly and writes output", {
  temp_dir <- file.path(tempdir(), "test_mh_overlay")
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }
  dir.create(temp_dir, showWarnings = FALSE)
  base_rast <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10, crs="EPSG:32632")
  n_cells_dummy <- ncell(base_rast)
  r1 <- base_rast
  values(r1) <- sample(c(1, 2, 3, NA), size=n_cells_dummy, replace=TRUE, prob = c(0.2, 0.5, 0.3, 0.2))
  file1 <- file.path(temp_dir, "raster_A.tif")
  writeRaster(r1, file1, overwrite=TRUE)
  r2 <- base_rast
  values(r2) <- sample(c(1, 2, 3, NA), size=n_cells_dummy, replace=TRUE, prob = c(0.2, 0.5, 0.3, 0.2))
  file2 <- file.path(temp_dir, "raster_B.tif")
  writeRaster(r2, file2, overwrite=TRUE)
  r3 <- base_rast
  values(r3) <- sample(c(1, 2, NA), size=n_cells_dummy, replace=TRUE, prob = c(0.4, 0.3, 0.3))
  file3 <- file.path(temp_dir, "raster_C.tif")
  writeRaster(r3, file3, overwrite=TRUE)

  categories_to_count <- c(1, 2, 3)
  output_file_name <- "test_category_counts.tif"
  output_path_full <- file.path(temp_dir, output_file_name)

  message("\nRunning mh_overlay test...")
  output_raster_stack <- mh_overlay(
    folder_path = temp_dir,
    output_filename = output_file_name,
    category_values = categories_to_count
  )
  message("mh_overlay call finished.")
  expect_true(file.exists(output_path_full),
              info = paste("Output file", output_file_name, "was not created in", temp_dir))
  expect_true(inherits(output_raster_stack, "SpatRaster"),
              info = "Function did not return a SpatRaster object")
  expect_equal(nlyr(output_raster_stack), length(categories_to_count),
               info = paste("Output raster stack should have", length(categories_to_count), "layers"))
  expected_layer_names <- paste0("Count_Category", categories_to_count)
  expect_equal(names(output_raster_stack), expected_layer_names,
               info = "Output layer names are not as expected")
})
testthat::teardown({
  message("\nCleaning up temporary directory: ", temp_dir)
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }
})
