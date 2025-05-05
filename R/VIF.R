#' @title Filter SpatRaster Layers based on Variance Inflation Factor (VIF)
#'
#' @description
#' This function iteratively filters layers from a `SpatRaster` object by removing the one with the highest Variance Inflation Factor (VIF) that exceeds a specified threshold (`th`). The process continues until all remaining layers have a VIF below the threshold or only one layer remains. VIF calculation is performed on the raster data converted to a `data.frame` after removing rows containing `NA` values in any column.
#'
#' @param x A `SpatRaster` object containing the layers (variables) to filter. Must be a valid `SpatRaster` object with multiple layers.
#' @param th A numeric value specifying the Variance Inflation Factor (VIF) threshold. Layers whose VIF exceeds this threshold are candidates for removal in each step. Defaults to 10.
#'
#' @details
#' The Variance Inflation Factor (VIF) quantifies the severity of multicollinearity in a linear regression analysis. A high VIF for a variable indicates that this variable is strongly correlated with other predictor variables. In the context of environmental variable selection, a high VIF suggests redundancy of information among variables.
#'
#' The function follows a common iterative procedure to reduce multicollinearity:
#' 1.  The input `SpatRaster` is converted to a `data.frame`.
#' 2.  Rows containing any `NA` values are removed. All VIF calculations are based on this cleaned dataset.
#' 3.  The VIF is calculated for each remaining variable.
#' 4.  If the highest VIF is greater than the threshold `th`, the corresponding variable is removed from the dataset.
#' 5.  Steps 3 and 4 are repeated until the highest VIF among the remaining variables is less than or equal to `th`, or until only one variable remains.
#'
#' During the process, the function prints messages to the console indicating which variables are being evaluated and which, if any, are removed.
#'
#' The internal function for calculating VIF includes checks to handle problematic situations such as columns with zero or near-zero variance and cases of perfect collinearity, which could otherwise cause errors in VIF calculation.
#'
#' @return
#' A `SpatRaster` object containing only the layers from the input `SpatRaster` that were retained by the VIF filter. The order of the remaining layers is kept according to the original `SpatRaster`. If all layers are excluded during the process, an empty `SpatRaster` is returned.
#'
#' @importFrom terra SpatRaster as.data.frame subset
#' @importFrom stats cov
#'
#' @examples
#' # requires terra package
#' library(terra)
#'
#' # Create a simple SpatRaster with 3 layers, one highly correlated
#' r <- rast(ncols=5, nrows=5) # Small raster
#' values(r) <- cbind(1:ncell(r), (1:ncell(r))*1.5 + rnorm(ncell(r), 0, 0.1), runif(ncell(r)))
#' names(r) <- c("v1", "v_corr", "v_rand")
#'
#' # Apply the VIF filtering function (replace 'your_vif_filter_function'
#' # with the actual name of the function)
#' # filtered_r <- your_vif_filter_function(x = r, th = 5)
#'
#' # Print the result (replace 'filtered_r' as needed)
#' # print(filtered_r)
#'
#' @export # Export the function if it's part of a package
#'
# Followed by the actual function code:
# your_vif_filter_function <- function(x, th = 10) {
#   # ... function implementation ...
# }

vif_filter <- function(x, th = 10) {
  if (!inherits(x, 'SpatRaster')) {
    stop("Input 'x' must be a SpatRaster object to return a filtered raster.")
  }
  original_raster <- x
  x_df <- terra::as.data.frame(x, na.rm = TRUE)
  if (nrow(x_df) == 0 || ncol(x_df) == 0) {
    warning("Data frame is empty after removing NAs. Cannot perform VIF calculation. Returning an empty SpatRaster.")
    return(original_raster[[character(0)]])
  }

  original_cor_matrix <- NULL
  if (ncol(x_df) > 1) {
    original_cor_matrix <- round(cor(x_df, method = "pearson"), 4)
  } else {
    original_cor_matrix <- "Correlation matrix not applicable (less than 2 original variables after removing NAs)."
  }

  calc_vif <- function(df) {
    if (ncol(df) <= 1) {
      return(numeric(0))
    }
    variances <- apply(df, 2, var, na.rm = TRUE)
    cols_zero_var <- names(variances[variances < .Machine$double.eps^0.5])
    if (length(cols_zero_var) > 0) {
      warning("Removing columns with zero or near-zero variance during VIF calculation:", paste(cols_zero_var, collapse = ", "))
      df <- df[, !(colnames(df) %in% cols_zero_var), drop = FALSE]
      if (ncol(df) <= 1) {
        return(numeric(0))
      }
    }
    vif_values <- sapply(1:ncol(df), function(i) {
      model <- try(lm(as.formula(paste(names(df)[i], "~ .")), data = df), silent = TRUE)
      if (inherits(model, "try-error") || is.null(summary(model)$r.squared) || is.na(summary(model)$r.squared) || summary(model)$r.squared >= 1) {
        return(Inf)
      }
      vif <- 1 / (1 - summary(model)$r.squared)
      if (is.infinite(vif)) {
        return(Inf)
      }
      return(vif)
    })
    names(vif_values) <- colnames(df)
    return(vif_values)
  }

  exc <- character(0)
  kept_vars <- colnames(x_df)

  while (length(kept_vars) > 1) {
    df_subset <- x_df[, kept_vars, drop = FALSE]
    v <- calc_vif(df_subset)
    if (length(v) == 0 || all(v < th)) {
      break
    }

    max_v_val <- max(v, na.rm = TRUE)
    if (is.infinite(max_v_val)) {
      ex <- names(v)[which(is.infinite(v))[1]]
    } else {
      ex <- names(v)[which.max(v)]
    }

    if (ex %in% exc) {
      warning("Variable ", ex, " with max VIF already in excluded list. Breaking loop.")
      break
    }

    exc <- c(exc, ex)
    kept_vars <- kept_vars[!(kept_vars %in% ex)]
  }

  final_vif_data <- NULL
  kept_df_subset <- x_df[, kept_vars, drop = FALSE]

  if (length(kept_vars) > 0) {
    if (length(kept_vars) > 1) {
      final_vif_values <- round(calc_vif(kept_df_subset),4)
      if (length(final_vif_values) > 0) {
        final_vif_data <- data.frame(VIF = final_vif_values)
      } else {
        final_vif_data <- "Could not calculate VIFs for kept variables (e.g., perfect collinearity remaining or insufficient variables after zero-variance removal)."
      }
    } else {
      final_vif_data <- "Only one variable kept. VIF calculation not applicable."
    }
  } else {
    final_vif_data <- "No variables kept."
  }

  cat("--- VIF Filtering Summary ---\n")
  cat("VIF filtering completed.\n")
  cat("Kept layers:", paste(kept_vars, collapse = ", "), "\n")
  cat("Excluded layers:", paste(exc, collapse = ", "), "\n")

  cat("\nPearson correlation matrix of original data:\n")
  if (is.matrix(original_cor_matrix)) {
    print(original_cor_matrix)
  } else {
    cat(original_cor_matrix, "\n")
  }

  cat("\nFinal VIF values for kept variables:\n")
  if (is.data.frame(final_vif_data)) {
    print(final_vif_data)
  } else {
    cat(final_vif_data, "\n")
  }

  cat("----------------------------\n")

  if (length(kept_vars) == 0) {
    warning("All variables were excluded. Returning an empty SpatRaster.")
    return(original_raster[[character(0)]])
  }

  result_raster <- subset(original_raster, kept_vars)

  return(result_raster)
}
