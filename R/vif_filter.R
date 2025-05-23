#' @title Filter SpatRaster Layers based on Variance Inflation Factor (VIF)
#'
#' @description This function iteratively filters layers from a `SpatRaster` object by removing the one with the highest Variance Inflation Factor (VIF) that exceeds a specified threshold (`th`).
#'
#' @param x A `SpatRaster` object containing the layers (variables) to filter.
#'   Must contain two or more layers.
#' @param th A numeric value specifying the Variance Inflation Factor (VIF)
#'   threshold. Layers whose VIF exceeds this threshold are candidates for
#'   removal in each iteration (default: 10).
#'
#' @return A [SpatRaster] object containing onlythe layers from the input
#' `x` that were retained by the VIF filtering process. The layers are returned
#' in their original order. If no layers meet the VIF threshold criterion
#' (all are excluded) or if the input becomes empty after removing NA values, an empty `SpatRaster` object is returned.
#'
#' @details This function implements a common iterative procedure to reduce multicollinearity among raster layers by removing variables with high Variance Inflation Factor (VIF).
#' The VIF for a specific predictor indicates how much the variance of its estimated coefficient is inflated due to its linear relationships with all other predictors in the model.
#' Conceptually, it is based on the proportion of variance that predictor shares with the other independent variables.
#' A high VIF value suggests a high degree of collinearity with other predictors (values exceeding `5` or `10` are often considered problematic; see O'Brien, 2007).
#' In this context, the function provides the Pearson correlation matrix between all initial variables, and a threshold parameter (`th`) is available to specify the critical VIF value for addressing multicollinearity.
#'
#' Here are the key steps:
#' \enumerate{
#'   \item Convert the input `SpatRaster` (`x`) to a `data.frame`.
#'   \item Remove rows containing any `NA` values across all variables from the `data.frame`.
#'   \item In each iteration, calculate the Variance Inflation Factor (VIF) for all variables currently remaining in the dataset.
#'   \item Identify the variable with the highest VIF among the remaining variables.
#'   \item If this highest VIF value is greater than the specified threshold (`th`), remove the variable with the highest VIF from the dataset, and the loop continues with the remaining variables.
#'   \item This iterative process (steps 3-5) repeats until the highest VIF among the remaining variables is less than or equal to `th`, or until only one variable remains in the dataset.
#' }
#'
#' Finally, the function returns a new `SpatRaster` object containing only the variables that were kept.
#' It also prints a summary including:
#' \itemize{
#' \item The original Pearson's correlation matrix between all initial variables.
#' \item The lists of variables that were kept and those that were excluded.
#' \item The final VIF values for the variables that were retained after the filtering process.
#' }
#'
#' The internal VIF calculation includes checks to handle potential numerical
#' instability, such as columns with zero or near-zero variance and cases of
#' perfect collinearity among variables, which could otherwise lead to errors
#' (e.g., infinite VIFs or issues with matrix inversion). Variables identified
#' as having infinite VIF due to perfect collinearity are prioritized for removal.
#'
#' References:
#' O’brien, R.M. A Caution Regarding Rules of Thumb for Variance Inflation Factors. Qual Quant 41, 673–690 (2007). https://doi.org/10.1007/s11135-006-9018-6
#'
#' @importFrom terra as.data.frame subset
#' @importFrom stats cov var lm as.formula cor
#' @importFrom utils packageVersion
#'
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
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
#' names(r_clim) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
#' terra::crs(r_clim) <- "EPSG:4326"
#' terra::plot(r_clim)
#' r_clim_filtered <- vif_filter(r_clim, th = 5)
#' terra::plot(r_clim_filtered)
#'
#'}
#' @export
#'
vif_filter <- function(x, th = 5) {
  if (!inherits(x, 'SpatRaster')) {
    stop("Input 'x' must be a SpatRaster object to return a filtered raster.")
  }
  if (!is.numeric(th)) {
    stop("Parameter 'th' must be numeric.")
  }
  original_raster <- x
  x_df <- terra::as.data.frame(x, na.rm = TRUE)
  if (nrow(x_df) == 0 || ncol(x_df) == 0) {
    warning(
      "Data frame is empty after removing NAs. Cannot perform VIF calculation. Returning an empty SpatRaster.")
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
      warning(
        "Removing columns with zero or near-zero variance during VIF calculation:",
        paste(cols_zero_var, collapse = ", ")
      )
      df <- df[, !(colnames(df) %in% cols_zero_var), drop = FALSE]
      if (ncol(df) <= 1) {
        return(numeric(0))
      }
    }
    vif_values <- sapply(1:ncol(df), function(i) {
      model <- try(stats::lm(as.formula(paste(names(df)[i], "~ .")), data = df), silent = TRUE)
      if (inherits(model, "try-error") ||
          is.null(summary(model)$r.squared) ||
          is.na(summary(model)$r.squared) ||
          summary(model)$r.squared >= 1) {
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
      warning("Variable ",
              ex,
              " with max VIF already in excluded list. Breaking loop.")
      break
    }
    exc <- c(exc, ex)
    kept_vars <- kept_vars[!(kept_vars %in% ex)]
  }
  final_vif_data <- NULL
  kept_df_subset <- x_df[, kept_vars, drop = FALSE]
  if (length(kept_vars) > 0) {
    if (length(kept_vars) > 1) {
      final_vif_values <- round(calc_vif(kept_df_subset), 4)
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
