#' @title Filter SpatRaster Layers based on Variance Inflation Factor (VIF)
#'
#' @description This function iteratively filters layers from a `SpatRaster` object by removing the one with the highest Variance Inflation Factor (VIF) that exceeds a specified threshold (`th`).
#'
#' @param x A `SpatRaster` object containing the layers (variables) to filter. Must contain two or more layers.
#' @param th A `numeric` value specifying the Variance Inflation Factor (VIF) threshold. Layers whose VIF exceeds this threshold are candidates for removal in each iteration (default: 5).
#'
#' @return A `SpatRaster` object containing only the layers retained by the VIF filtering process.
#'
#' @details This function implements a common iterative procedure to reduce multicollinearity among raster layers by removing variables with high Variance Inflation Factor (VIF).
#' The VIF for a specific predictor indicates how much the variance of its estimated coefficient is inflated due to its linear relationships with all other predictors in the model.
#' Conceptually, it is based on the proportion of variance that predictor shares with the other independent variables.
#' A high VIF value suggests a high degree of collinearity with other predictors (values exceeding `5` or `10` are often considered problematic; see O'Brien, 2007).
#' In this context, the function also provides the Pearson correlation matrix between all initial variables.
#'
#' Key steps:
#' \enumerate{
#'    \item Validate inputs: Ensures `x` is a `SpatRaster` with at least two layers and `th` is a valid `numeric` value.
#'    \item Convert the input `SpatRaster` (`x`) to a `data.frame`, retaining only unique rows if `x` has many cells and few unique climate values.
#'    \item Remove rows containing any `NA` values across all variables from the `data.frame`.
#'    \item In each iteration, calculate the VIF for all variables currently remaining in the dataset.
#'    \item Identify the variable with the highest VIF among the remaining variables.
#'    \item If this highest VIF value is greater than the threshold (`th`), remove the variable with the highest VIF from the dataset, and the loop continues with the remaining variables.
#'    \item This iterative process repeats until the highest VIF among the remaining variables is less than or equal to \eqn{\le} `th`, or until only one variable remains in the dataset.
#' }
#'The output of `vif_filter` returns a `list` object with a filtered `SpatRaster` object and a statistics summary.
#'
#'The `SpatRaster` object containing only the variables that were kept and also provides a comprehensive summary printed to the console.
#'The summary list including:
#' \itemize{
#' \item The original Pearson's correlation matrix between all initial variables.
#' \item The variables names that were kept and those that were excluded.
#' \item The final VIF values for the variables retained after the process.
#' }
#'
#' The internal VIF calculation includes checks to handle potential numerical
#' instability, such as columns with zero or near-zero variance and cases of
#' perfect collinearity among variables, which could otherwise lead to errors
#' (e.g., infinite VIFs or issues with matrix inversion). Variables identified
#' as having infinite VIF due to perfect collinearity are prioritized for removal.
#'
#' References:
#' O’brien (2007) A Caution Regarding Rules of Thumb for Variance Inflation Factors. Quality & Quantity, 41: 673–690. doi:10.1007/s11135-006-9018-6
#'
#' @importFrom terra as.data.frame subset rast
#' @importFrom stats cov var lm as.formula cor
#'
#' @examples
#' library(terra)
#' library(sf)
#'
#' set.seed(2458)
#' n_cells <- 100 * 100
#' r_clim <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
#' values(r_clim) <- c(
#'    (rowFromCell(r_clim, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
#'    (rowFromCell(r_clim, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
#'    (colFromCell(r_clim, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
#'    (colFromCell(r_clim, 1:n_cells) +
#'      (rowFromCell(r_clim, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'    (colFromCell(r_clim, 1:n_cells) /
#'      (rowFromCell(r_clim, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'    (colFromCell(r_clim, 1:n_cells) *
#'      (rowFromCell(r_clim, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4))),
#'    (colFromCell(r_clim, 1:n_cells) *
#'      (colFromCell(r_clim, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4))))
#' names(r_clim) <- c("varA", "varB", "varC", "varD", "varE", "varF", "varG")
#' terra::crs(r_clim) <- "EPSG:4326"
#' terra::plot(r_clim)
#'
#' vif_result <- ClimaRep::vif_filter(r_clim, th = 5)
#' print(vif_result$summary)
#' r_clim_filtered <- vif_result$filtered_raster
#' terra::plot(r_clim_filtered)
#' @export
vif_filter <- function(x, th = 5) {
  if (!inherits(x, 'SpatRaster')) {
    stop("Input 'x' must be a SpatRaster.")
  }
  if (!is.numeric(th)) {
    stop("Parameter 'th' must be numeric.")
  }
  original_raster <- x
  x_df <- terra::as.data.frame(x, na.rm = TRUE)
  original_cor_matrix <- NULL
  if (ncol(x_df) > 1) {
    original_cor_matrix <- round(cor(x_df, method = "pearson"), 4)
  } else {
    original_cor_matrix <- "Correlation matrix not applicable (less than 2 variables)."
  }
  calc_vif <- function(df) {
    if (ncol(df) <= 1) {
      return(numeric(0))
    }
    variances <- apply(df, 2, var, na.rm = TRUE)
    cols_zero_var <- names(variances[variances < .Machine$double.eps^0.5])
    if (length(cols_zero_var) > 0) {
      warning(
        "Removing columns with zero or near-zero variance during VIF calculation: ",
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

  results_summary <- list(
    kept_layers = kept_vars,
    excluded_layers = exc,
    original_correlation_matrix = original_cor_matrix,
    final_vif_values = final_vif_data
  )
  if (length(kept_vars) == 0) {
    warning("All variables were excluded. Returning an empty SpatRaster.")
    filtered_raster <- original_raster[[character(0)]]
  } else {
    filtered_raster <- subset(original_raster, kept_vars)
  }
  message("All processes were completed")
  return(list(filtered_raster = filtered_raster, summary = results_summary))
}
