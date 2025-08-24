#' @title Filter SpatRaster layers based on Variance Inflation Factor (VIF)
#'
#' @description This function iteratively filters layers from a `SpatRaster` object by removing the one with the highest Variance Inflation Factor (VIF) that exceeds a specified threshold (`th`).
#'
#' @param x A `SpatRaster` object containing the layers (variables) to filter. Must contain two or more layers.
#' @param th A `numeric` value specifying the Variance Inflation Factor (VIF) threshold. Layers whose VIF exceeds this threshold are candidates for removal in each iteration (default: 5).
#'
#' @return A `list` object containing the filtered `SpatRaster` and a summary of the filtering process.
#'
#' @details This function implements a common iterative procedure to reduce multicollinearity among raster layers by removing variables with a high Variance Inflation Factor (VIF).
#' The VIF for a specific predictor indicates how much the variance of its estimated coefficient is inflated due to its linear relationships with all other predictors in the model. A high VIF value suggests a high degree of collinearity with other predictors (values exceeding `5` or `10` are often considered problematic; see O'Brien, 2007; Legendre & Legendre, 2012).
#'
#' **The filtering process is fully automated and robust:**
#'
#' \enumerate{
#' \item Validates the input and converts the `SpatRaster` to a `data.frame` for calculations.
#' \item In each step, the function attempts to calculate VIF efficiently using matrix inversion. If perfect collinearity is detected (resulting in a singular matrix that cannot be inverted), the function automatically switches to a more robust method based on linear regressions to handle the situation without an error.
#' \item The function identifies the variable with the highest VIF among the remaining variables. If its VIF is greater than the threshold (`th`), that variable is removed. The process repeats until all remaining variables are below the threshold or until only one variable remains.
#' }
#'
#' The output is a `list` containing two main components:
#' \itemize{
#' \item `SpatRaster` object with the variables that were retained after the filtering process.
#' \item A list with a detailed summary of the process, including the names of the kept and excluded variables, the original Pearson's correlation matrix, and the final VIF values for the retained variables.
#' }
#' The internal VIF calculation includes checks to handle potential numerical instability, such as columns with zero or near-zero variance and cases of perfect collinearity among variables,
#' which could otherwise lead to errors (e.g., infinite VIFs). Variables identified as having infinite VIF due to perfect collinearity are prioritized for removal.
#'
#' References:
#' O’Brien, R. M. (2007). A caution regarding rules of thumb for variance inflation factors. Quality & Quantity, 41(5), 673–690. https://doi.org/10.1007/s11135-006-9018-6
#' Legendre, P., & Legendre, L. (2012). Interpretation of ecological structures. In P. Legendre & L. Legendre (Eds.), *Developments in Environmental Modelling* (Vol. 24, pp. 521-624). Elsevier. https://doi.org/10.1016/B978-0-444-53868-0.50010-1
#'
#' @importFrom terra as.data.frame subset rast
#' @importFrom stats cov var lm as.formula cor
#'
#' @examples
#' library(terra)
#'
#' set.seed(2458)
#' n_cells <- 100 * 100
#' r_clim <- terra::rast(ncols = 100, nrows = 100, nlyrs = 7)
#' values(r_clim) <- c(
#'   (rowFromCell(r_clim, 1:n_cells) * 0.2 + rnorm(n_cells, 0, 3)),
#'   (rowFromCell(r_clim, 1:n_cells) * 0.9 + rnorm(n_cells, 0, 0.2)),
#'   (colFromCell(r_clim, 1:n_cells) * 0.15 + rnorm(n_cells, 0, 2.5)),
#'   (colFromCell(r_clim, 1:n_cells) +
#'     (rowFromCell(r_clim, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'   (colFromCell(r_clim, 1:n_cells) /
#'     (rowFromCell(r_clim, 1:n_cells)) * 0.1 + rnorm(n_cells, 0, 4)),
#'   (colFromCell(r_clim, 1:n_cells) *
#'     (rowFromCell(r_clim, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4))),
#'   (colFromCell(r_clim, 1:n_cells) *
#'     (colFromCell(r_clim, 1:n_cells) + 0.1 + rnorm(n_cells, 0, 4))))
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
  calc_vif_robust <- function(df) {
    if (ncol(df) <= 1) return(numeric(0))
    tryCatch({
      cor_matrix <- cor(df, use = "pairwise.complete.obs")
      vif_values <- diag(solve(cor_matrix))
      names(vif_values) <- colnames(df)
      return(vif_values)
    },
    error = function(e) {
      if (grepl("singular", e$message)) {
        message("Perfect collinearity detected. Switching to the regression method for VIF calculation.")
        vif_values <- sapply(names(df), function(name) {
          model <- try(stats::lm(as.formula(paste(name, "~ .")), data = df), silent = TRUE)
          if (inherits(model, "try-error") || summary(model)$r.squared >= 1) {
            return(Inf)
          }
          1 / (1 - summary(model)$r.squared)
        })
        names(vif_values) <- names(df)
        return(vif_values)
      } else {
        stop(e)
      }
    })
  }
  message("Starting iterative VIF filtering process")
  kept_vars <- colnames(x_df)
  excluded_vars <- character(0)
  while (length(kept_vars) > 1) {
    df_subset <- x_df[, kept_vars, drop = FALSE]
    v <- calc_vif_robust(df_subset)
    if (any(is.infinite(v))) {
      max_vif_name <- names(v)[which(is.infinite(v))[1]]
      kept_vars <- setdiff(kept_vars, max_vif_name)
      excluded_vars <- c(excluded_vars, max_vif_name)
      message(paste(" - Removing variable with infinite VIF:", max_vif_name))
      next
    }
    if (all(v < th)) {
      break
    }
    max_vif_name <- names(v)[which.max(v)]
    if (v[max_vif_name] > th) {
      kept_vars <- setdiff(kept_vars, max_vif_name)
      excluded_vars <- c(excluded_vars, max_vif_name)
      message(paste(" - Removing variable with highest VIF (", round(v[max_vif_name], 2), "):", max_vif_name))
    } else {
      break
    }
  }
  message("VIF filtering process completed")
  final_vif_values <- NULL
  if (length(kept_vars) > 1) {
    final_vif_values <- round(calc_vif_robust(x_df[, kept_vars, drop = FALSE]), 4)
  } else if (length(kept_vars) == 1) {
    final_vif_values <- "Only one variable was kept. VIF calculation is not applicable."
  } else {
    final_vif_values <- "No variables were kept."
  }
  original_cor_matrix <- round(cor(x_df, method = "pearson"), 4)
  results_summary <- list(
    kept_layers = kept_vars,
    excluded_layers = excluded_vars,
    original_correlation_matrix = original_cor_matrix,
    final_vif_values = final_vif_values
  )
  filtered_raster <- if (length(kept_vars) > 0) subset(original_raster, kept_vars) else original_raster[[character(0)]]
  return(list(filtered_raster = filtered_raster, summary = results_summary))
}
