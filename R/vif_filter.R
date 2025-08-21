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
#' A high VIF value suggests a high degree of collinearity with other predictors (values exceeding `5` or `10` are often considered problematic; see O'Brien, 2007).
#'
#' Key steps:
#' \enumerate{
#' \item Validate inputs: Ensures `x` is a `SpatRaster` with at least two layers and `th` is a valid `numeric` value.
#' \item Convert the input `SpatRaster` (`x`) to a `data.frame` and remove rows with missing values (`NA`).
#' \item **Iterative Filtering Process:**
#' \itemize{
#' \item In each step, the VIF for each remaining variable is calculated by fitting a linear regression model where that variable is the response and all other remaining variables are the predictors. The VIF is then computed from the model's R-squared value using the formula \eqn{VIF = 1 / (1 - R^2)}.
#' \item The variable with the highest VIF is identified.
#' \item If this highest VIF value is greater than the threshold (`th`), that variable is removed.
#' \item This process repeats until the highest VIF among the remaining variables is less than or equal to `th`, or until only one variable remains.
#' }
#' \item The function handles cases of perfect collinearity or variables with zero variance by prioritizing them for removal, preventing numerical instability.
#' }
#'The function returns a `list` with two main components:
#' \itemize{
#' \item `filtered_raster`: A `SpatRaster` object containing only the layers that were kept.
#' \item `summary`: A list containing the original Pearson's correlation matrix, the names of the kept and excluded variables, and the final VIF values for the retained variables.
#' }
#'The internal VIF calculation includes checks to handle potential numerical instability, such as columns with zero or near-zero variance and cases of perfect collinearity among variables,
#'which could otherwise lead to errors (e.g., infinite VIFs). Variables identified as having infinite VIF due to perfect collinearity are prioritized for removal.
#'
#' References:
#' O’brien (2007) A Caution Regarding Rules of Thumb for Variance Inflation Factors. Quality & Quantity, 41: 673–690. doi:10.1007/s11135-006-9018-6
#' Legendre & Legendre (2012) Numerical ecology (3rd ed.). Elsevier.
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
  if (ncol(x_df) <= 1) {
    warning("Less than two variables present. VIF filtering not applicable.")
    return(list(
      filtered_raster = original_raster,
      summary = list(
        kept_layers = names(original_raster),
        excluded_layers = character(0),
        original_correlation_matrix = "Correlation matrix not applicable (less than 2 variables).",
        final_vif_values = "Only one variable kept. VIF calculation not applicable."
      )
    ))
  }
  calc_vif_matrix_method <- function(df) {
    if (ncol(df) <= 1) return(numeric(0))
    cor_matrix <- cor(df, use = "pairwise.complete.obs")
    if (det(cor_matrix) == 0) {
      warning("Perfect collinearity detected. Cannot invert the correlation matrix.")
      return(rep(Inf, ncol(df)))
    }
    vif_values <- diag(solve(cor_matrix))
    names(vif_values) <- colnames(df)
    return(vif_values)
  }
  kept_vars <- colnames(x_df)
  excluded_vars <- character(0)
  while (length(kept_vars) > 1) {
    df_subset <- x_df[, kept_vars, drop = FALSE]
    v <- calc_vif_matrix_method(df_subset)
    if (all(v < th) || any(is.infinite(v))) {
      break
    }
    max_vif_name <- names(v)[which.max(v)]
    kept_vars <- setdiff(kept_vars, max_vif_name)
    excluded_vars <- c(excluded_vars, max_vif_name)
  }
  final_vif_values <- NULL
  if (length(kept_vars) > 1) {
    final_vif_values <- round(calc_vif_matrix_method(x_df[, kept_vars, drop = FALSE]), 4)
  } else if (length(kept_vars) == 1) {
    final_vif_values <- "Only one variable kept. VIF calculation not applicable."
  } else {
    final_vif_values <- "No variables kept."
  }
  original_cor_matrix <- round(cor(x_df, method = "pearson"), 4)
  results_summary <- list(
    kept_layers = kept_vars,
    excluded_layers = excluded_vars,
    original_correlation_matrix = original_cor_matrix,
    final_vif_values = final_vif_values)
  filtered_raster <- if (length(kept_vars) > 0) subset(original_raster, kept_vars) else original_raster[[character(0)]]
  message("All processes were completed")
  return(list(filtered_raster = filtered_raster, summary = results_summary))
}
