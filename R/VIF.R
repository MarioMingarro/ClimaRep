#' Filter multicollinear variables based on the Variance Inflation Factor (VIF)
#'
#' This function takes a terra SpatRaster object or a data frame and iteratively removes
#' variables with the highest VIF until all VIFs are below a specified threshold.
#'
#' @param x A terra SpatRaster object or a data frame. If it is a SpatRaster, it
#'   will be internally converted to a data frame.
#' @param th Threshold value for the VIF. Variables with a VIF greater than this
#'   value will be removed. The default value is 10.
#'
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{variables}: A character vector with the names of the
#'       remaining variables after filtering.
#'     \item \code{excluded}: A character vector with the names of the
#'       variables that were excluded due to multicollinearity.
#'     \item \code{corMatrix}: A Pearson correlation matrix between the
#'       remaining variables.
#'     \item \code{results}: A data frame with the names of all original
#'       variables and their calculated VIF values (the last calculated value
#'       before exclusion, if any).
#'   }
#'
#' @examples
#' # Example with a data frame (you can create one to test)
#' df <- data.frame(
#'   var1 = rnorm(100),
#'   var2 = 2 * rnorm(100) + 0.5 * df$var1,
#'   var3 = rnorm(100) + 0.2 * df$var1,
#'   var4 = rnorm(100)
#' )
#' result_df <- vif_filter(df, th = 5)
#' print(result_df)
#'
#' # Example with a terra SpatRaster object (you need to have a SpatRaster object in your environment)
#' \dontrun{
#'   library(terra)
#'   r <- rast(matrix(rnorm(100), 10, 10))
#'   s <- rast(list(r, 2*r + 0.1, r + 0.5, -r))
#'   names(s) <- c("layer1", "layer2", "layer3", "layer4")
#'   result_raster <- vif_filter(s, th = 5)
#'   print(result_raster)
#' }
#'
#' @importFrom terra as.data.frame
#' @importFrom stats as.formula lm summary cor
#'
#' @export
vif_filter <- function(x, th = 10) {

  # Function to calculate the VIF of all variables
  calc_vif <- function(df) {
    vif_values <- sapply(1:ncol(df), function(i) {
      formula <- as.formula(paste(names(df)[i], "~ ."))  # Formula of lm with all other variables
      model <- lm(formula, data = df)
      return(1 / (1 - summary(model)$r.squared))
    })
    names(vif_values) <- colnames(df)
    return(vif_values)
  }

  # Convert the terra SpatRaster object to a data frame
  if (inherits(x, 'SpatRaster')) {
    x <- terra::as.data.frame(x, na.rm = TRUE)
  }


  # Iteratively remove multicollinear variables
  exc <- character(0)  # List of excluded variables
  while (TRUE) {
    v <- calc_vif(x)  # Calculate the VIF
    if (max(v) < th) break  # End if no VIF is above the threshold
    ex <- names(v)[which.max(v)]  # Variable with the highest VIF
    exc <- c(exc, ex)  # Add to the list of excluded variables
    x <- x[, !(colnames(x) %in% ex)]  # Remove variable
  }

  # Create list of results
  result <- list(
    variables = colnames(x),
    excluded = exc,
    corMatrix = cor(x, method = "pearson"),
    results = data.frame(Variables = names(v), VIF = v)
  )


  print(result)
  return(result)
}
