#' @title Filtra capas de un SpatRaster basadas en el Factor de Inflación de Varianza (VIF)
#'
#' @description
#' Esta función filtra iterativamente las capas de un objeto \code{SpatRaster} eliminando aquellas con el VIF más alto que superan un umbral especificado, hasta que todas las capas restantes tengan un VIF por debajo del umbral. El cálculo de VIF se realiza después de eliminar filas con valores \code{NA}.
#'
#' @param x Un objeto \code{SpatRaster} que contiene las capas a filtrar. Debe ser un objeto \code{SpatRaster} válido.
#' @param th Un valor numérico que especifica el umbral del VIF. Las capas con un VIF superior a este umbral son candidatas a ser eliminadas. El valor por defecto es 10.
#'
#' @details
#' La función convierte el \code{SpatRaster} a un \code{data.frame}, eliminando cualquier fila que contenga \code{NA} en *cualquiera* de las columnas antes de calcular el VIF. El proceso de filtrado es iterativo: en cada paso, se calcula el VIF para las variables restantes y se elimina la variable con el VIF más alto si éste excede el umbral \code{th}. Este proceso continúa hasta que todas las variables restantes tengan un VIF por debajo del umbral o solo quede una variable. La función imprime un resumen del proceso en la consola, incluyendo las capas conservadas y excluidas, la matriz de correlación original (calculada después de eliminar NAs) y los VIF finales de las capas conservadas. La función interna para calcular VIF incluye comprobaciones para manejar columnas con varianza cero o casi cero y casos de colinealidad perfecta.
#'
#' @return
#' Un objeto \code{SpatRaster} que contiene solo las capas que pasaron el filtro del VIF. Si todas las capas son excluidas durante el proceso, devuelve un \code{SpatRaster} vacío.
#'
#' @importFrom terra SpatRaster as.data.frame subset
#'
#' @examples
#' # requires terra package
#' library(terra)
#'
#'
#'
#'
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
