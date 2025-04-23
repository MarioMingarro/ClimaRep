#' Filter multicollinear variables based on Variance Inflation Factor (VIF)
#'
#' This function iteratively removes variables with the highest VIF from a
#' SpatRaster or data frame until all remaining VIFs are below a threshold.
#' It provides a detailed summary of the filtering process, including accurate
#' VIF values for all original variables. If the input is a SpatRaster and
#' an export directory is specified, the selected layers are saved individually.
#'
#' @param x A terra SpatRaster object or a data frame. If SpatRaster, it's
#'   converted to data frame for VIF calculation.
#' @param th Threshold value for the VIF. Variables with VIF >= th are removed.
#'   Default is 10.
#' @param export_dir Optional directory path to save selected raster layers
#'   as GeoTIFFs if input `x` is a SpatRaster. Default is NULL (no export).
#'
#' @return A list with two main components:
#'   \itemize{
#'     \item \code{filtered_object}: The input object (`SpatRaster` or `data.frame`)
#'           containing only the selected variables. Returns NULL if input was
#'           SpatRaster and no variables were selected.
#'     \item \code{summary_results}: A list containing:
#'           \itemize{
#'             \item \code{variables}: Character vector of selected variable names.
#'             \item \code{excluded}: Character vector of excluded variable names.
#'             \item \code{correlation_matrix}: Pearson correlation matrix of
#'                   selected variables (NULL if 1 or 0 selected).
#'             \item \code{vif_results}: Data frame of all original variables
#'                   and their final VIFs when removed or process stopped.
#'           \item \code{threshold}: The VIF threshold used.
#'           }
#'   }
#' @details The function calculates VIF for all variables in the current set.
#'   If the maximum VIF is above the threshold, the variable with the highest
#'   VIF is removed. This process repeats until all remaining VIFs are below
#'   the threshold or only one variable remains. VIFs are calculated based on
#'   a data frame derived from the input, removing rows with NA in any relevant column.
#'
#' @importFrom terra as.data.frame ncell rast subset writeRaster
#' @importFrom stats as.formula lm summary cor na.omit
#' @importFrom utils capture.output
#' @importFrom methods inherits
#'
#' @export
vif_filter <- function(x, th = 10, export_dir = NULL) {

  # --- 1. Preparación ---
  is_raster_input <- inherits(x, 'SpatRaster')
  original_input <- x # Keep a copy of the original object

  # Convert to data frame and remove NA rows for VIF calculation
  if (is_raster_input) {
    # Use na.rm = TRUE to remove cells with NA in any layer relevant to this analysis
    df_for_vif <- terra::as.data.frame(x, na.rm = TRUE)
  } else {
    # Assume data frame input is ready, but remove NA rows just in case
    df_for_vif <- na.omit(x)
  }

  # Handle case with no valid data
  if (nrow(df_for_vif) == 0 || ncol(df_for_vif) == 0) {
    warning("No complete cases or no variables found for VIF calculation. Returning original object with empty results.")
    summary_results <- list(
      variables = character(0),
      excluded = if(ncol(original_input)>0) names(original_input) else colnames(original_input),
      correlation_matrix = NULL,
      vif_results = data.frame(Variables = if(ncol(original_input)>0) names(original_input) else colnames(original_input), VIF = NA),
      threshold = th
    )
    print(summary_results) # Print summary even on error
    return(list(filtered_object = x, summary_results = summary_results))
  }

  # Store original variable names
  original_vars <- colnames(df_for_vif)
  df_filtered <- df_for_vif # Working copy for filtering

  # Check if enough data points for initial VIF calculation
  if (nrow(df_filtered) < ncol(df_filtered) + 1 && ncol(df_filtered) > 1) {
    warning(paste0("Fewer data points (", nrow(df_filtered), ") than variables (", ncol(df_filtered), "+1) for initial VIF calculation. VIFs may be unreliable or infinite. Proceeding with caution."))
  }


  # --- 2. Cálculo y Filtrado Iterativo ---

  calc_vif <- function(df) {
    if (ncol(df) <= 1) { # VIF not calculated for 0 or 1 variable
      return(numeric(0))
    }
    # Need at least k+1 observations for k predictors
    if (nrow(df) < ncol(df) + 1) {
      # Return Inf for all variables if not enough data for model fitting
      return(setNames(rep(Inf, ncol(df)), colnames(df)))
    }

    vif_values <- sapply(1:ncol(df), function(i) {
      formula_str <- paste(names(df)[i], "~ .")
      # Use capture.output to suppress potential warnings from lm model fitting itself
      # Suppress also warning about perfect fit (R-squared=1) which makes VIF Inf
      model_fit <- tryCatch({
        suppressWarnings(lm(as.formula(formula_str), data = df))
      }, error = function(e) {
        # If lm itself errors out (e.g., perfect collinearity that wasn't caught)
        return(NULL)
      })

      if (is.null(model_fit) || !inherits(model_fit, "lm")) {
        return(Inf) # Indicate perfect collinearity or model failure
      }

      r_squared <- summary(model_fit)$r.squared

      if (is.na(r_squared) || r_squared >= 1) {
        return(Inf) # VIF is infinite if R-squared is 1
      }

      vif <- 1 / (1 - r_squared)
      return(vif)
    })
    names(vif_values) <- colnames(df)
    return(vif_values)
  }

  excluded_vars <- character(0)
  all_vif_at_decision <- list() # Store VIF for each var when decision is made

  # Continue filtering as long as there are variables and max VIF >= threshold
  while (ncol(df_filtered) > 1) {
    v <- calc_vif(df_filtered)

    # If no VIFs calculated (<=1 var) or all are finite and below threshold, stop
    if (length(v) == 0 || (all(is.finite(v)) && max(v, na.rm = TRUE) < th)) {
      # Store VIFs for the remaining variables
      if(length(v) > 0) all_vif_at_decision[names(v)] <- v
      break # Exit loop
    }

    # Find variable with max VIF (prioritize Inf)
    if (any(is.infinite(v))) {
      var_to_remove <- names(v)[which(is.infinite(v))[1]] # Take the first one if multiple Inf
      all_vif_at_decision[[var_to_remove]] <- Inf # Store Inf VIF
    } else {
      var_to_remove <- names(v)[which.max(v)]
      all_vif_at_decision[[var_to_remove]] <- v[var_to_remove] # Store the max VIF
    }

    # Remove the variable
    excluded_vars <- c(excluded_vars, var_to_remove)
    df_filtered <- df_filtered[, !(colnames(df_filtered) %in% var_to_remove), drop = FALSE]

    # If only one variable remains after removal, store its VIF (0) and break
    if (ncol(df_filtered) == 1) {
      all_vif_at_decision[[colnames(df_filtered)]] <- 0 # VIF for single variable in this context
      break
    }
  }

  # Handle case where loop broke with 0 or 1 variable initially
  if (ncol(df_filtered) <= 1 && !(colnames(df_filtered) %in% names(all_vif_at_decision)) && ncol(df_filtered) > 0) {
    all_vif_at_decision[[colnames(df_filtered)]] <- 0 # Assign 0 if it was the only variable and loop didn't run/store
  }


  # --- 3. Recopilar Resultados ---

  selected_vars <- colnames(df_filtered)

  # Ensure all original variables have an entry in all_vif_at_decision
  # Variables that were not processed by calc_vif (e.g., only 1 initially) might be missing.
  # Variables perfectly collinear might have caused early errors.
  # Let's create the final VIF df ensuring all original names are present.
  final_vif_df <- data.frame(Variables = original_vars)
  final_vif_df$VIF <- sapply(original_vars, function(var) {
    # Get the stored VIF, default to NA if not found (shouldn't happen if logic is perfect, but for safety)
    val <- all_vif_at_decision[[var]]
    if (is.null(val)) NA else val
  })
  # Ensure order matches original_vars
  final_vif_df <- final_vif_df[match(original_vars, final_vif_df$Variables), ]


  correlation_matrix <- if (ncol(df_filtered) > 1) cor(df_filtered, method = "pearson") else NULL

  summary_results <- list(
    variables = selected_vars,
    excluded = excluded_vars,
    correlation_matrix = correlation_matrix,
    vif_results = final_vif_df, # This is the accurate VIF data frame
    threshold = th
  )

  # --- 4. Exportar Rásteres Seleccionados (si aplica) ---

  filtered_object <- NULL # Will store the filtered SpatRaster or data.frame

  if (is_raster_input) {
    if (length(selected_vars) > 0) {
      # Subset the original raster stack using the selected layer names
      filtered_object <- terra::subset(original_input, selected_vars)

      # Export selected rasters if export_dir is specified
      if (!is.null(export_dir)) {
        message("Exporting selected raster layers to: ", export_dir)
        # Create directory if it doesn't exist
        if (!dir.exists(export_dir)) {
          dir.create(export_dir, recursive = TRUE)
        }

        # Loop through layers and export each one
        export_success <- TRUE
        for (var_name in names(filtered_object)) {
          file_path <- file.path(export_dir, paste0(var_name, ".tif"))
          single_layer <- terra::subset(filtered_object, var_name) # Select single layer by name
          tryCatch({
            terra::writeRaster(single_layer, file_path, overwrite = TRUE, gdal=c("COMPRESS=LZW"))
          }, error = function(e) {
            warning("Failed to export raster layer '", var_name, "' to '", file_path, "': ", e$message)
            export_success <- FALSE
          })
        }
        if(export_success) message("Export complete.")
      }
    } else {
      # If no variables selected, filtered_object remains NULL (as initialized)
      warning("No variables selected after VIF filtering. No raster object returned.")
    }

  } else { # Input was a data frame
    filtered_object <- df_filtered # df_filtered already contains the selected columns
    if (length(selected_vars) == 0) {
      warning("No variables selected after VIF filtering. Returning an empty data frame.")
    }
  }

  # --- 5. Imprimir Resumen y Retornar Resultados ---

  cat("\n--- VIF Filtering Summary ---\n")
  print(summary_results)
  cat("-----------------------------\n")


  # Return a list containing the filtered object and the summary results
  return(list(filtered_object = filtered_object, summary_results = summary_results))
}

# --- Cómo usar la función mejorada ---

# 1. Cargar tus datos (SpatRaster o data frame)
# library(terra)
# present_climatic_variables <- rast("path/to/your/raster_stack.tif")
# df_variables <- as.data.frame(present_climatic_variables, na.rm = TRUE)


# 2. Llamar a la función
#    Puedes especificar el umbral VIF (th) y el directorio de exportación (export_dir)
#    Si input es SpatRaster y export_dir es NULL, no exporta.
#    Si input es data.frame, export_dir es ignorado.

# Ejemplo con SpatRaster y exportación:
# resultado <- vif_filter(present_climatic_variables, th = 5, export_dir = "selected_climatic_variables")

# Ejemplo con SpatRaster sin exportación:
# resultado <- vif_filter(present_climatic_variables, th = 5)

# Ejemplo con data frame:
# resultado <- vif_filter(df_variables, th = 5)


# 3. Trabajar con los resultados

# El resultado es una lista con dos componentes:
# resultado$filtered_object  # Este es el SpatRaster o data frame filtrado
# resultado$summary_results  # Esta es la lista con el resumen (variables, excluidas, VIFs, etc.)

# Para plotear el SpatRaster filtrado (si la entrada fue raster y quedaron variables):
# if (inherits(resultado$filtered_object, "SpatRaster")) {
#   plot(resultado$filtered_object)
# }

# Para ver las variables seleccionadas y sus VIFs:
# print(resultado$summary_results$variables)
# print(resultado$summary_results$vif_results)
