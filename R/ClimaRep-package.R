#' @keywords internal
#' @aliases ClimaRep
"_PACKAGE"

#' @title ClimaRep: Prueba Minima
#' @description
#' ## Overview
#' The primary goal of `ClimaRep` is to quantify how well the climate within specific polygons (`sf`) represents the broader climate space defined by climate variables (`SpatRaster`) within a study area (`sf`).
#' It also provides functions to evaluate how this representativeness changes under projected future climate conditions.
#'
#' ## Key Features
#' The package includes functions for:
#' * Filtering raster climate variables to reduce multicollinearity (`vif_filter`).
#' * Estimating current climate representativeness (`mh_rep`).
#' * Estimating changes in climate representativeness under future climate projections (`mh_rep_ch`).
#' * Estimating climate representativeness overlay (`mh_overlay`).
#'
#' ## More Details
#' https://github.com/MarioMingarro/ClimaRep
#' @name ClimaRep-package
NULL
