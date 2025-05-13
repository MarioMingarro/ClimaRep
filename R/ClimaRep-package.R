#' @keywords internal
#' @docType package
#' @name ClimaRep-package
#' @aliases ClimaRep
#' @title Estimating Climate Representativeness
#' @description
#' The `ClimaRep` package offers tools to estimate the climate representativeness of defined areas and quantifies and analyzes its transformation under future climate change scenarios.
#' @details
#' ## Overview
#' The primary goal of `ClimaRep` is to quantify how well the climate within specific polygons (`sf`) represents the broader climate space defined by climate variables (`SpatRaster`) within a study area (`sf`).
#' It also provides functions to evaluate how this representativeness changes under projected future climate conditions.
#'
#' ## Key Features
#' The package includes functions for:
#' * Filtering raster climate variables to reduce multicollinearity (`\code{\link{vif_filter}}`).
#' * Estimating current climate representativeness (`\code{\link{mh_rep}}`).
#' * Estimating changes in climate representativeness under future climate projections (`\code{\link{mh_rep_ch}}`).
#' * Estimating climate representativeness overlay (`\code{\link{mh_overlay}}`).
#'
#' ## More Details
#' https://github.com/MarioMingarro/ClimaRep
#'
