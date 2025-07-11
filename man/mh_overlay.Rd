% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mh_overlay.R
\name{mh_overlay}
\alias{mh_overlay}
\title{Overlay Mahalanobis-based Climate Representativeness Classifications}
\usage{
mh_overlay(folder_path)
}
\arguments{
\item{folder_path}{\code{character}. The path to the directory containing the classification rasters (\code{.tif}) generated by \code{mh_rep} or \code{mh_rep_ch}. These rasters should primarily
contain the categories: \code{1} (Retained/Represented), \code{2} (Lost), and \code{3} (Novel).
Category \code{0} (Non-represented) will be ignored for the RGB output.}
}
\value{
Writes the multi-layered (\code{ClimaRep_overlay.tif}) outputs to disk in a new \code{overlay} subfolder within the \code{folder_path}.
When \code{mh_rep_ch} results are used, the output layers consistently represent counts for \strong{Lost} (Red), \strong{Retained} (Green), and \strong{Novel} (Blue) categories across all input rasters. Designed for direct RGB plotting.
When \code{mh_rep} results are used, the output layers consistently represent counts for \strong{Represented} categories across all input rasters.
}
\description{
Combines multiple single-layer rasters (\code{tif}), outputs from \code{mh_rep} or \code{mh_rep_ch} for different input polygons, into a multi-layered \code{SpatRaster}.

This function handles inputs from both \code{mh_rep} (which primarily contains \strong{Represented} areas) and \code{mh_rep_ch} (which includes \strong{Retained}, \strong{Lost}, and \strong{Novel} areas). The output layers consistently represent counts of each input.
}
\details{
This function streamlines the aggregation of Climate Representativeness classifications. It is designed to work with outputs from both \code{mh_rep} and \code{mh_rep_ch}.

For each of the three key categories (Lost, Retained/Represented, Novel), the function:
\enumerate{
\item Identifies and reads all \code{.tif} files within the \code{folder_path}.
\item For each input raster, it creates a binary layer: \code{1} if the cell's value matches the target category (e.g., \code{2} for 'Lost'), and \code{0} otherwise.
\item Sums these binary layers to generate a cumulative count for that specific category at each grid cell.
}

The three resulting count layers (Lost, Retained, Novel) are then consistently stacked in the following order:
\itemize{
\item First layer (Red): Cumulative count of \strong{Lost}.
\item Second layer (Green): Cumulative count of \strong{Retained}.
\item Third layer (Blue): Cumulative count of \strong{Novel}.
}
This fixed order ensures that the output \code{SpatRaster} is immediately ready for direct RGB visualization using \code{terra::plotRGB()}, where the color mixtures will intuitively reflect
the spatial agreement of these change types.

The output \code{SpatRaster} contains raw counts. While \code{terra::plotRGB()} often handles stretching for visualization, users might normalize these counts manually (e.g., to 0-255) for finer control over visual contrast.

A new subfolder named \verb{overlay/} will be created within the \code{folder_path}. The resulting three-layered RGB will be saved as \code{ClimaRep_overlay.tif} inside this new \verb{overlay/} subfolder.
}
\examples{
ClimaRep_overlay <- ClimaRep::mh_overlay(folder_path = system.file("extdata", package = "ClimaRep"))
terra::plotRGB(ClimaRep_overlay)
terra::plot(ClimaRep_overlay)
}
