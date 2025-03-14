% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3dGroupedCircularBoxplot.R
\name{GroupedCircularBoxplot.3D}
\alias{GroupedCircularBoxplot.3D}
\title{Create grouped circular boxplots around the 3D torus surface.}
\usage{
GroupedCircularBoxplot.3D(
  data_in,
  template = "degrees",
  marg = "large",
  shrink = 1,
  H = FALSE,
  constant = "optimal",
  lwd = 2,
  plot_cols = RColorBrewer::brewer.pal(8, "Set2"),
  line_cols = RColorBrewer::brewer.pal(8, "Dark2"),
  arrow_cols = line_cols,
  draw_arrow = TRUE,
  minimal = FALSE,
  scale_widths = FALSE
)
}
\arguments{
\item{data_in}{List of circular objects to be plotted}

\item{template}{One of "degrees" (default), "radians" and "geographics" with eight directions}

\item{marg}{character specifying the baseline 2D circle region as either "large" (default) or "small"}

\item{shrink}{Numeric specifying the factor by which to scale the baseline circular plot. Numbers less than 1 will increase the size.}

\item{H}{Logical indicating if each data point should be drawn outside the hinges of the boxplot}

\item{constant}{Numeric specifying the multiplicative factor determining how far whiskers extend from box. A value of
"optimal" will choose values based on a von Mises distribution (see Buttarazzi et al. 2018)}

\item{lwd}{Scalar indicating the width of the quartile boxplot lines if \code{minimal = TRUE}}

\item{plot_cols}{Vector with the same length as \code{data_in}, specifying the color of the boxplot}

\item{line_cols}{Vector with the same length as \code{data_in}, specifying the color of the median lines}

\item{arrow_cols}{Vector with the same length as \code{data_in}, specifying the color of the arrows pointing to the
medians. Default value is the same as line_cols}

\item{draw_arrow}{Logical specifying if arrows pointing to each median should be drawn}

\item{minimal}{Logical. If true, a quartile boxplot is created with simple colored line and point will replace the box and median line. Radial lines
at fence points will also not be drawn}

\item{scale_widths}{Logical, should the width of each boxplot be scaled based on (the square root of) it's distance from the center?}
}
\description{
Given a named list of circular objects, create a grouped circular boxplot.
}
\examples{
library(circular)
library(GroupedCircularBoxplot)
set.seed(123)
data <- list(
    x = rvonmises(100, circular(pi), 5),
    y = rvonmises(100, circular(pi/2), 2.5),
    z = rvonmises(100, circular(7*pi/4), 8)
)
GroupedCircularBoxplot.3D(data)
}
\author{
Josh Berlinski

Fan Dai
}
