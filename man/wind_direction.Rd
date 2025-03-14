% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{wind_direction}
\alias{wind_direction}
\title{Wind direction}
\format{
An object of class \code{list} of length 2.
}
\source{
Iowa State University's Iowa Environmental Mesonet \href{https://mesonet.agron.iastate.edu/request/download.phtml}
}
\usage{
data(wind_direction)
}
\description{
Measurements of wind direction from Ames, Iowa, USA and
Jamshedpur, India. Data were obtained for over the entire
2023 year. The list contains two data frames, \code{AMW} and \code{VEJS},
corresponding to Ames and Jamshedpur, respectively. The Ames data
were collected every five minutes, while the Jamshedpur data are
half-hourly.
}
\details{
Each datframe contains three columns, the first corresponding to
the code of the weather station \code{name}, the second \code{valid} which
contains the time the mesurement was taken, and \code{drct}, which is
the direction from which the wind was measured in degrees from true
north, clockwise.

Note that a direction of 0 degrees indicates no measurable wind, and
thus does not contain useful information on the wind direction. A value
of 360 degrees corresponds to true north.
}
\keyword{datasets}
