#' Phases of clock genes
#'
#' Data from an experiment performed on mice to determine compare
#' the peak phase times of circadian clock genes under ad-libitum
#' and time-restricted feeding regimines. The two imported data frames
#' are `alf` and `trf` corresponding to ad-libitum and time-restricted
#' results, respectively. The observations are measured in degrees,
#' corresponding to time on a 24-hour zeigeber time clock cycle.
#'
#' @docType data
#'
#' @usage data(deota)
#'
#' @references Deota S, Lin T, Chaix A, Williams A, Le H, Calligara H,
#' Ramasamy R, Huang L, Panda S (2023) Diurnal transcriptome landscape
#' of a multi-tissue response to time-restricted feeding in mammals.
#' Cell Metabolism 35(1):150-164.e4
#'
#' @source \href{https://doi.org/10.1016/j.cmet.2022.12.006}
"deota"

#' Wind direction
#'
#' Measurements of wind direction from Ames, Iowa, USA and
#' Jamshedpur, India. Data were obtained for over the entire
#' 2023 year. The list contains two data frames, `AMW` and `VEJS`,
#' corresponding to Ames and Jamshedpur, respectively. The Ames data
#' were collected every five minutes, while the Jamshedpur data are
#' hourly.
#'
#' Each datframe contains three columns, the first corresponding to
#' the code of the weather station `name`, the second `valid` which
#' contains the time the mesurement was taken, and `drct`, which is
#' the direction from which the wind was measured in degrees from true
#' north, clockwise.
#'
#' Note that a direction of 0 degrees indicates no measurable wind, and
#' thus does not contain useful information on the wind direction. A value
#' of 360 degrees corresponds to true north.
#'
#' @docType data
#'
#' @usage data(wind_direction)
#'
#'
#' @source Iowa State University's Iowa Environmental Mesonet \href{https://mesonet.agron.iastate.edu/request/download.phtml}
"wind_direction"
