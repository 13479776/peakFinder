#' parameters
#'
#' A dataset containing peak picking parameters.
#'
#' The contents are as follows:
#' \itemize{
#'   \item method.
#'   \item ppm.
#'   \item peakwidth.
#'   \item snthresh.
#'   \item prefilter.
#'   \item mzCenterFun.
#'   \item intergrate.
#'   \item noise.
#'   \item adjustTtime
#'   \item group.bw
#'   \item group.mzwid
#'   \item group.minfac
#'   \item group.minsamp
#'   \item group.binsize
#' }




#' @importFrom utils read.csv
#' @importFrom utils capture.output
#' @importFrom utils data
#' @importFrom utils capture.output
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom stats median
#' @importFrom methods new
#' @importFrom methods setClass
#' @importFrom graphics plot
#' @importFrom methods validObject
#' @importFrom magrittr %>%
#' @importFrom MSnbase centroided
#' @importFrom signal butter
#'
.onAttach <- function(...) {
  packageStartupMessage("\nThis is 'peakFinder version 0.0.1' ")
  # statTarget::statTargetGUI()
}



