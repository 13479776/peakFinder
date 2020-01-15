#' @name ipo.param
#' @title Paramters of peak extraction
#' @description  To read the parameter used in the peakTarget function
#' @param resultPP xx
#' @param resultRG xx
#' @param pf.param xx
#' @examples
#' \dontrun{
#' library(IPO)
#' library(peakFinder)
#' msdata <- list.files(getwd(),"mzML", full.names=TRUE)
#' paramsPP <- getDefaultXcmsSetStartingParams()
#' paramsPP$mzdiff <- -0.001
#' paramsPP$ppm <- c(5,15)
#' paramsPP$min_peakwidth <- c(7,14)
#' paramsPP$max_peakwidth <- c(20,60)
#' paramsPP$noise <- 100000
#' paramsPP$snthresh <- c(2,15)
#' resultPP <- optimizeXcmsSet(msdata[2:3], paramsPP, subdir="test")
#' paramsRG <- getDefaultRetGroupStartingParams()
#' paramsRG$gapInit <- 0.2
#' paramsRG$profStep <- 1
#' paramsRG$minfrac <- 0.75
#' resultRG <- optimizeRetGroup(resultPP$best_settings$xset, paramsRG, nSlaves=2)
#' # ipo----peakFinder
#' inhouseParam <- ipo.param(resultPP,resultRG,pf.param= NULL)
#' peak <- peakFinder(platformPara = inhouseParam,bpPARAM = T)
#' xy <- featureChroma(peak,expandRt = 10,bpPARAM = T)
#' peak.table(peak,fileName = "peakMatrix")
#' peakReport(subdir = "chromatography",object = xy, peakInter = NULL,expandRT = 10)

#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export


ipo.param <- function(resultPP,resultRG,pf.param = NULL) {


  if(is.null(pf.param)){
   path <- system.file("parameters",package = "peakFinder")
   param <- list.files(path,".pk.param.csv",recursive=FALSE,full.names=TRUE)
   inhouseParam <- readPKparam(param)
   pf.param <- inhouseParam
  } else {
    pf.param = pf.param
  }

  #
  ipo.x <- resultPP
  ipo.y <- resultRG

  pf.param$ppm <- ipo.x$best_settings$parameters$ppm
  pf.param$peakwidth <- c(ipo.x$best_settings$parameters$min_peakwidth,ipo.x$best_settings$parameters$max_peakwidth)
  pf.param$snthresh <- ipo.x$best_settings$parameters$snthresh
  pf.param$prefilter <- c(ipo.x$best_settings$parameters$prefilter, ipo.x$best_settings$parameters$value_of_prefilter)
  pf.param$noise <- ipo.x$best_settings$parameters$noise
  pf.param$mzCenterFun <- ipo.x$best_settings$parameters$mzCenterFun
  pf.param$integrate <- ipo.x$best_settings$parameters$integrate
  #
  pf.param$group.bw <- ipo.y$best_settings$bw
  pf.param$group.mzwid <- ipo.y$best_settings$mzwid
  pf.param$group.minfrac <- ipo.y$best_settings$minfrac
  return(pf.param)
}
