#' @name peak_inHouse
#' @title 'in-house' paramters of peak extraction
#' @description  To edit the parameter of 'in-house' in the peakTarget function
#' @examples
#' \dontrun{
#' peak_inHouse()
#' }
#' @keywords peaks, peakTarget, peak_inHouse
#' @export

readPKparam <- function(x){
  Para <- read.csv("inhouse.pk.param.csv",header = TRUE,stringsAsFactors= FALSE)
  #Para <- t(Para)
  rownames(Para) <- Para[,1]; Para <- t(Para)
  Para <- as.data.frame(Para[-1,],stringsAsFactors= FALSE)
  targetPara <- list()
  targetPara$ppm <- Para$ppm[1]
  targetPara$peakwidth <- c(Para$peakwidth[1],Para$peakwidth[2])
  targetPara$snthresh <- Para$snthresh[1]
  targetPara$prefilter <- c(Para$prefilter[1],Para$prefilter[2])
  targetPara$integrate <- Para$integrate[1]
  targetPara$mzdiff <- Para$mzdiff[1]
  targetPara$noise <- Para$noise[1]
  targetPara$adjustRtime.binSize <- Para$adjustRtime.binSize[1]
  targetPara$group.bw <- Para$group.bw[1]
  targetPara$group.mzwid <- Para$group.mzwid[1]
  targetPara$group.minfrac <- Para$group.minfrac[1]
  targetPara$group.minsamp <- Para$group.minsamp[1]
  targetPara$group.binsize <- Para$group.binsize[1]
  targetPara <- lapply(targetPara, as.numeric)
  targetPara$method <- Para$method[1]
  targetPara$mzCenterFun <- Para$mzCenterFun[1]



  paralink <- system.file('data',package = 'peakTarget')
  save(targetPara,file = paste(paralink,"/inhouse.RData",sep = ""))
}



if(class(platformPara) == "character") {

  cat(paste("\nplatformPara selected:", platformPara), sep = "\n")


  if(!platformPara %in% c("UOrbi","uorbi","HOrbi","horbi",
                          "UTOF","utof","HTOF","htof")){
    stop(platformPara," ?? No platformPara found!!!")
  }

  if(platformPara %in% c("UOrbi","uorbi"))

  {
    data("xcms.UPLC.Orbi",package = "peakTarget")
  }

  if(platformPara %in% c("HOrbi","horbi"))

  {
    data("xcms.HPLC.Orbi",package = "peakTarget")
  }
  if(platformPara %in% c("UTOF","utof"))

  {
    data("xcms.UPLC.TOF",package = "peakTarget")
  }

  if(platformPara %in% c("HTOF","htof"))

  {
    data("xcms.HPLC.TOF",package = "peakTarget")
  }
  targetPara
}

