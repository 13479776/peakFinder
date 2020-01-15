#' @name readPKparam
#' @title Paramters of peak extraction
#' @description  To read the parameter used in the peakTarget function
#' @param x a template with '.pk.param.xls' format. if files is NULL, a file selection dialog appears.
#' @examples
#' \dontrun{
#' readPKparam()
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export
readPKparam <- function(x = NULL){

  if(is.null(x)){
    dataPath <- file.choose()
    Para <- read.table(dataPath,header = TRUE,stringsAsFactors= FALSE,sep="\t")
    } else {
    Para <- read.table(x,header = TRUE,stringsAsFactors= FALSE)
    }
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
  targetPara$fitgauss <- Para$fitgauss[1]
  targetPara$verboseColumns <- Para$verboseColumns[1]
  targetPara$method <- Para$method[1]
  targetPara$mzCenterFun <- Para$mzCenterFun[1]
  return(targetPara)
}

#' @name writePKparam
#' @title Paramters of peak extraction
#' @description  To read the parameter used in the peakTarget function
#' @param object a template with '.pk.param.xls' format
#' @param paramfile xx
#' @examples
#' \dontrun{
#' writePKparam()
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export

writePKparam <- function(object,paramfile){
  if(class(object) == "peakTarget") {
  x <- object@inputPara
  temp <- as.data.frame(do.call(rbind,x), stringsAsFactors=FALSE)
  temp <- cbind(rownames(temp),temp)
  colnames(temp) <- c("id","start","end")
  temp$end <- "NULL"
  temp["peakwidth","end"] <- x$peakwidth[2]
  temp["prefilter","end"] <- x$prefilter[2]

  write.table(temp, paste(paramfile,".pk.param.xls",sep = ""), row.names = FALSE,quote = FALSE, sep = "\t")

  } else {stop("peakTarget object is required!")}
}
