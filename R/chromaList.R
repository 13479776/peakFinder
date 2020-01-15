#' @name chromaProfile
#' @title chromaProfile
#' @description  chromaProfile for merictable function
#' @param xy a object of featureChroma
#' @examples
#' \dontrun{
#' chrs_raw <- chromaProfile()
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export

chromaProfile <- function(xy) {

  #xy is featureChroma objects
  metaID <- list()
  for(i in 1:dim(xy)[1]) {
    feature <- chromaList(xy,i)
    metaID[i] <- list(lapply(split(feature,feature$sampleID),peakgroup))
  }
  return(metaID)

  ####
  #functions
  ####
  peakObj <- setClass("peakObj",slots =
                        c(time = "numeric",
                          sig = "data.frame",
                          area = "character")
  )
  peakgroup <- function(x){
    if(!class(x) == "data.frame") stop("data.frame")
    time <- x$rtime
    sig <- x$intensity
    #if(length(sig[!is.na(sig)]) <= 5) {sig[is.na(sig)] <- -9999999}
    #if(length(sig[!is.na(sig)]) >5 & length(sig[!is.na(sig)]) <= 10) {
    #  sig[is.na(sig)] <- -0.0000099}
    #if(length(sig[!is.na(sig)]) > 10) {
    #  sig[is.na(sig)] <- min(sig,na.rm = TRUE)/2
    #}
    peak.sig <- as.data.frame(cbind(sig,sig))
    area <- as.character(x$name)
    peakgroup <- peakObj(
      time = time,
      sig = peak.sig,
      area = area)
    return(peakgroup)
  }

}

#' @name chromaList
#' @title chromaList
#' @description  chromaList
#' @param x a object of featureChroma
#' @param n ID of feature
#' @examples
#' \dontrun{
#' chrs_raw <- chromaList()
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export

chromaList <-
  function(x, n){
    pk <- x[n,]
    rtint <- list()
    for(i in 1:length(pk)){
      Sname <- pk@phenoData@data$sampleNames[i]
      IDrt <- round(median(c(pk@featureData@data[1,3],pk@featureData@data[1,4])),0)
      IDmz <- round(median(c(pk@featureData@data[1,1],pk@featureData@data[1,2])),0)
      mzID <- paste("M",IDmz, "T",IDrt,sep="")

      rtint[i] <- list(cbind(
        data.frame(cbind(rtime(pk[,i]),intensity(pk[,i]))),
        sampleID = as.character(Sname),
        name = as.character(mzID)))
    }
    rtintO <- do.call(rbind,rtint)
    colnames(rtintO) <- c("rtime","intensity","sampleID","name")
    return(rtintO)
  }


