#' @name chromatogramReport
#' @title chromatogramReport
#' @description  peakReport
#' @param subdir directory
#' @param object a object of featureChroma
#' @param peakInter xx
#' @param expandRT rt
#' @param bpPARAM xx
#' @param LCPUs the number of logical CPU cores on the current host.
#' @examples
#' \dontrun{
#' peakReport(subdir, object, expandRT = 5)
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export

chromatogramReport <- function(subdir, object, peakInter = NULL, expandRT = 30,bpPARAM = FALSE,LCPUs = 1L) {

  if(!class(object) %in% c("Chromatograms","XChromatograms")) {stop("Chromatograms or XChromatograms format only")}


  if(bpPARAM == TRUE) {

    #core <- detectCores(logical = FALSE)

    nThreads <- parallel::detectCores(logical = TRUE)

    if((nThreads) <= 0) stop("cluster not started; no workers specified; try to decrease saveLCPUs ")

    if(LCPUs >= nThreads) stop(paste("The LCPUs is too large! Please do not more than ",nThreads, " threads\n",sep = ""))


    if(is_macos()) {
      mParam <- BiocParallel::MulticoreParam(LCPUs)
      BiocParallel::register(BiocParallel::bpstart(mParam))
      cat("ios done\n")
    }
    if(is_windows()) {
      mParam <- BiocParallel::SnowParam(workers = LCPUs, type = "SOCK")
      BiocParallel::register(BiocParallel::bpstart(mParam))
      cat("win done\n") }
  } else {
    mParam <- BiocParallel::SerialParam()
    BiocParallel::register(BiocParallel::bpstart(mParam), default=TRUE)
  }


if(is.null(peakInter)) {

  dir.create(subdir)
  path <- paste(getwd(),"/",subdir,sep = "")
  fok <- 1:dim(xy)[1]
  #message("Total:",length(fok)," peaks")
  yyyy <- BiocParallel::bplapply(fok,po,object = object,pathok = path,expandRT= expandRT,BPPARAM=mParam)
  BiocParallel::bpstop(mParam)
} else {

   peakInter <- peakInter
   ok <- grep("TRUE",peakInter$jugx)
   fok <- grep("FALSE",peakInter$jugx)
# set directory
dir.create(subdir)
path <- paste(getwd(),"/",subdir,sep = "")
dir.create(paste(path,"/","authentic",sep = ""))
pathok <- paste(path,"/","authentic",sep = "")
dir.create(paste(path,"/","removed",sep = ""))
pathdok <-paste(path,"/","removed",sep = "")

#
cat(paste(" authentic ", length(ok), " peaks...", sep = ""))
#message("Total:",length(ok)," peaks")
xxxx <- BiocParallel::bplapply(ok,po,object = object,pathok = pathok,expandRT= expandRT,BPPARAM=mParam)
#message("removed peaks...")
cat("\n",paste("removed ", length(fok), " peaks...", sep = ""))
#message("Total:",length(fok)," peaks")
yyyy <- BiocParallel::bplapply(fok,po,object = object,pathok = pathdok,expandRT= expandRT,BPPARAM=mParam)
BiocParallel::bpstop(mParam)

}
  message("finshed!")

}
##
#' @name peakInter
#' @title peakInter
#' @description  peakInter
#' @param object directory
#' @param divg a object of featureChroma
#' @param topIndex xx
#' @examples
#' \dontrun{
#' peakReport(subdir, object, expandRT = 5)
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export
peakInter <- function(object,divg, topIndex) {

  if(!class(object) %in% c("Chromatograms","XChromatograms")) {stop("Chromatograms or XChromatograms  format only")}
  if(!class(divg) == "factor") {stop("divg shou be factor")}

  index <- paste("F",(rep(1:dim(object)[1],each = dim(object)[2])),sep = "")
  if(length(divg) != length(index)) { stop("length problems")}

  com <- as.data.frame(cbind(as.character(divg),index),stringsAsFactors = F)

  #
  #require(dplyr)
  jug <- function(x,topIndex) if("TRUE" %in% c(topIndex %in% x))
  {return(TRUE)} else
  {return(FALSE)}

  com$index <- factor(com$index, as.character(unique(com$index)))
  new <- com %>%
    dplyr::group_by(index) %>%
    dplyr::summarise(counts = n(), hit = paste(V1,collapse="-"),
                     jugx =jug(V1,topIndex= topIndex)
    )

  return(new)
}

po <- function(m,object,pathok,expandRT,leg.pos = "bottom"){
  #plot_list1 = list()
  dat <- chromaList(object, m)
  ID <- levels(dat$name)[1]
  plotChroma(object,m,expandRT = expandRT, facet_wrap = F, leg.pos = leg.pos)

  cat(m,"...")
  ggplot2::ggsave(paste(pathok,"/", ID,"_",m,".pdf",sep = ""),width = 6, height = 5)
}


#' @name peakReport
#' @title peakReport
#' @description  peakReport
#' @param fileName file name
#' @param object a object of peakTarget
#' @return Two tables with .ok.param.tsv and .tsv format
#' @examples
#' \dontrun{
#' peak.table(object,fileName)
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export
peakReport <- function(object,fileName) {

  if(!class(object) == "peakTarget") stop("peakTarget class is required!")
  cot <- object@peakTable
  cot1 <- cot[,1:11]
  cot2 <- as.matrix(cot[,12:dim(cot)[2]])
  meanc <- apply(cot2,1, mean, na.rm = TRUE)
  sdc <- apply(cot2,1, sd, na.rm = TRUE)
  cotnew <- data.frame(cot1,meanArea=meanc,sdArea=sdc,cot2)


  table <- write.table(cotnew, paste(fileName,".tsv",sep = ""), quote = FALSE, sep = "\t", col.names = NA)
  param <- writePKparam(object,fileName)
  cat("peak.table:: Peak Matrix was exported!")
}
