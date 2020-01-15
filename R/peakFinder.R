#' @name peakFinder
#' @title Easy peaks extraction from mass spectroemtry data
#' @description  Peaks extraction from high resulution mass spectrometry data,
#' Optimized parameters are provided for the data produced from UPLC-ORBIRAP,
#' HPLC-ORBITRAP, HPLC-TOF and UPLC-TOF platform.
#' @param files if files is NULL, a file selection dialog appears. Alternatively,
#' a character vectors, containing one of MS data files should be added.
#' @param platformPara The parameters of peaks extraction with centWave in XCMS.
#' A file with list format for in-house parameters.
#' @param bpPARAM xx
#' @param LCPUs  if possible, the number of logical CPUs (if bpPARAM TRUE)
#' @return A peakTarget object. The object includes platformPara information,
#' a peakTable (the output file format of diffreports in xcms package),
#' and a XCMSnExp objects. See the details at https://stattarget.github.io
#' @examples
#' \dontrun{
#' #read the in-house parameter file
#' path <- system.file("parameters",package = "peakFinder")
#' param <- list.files(path,".pk.param.xls",recursive=FALSE,full.names=TRUE)
#' inhouseParam <- readPKparam(param)
#' #perform the peakFinder
#' peakTable <- peakFinder(platformPara = inhouseParam)
#' #peakTable <- peakFinder(platformPara = inhouseParam,bpPARAM = T)
#' xy <- featureChroma(peakTable,expandRt = 5,bpPARAM = T)
#' peakReport(peakTable,fileName = "peakMatrix")
#' chromatogramReport(subdir = "chromatography",object = xy,expandRT = 5)
#' }
#' @keywords peaks, peakFinder, peak_inHouse
#' @export

peakFinder <- function(files = NULL, platformPara,
                       bpPARAM = FALSE,
                       LCPUs = 1L){
  #require(xcms)
  #require(CAMERA)
  #require(BiocParallel)

  cat("\nSelect one of MS data with '.mzML', '.mzXML' '.cdf' or '.CDF' format\n")
  cat("The MS data in the selected folder will be processed","\n\n")
  if(is.null(files)){
    dataPath <- file.choose()

  # show
  if(length(grep(".mzML$",basename(dataPath),fixed =FALSE)) >= 1L){
    cdffiles<-list.files(dirname(dataPath),".mzML$",recursive=FALSE,full.names=TRUE)
  }
  if(length(grep(".mzXML$",basename(dataPath),fixed =FALSE)) >= 1L){
    cdffiles<-list.files(dirname(dataPath),".mzXML$",recursive=FALSE,full.names=TRUE)
  }
  if(length(grep(".CDF$",basename(dataPath),fixed =FALSE)) >= 1L){
    cdffiles<-list.files(dirname(dataPath),".CDF$",recursive=FALSE,full.names=TRUE)
  }
  if(length(grep(".cdf$",basename(dataPath),fixed =FALSE)) >= 1L){
    cdffiles<-list.files(dirname(dataPath),".cdf$",recursive=FALSE,full.names=TRUE)
  }
    } else {

      dataPath = files

      if(length(grep(".mzML$",basename(dataPath),fixed =FALSE)) >= 1L){
        cdffiles<-list.files(dirname(dataPath),".mzML$",recursive=FALSE,full.names=TRUE)
      }
      if(length(grep(".mzXML$",basename(dataPath),fixed =FALSE)) >= 1L){
        cdffiles<-list.files(dirname(dataPath),".mzXML$",recursive=FALSE,full.names=TRUE)
      }
      if(length(grep(".CDF$",basename(dataPath),fixed =FALSE)) >= 1L){
        cdffiles<-list.files(dirname(dataPath),".CDF$",recursive=FALSE,full.names=TRUE)
      }
      if(length(grep(".cdf$",basename(dataPath),fixed =FALSE)) >= 1L){
        cdffiles<-list.files(dirname(dataPath),".cdf$",recursive=FALSE,full.names=TRUE)
      }
      #cdffiles <- cdffiles
      }

  if(length(cdffiles) == 0) {
    stop("No MS Data Found!")}
  #path_dir <- dirname(cdffiles)[1]
  cat(paste("MS Data Found:", cdffiles), sep = "\n")

  cat(paste("\nMS Data reading","( n = ",length(cdffiles),")"), sep = "\n")

  MSdat <- MSnbase::readMSData(cdffiles, mode = "onDisk",msLevel. = 1L)

  if(class(platformPara) == "list"){

    #cat(paste("\nplatformPara selected: 'in-house'"), sep = "\n")

    targetPara <- platformPara
  } else { stop("No platformParam found!")}


  cat("\nPeaks picking...\n")
  #BiocParallel
  if(bpPARAM == TRUE) {

    #core <- detectCores(logical = FALSE)

    nThreads <- parallel::detectCores(logical = TRUE)

    if((nThreads) <= 0) stop("cluster not started; no workers specified; try to decrease saveLCPUs ")

    if(LCPUs >= nThreads) stop(paste("The LCPUs is too large! Please do not more than ",nThreads, " threads\n",sep = ""))


    if(is_unix() | is_macos()) {
      mParam <- BiocParallel::MulticoreParam(LCPUs)
      BiocParallel::register(BiocParallel::bpstart(mParam))
      cat("Unix & IOS Sys\n")
    }
    if(is_windows()) {
      mParam <- BiocParallel::SnowParam(workers = LCPUs, type = "SOCK")
      BiocParallel::register(BiocParallel::bpstart(mParam))
      cat("Windows Sys\n") }
  } else {
    mParam <- BiocParallel::SerialParam()
    #mParam <- bpparam()
    BiocParallel::register(BiocParallel::bpstart(mParam), default=TRUE)
  }
  #register(mParam)
  paraCentWave <- xcms::CentWaveParam(
    ppm=targetPara$ppm, # key param
    peakwidth=targetPara$peakwidth, # key param
    snthresh=targetPara$snthresh,
    prefilter=targetPara$prefilter, # key param
    mzCenterFun=targetPara$mzCenterFun,
    integrate=targetPara$integrate,
    mzdiff=targetPara$mzdiff,
    noise=targetPara$noise, # key param
    fitgauss = as.logical(targetPara$fitgauss),
    verboseColumns = as.logical(targetPara$verboseColumns)

  )
  MSdat <- xcms::findChromPeaks(MSdat,  param = paraCentWave,
                                BPPARAM = mParam,
                                return.type = "XCMSnExp",
                                msLevel = 1L
  )
  MSdat <- xcms::adjustRtime(MSdat, param = xcms::ObiwarpParam(
    binSize = targetPara$adjustRtime.binSize,
  ), msLevel = 1L
  )
  MSdat <- xcms::groupChromPeaks(MSdat, param = xcms::PeakDensityParam(
    binSize = targetPara$adjustRtime.binSize,
    sampleGroups = rep(1, length(fileNames(MSdat))),
    bw=targetPara$group.bw,
    #mzwid=targetPara$group.mzwid,
    minFraction=targetPara$group.minfrac,
    minSamples=targetPara$group.minsamp))

  MSdat <- fillChromPeaks(MSdat, param = FillChromPeaksParam(
    expandMz = 0.5,
    expandRt = 10,
    ppm = 0,
    fixedMz = 0,
    fixedRt = 0)
    ,
    BPPARAM = mParam
  )


  groupval2 <- xcms::featureValues(MSdat, method = c("medret"), value = "into", intensity = "into", filled = TRUE)
  #
  pks <- xcms::chromPeaks(MSdat)

  mat <- do.call(rbind, lapply(xcms::featureDefinitions(MSdat)$peakidx, function(z) {
    pks_current <- pks[z, , drop = FALSE]
    c(pks_current[, c("rt")][[1]]/60L,
      range(pks_current[, c("rtmin", "rtmax")]/60L),
      pks_current[, c("mz")][[1]],
      range(pks_current[, c("mzmin", "mzmax")]),
      pks_current[, c("sn")][[1]]
    )
  }))
  colnames(mat) <- c("rtmed","rtmin", "rtmax","mzmed", "mzmin", "mzmax","sn")
  mat <- as.data.frame(mat,stringsAsFactors= FALSE)

  IDrt <- round((mat[,1]*60),0)
  IDmz <- round((mat[,4]),0)
  peaksID <- paste("M",IDmz, "T",IDrt,sep="")

  #peaksID <- xcms::groupnames(MSdat2)
  #rm(MSdat2)
  stat <- matrix(-1L,nrow = length(peaksID), ncol = 3)
  colnames(stat) <- c("fold","tstat","pvalue")
  tsv <- cbind(name = peaksID, stat, round(mat,4), groupval2, stringsAsFactors =FALSE)

  #
  BiocParallel::bpstop(mParam)

  peakTable <- new("peakTarget")
  peakTable@inputPara <- targetPara

  peakTable@peakTable <- tsv
  peakTable@peakChroma <- MSdat
  return(peakTable)
}



