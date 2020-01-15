#' @name peakTarget
#' @title Easy peaks extraction from mass spectroemtry data
#' @description  Peaks extraction from high resulution mass spectrometry data,
#' Optimized parameters are provided for the data produced from UPLC-ORBIRAP,
#' HPLC-ORBITRAP, HPLC-TOF and UPLC-TOF platform.
#' @param files if files is NULL, a file selection dialog appears. Alternatively,
#' a character vectors, containing one of MS data files should be added.
#' @param platformPara the parameters of peaks extraction with centWave in XCMS,
#' i.e. 'uorbi' for UPLC-ORBITRAP; 'horbi' for HPLC-ORBITRAP; 'utof' for UPLC-TOF;
#' 'htof' for HPLC-TOF; 'in-house' for in-house setting.
#' @param annotation Annotation of peaks or not.
#' @param isotopeRemove Removal of isotope peaks or not.
#' @param polarity 'positive' or 'negative' ion mode .
#' @return A peakTarget object. The object includes platformPara information,
#' a peakTable (the output file format of diffreports in xcms package),
#' and a XCMSnExp objects. See the details at https://stattarget.github.io
#' @examples
#' \dontrun{
#' peakTable <- peakTarget(platformPara = "uorbi", polarity = "positive")
#' }
#' @keywords peaks, peakTarget, peak_inHouse
#' @export

peakTarget <- function(files = NULL, platformPara, polarity="positive",
                       annotation=TRUE, isotopeRemove=TRUE, bpPARAM = TRUE){
    #require(xcms)
    #require(CAMERA)
    #require(BiocParallel)

    if(!platformPara %in% c("UOrbi","uorbi","HOrbi","horbi",
                            "UTOF","utof","HTOF","htof","in-house","IH","ih")){
      stop(platformPara," ?? No platformPara found!!!")
    }

    cat("\nSelect one of MS data with '.mzML', '.mzXML' and '.CDF' format\n")
    cat("The MS data in the selected folder will be processed","\n\n")
    if(is.null(files)){
      dataPath <- file.choose()
    } else {
      dataPath = files }

    # show
    cdffiles<-list.files(dirname(dataPath),".mzXML",recursive=FALSE,full.names=TRUE)
    if(length(cdffiles) == 0 ) {
      cdffiles<-list.files(dirname(dataPath),".mzML",recursive=FALSE,full.names=TRUE)
    }
    if(length(cdffiles) == 0 ) {
      cdffiles<-list.files(dirname(dataPath),".CDF",recursive=FALSE,full.names=TRUE)
    }
    if(length(cdffiles) == 0 ) {
      stop("None MS Data Found!")
    }

    cat(paste("MS Data Found:", cdffiles), sep = "\n")
    MSdat <- readMSData(cdffiles, mode = "onDisk")

    if(platformPara %in% c("UOrbi","uorbi"))

    {
      data("xcms.UPLC.Orbi")
    }

    if(platformPara %in% c("HOrbi","horbi"))

    {
      data("xcms.HPLC.Orbi")
    }
    if(platformPara %in% c("UTOF","utof"))

    {
      data("xcms.UPLC.TOF")
    }

    if(platformPara %in% c("HTOF","htof"))

    {
      data("xcms.HPLC.TOF")
    }

    if(platformPara %in% c("in-house","IH","ih"))

    {
      data("inhouse")
    }

    cat(paste("\nplatformPara selected:", platformPara), sep = "\n")


    cat("\nPeaks picking...\n")
    #BiocParallel
    if(bpPARAM == TRUE) {

      #core <- detectCores(logical = FALSE)

      nThreads <- detectCores(logical = TRUE)
      if(is_macos()) {
        mParam <- MulticoreParam(nThreads-1)
        register(bpstart(mParam))
        cat("ios done\n")
      }
      if(is_windows()) {
        mParam <- SnowParam(workers = nThreads-1, type = "SOCK")
        register(bpstart(mParam))
        cat("win done\n") }
    } else {
      mParam <- SerialParam()
      register(mParam, default=TRUE)
      }
    #register(mParam)

    paraCentWave <- CentWaveParam(
      ppm=targetPara$ppm, # key param
      peakwidth=targetPara$peakwidth, # key param
      snthresh=targetPara$snthresh,
      prefilter=targetPara$prefilter, # key param
      mzCenterFun=targetPara$mzCenterFun,
      integrate=targetPara$integrate,
      mzdiff=targetPara$mzdiff,
      noise=targetPara$noise # key param
      )
    MSdat <- findChromPeaks(MSdat,  param = paraCentWave,
      BPPARAM = mParam,
      return.type = "XCMSnExp",
      msLevel = 1L
    )
    MSdat <- adjustRtime(MSdat, param = ObiwarpParam(
                                binSize = targetPara$adjustRtime.binSize,
                                ), msLevel = 1L
                         )
    MSdat <- groupChromPeaks(MSdat, param = PeakDensityParam(
      binSize = targetPara$adjustRtime.binSize,
      sampleGroups = rep(1, length(fileNames(MSdat))),
      bw=targetPara$group.bw,
      #mzwid=targetPara$group.mzwid,
      minFraction=targetPara$group.minfrac,
      minSamples=targetPara$group.minsamp))


    groupval2 <- featureValues(MSdat, method = c("medret"), value = "into", intensity = "into", filled = TRUE)
    #
    pks <- chromPeaks(MSdat)

    mat <- do.call(rbind, lapply(featureDefinitions(MSdat)$peakidx, function(z) {
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
    peaksID <- groupnames(MSdat)
    stat <- matrix(-1L,nrow = length(peaksID), ncol = 3)
    colnames(stat) <- c("fold","tstat","pvalue")
    tsv <- cbind(name = peaksID, stat, round(mat,4), groupval2, stringsAsFactors =FALSE)





    # annotation
    if(annotation){
    cat("Peaks annotation...\n")
    annoMSdat <- suppressMessages(as(MSdat,"xcmsSet"))
    annoMSdat <-CAMERA::xsAnnotate(annoMSdat, polarity = polarity)
    suprr = capture.output(annoMSdat <-CAMERA:: groupFWHM(annoMSdat))
    suprr = capture.output(annoMSdat <- CAMERA::findIsotopesWithValidation(annoMSdat, ppm = targetPara$ppm,
                                            mzabs=abs(targetPara$mzdiff)))
    #annoMSdat <- findAdducts(annoMSdat,polarity=polarity)
    peakTableAnnoMSdat <- CAMERA::getPeaklist(annoMSdat)
    #annoTable <- cbind( peakTableAnnoMSdat[,dim(peakTableAnnoMSdat)[2]-2],
     #                   peakTableAnnoMSdat[,dim(peakTableAnnoMSdat)[2]-1],
     #                   peakTableAnnoMSdat[,dim(peakTableAnnoMSdat)[2]])

    if(isotopeRemove){
    iso1 <- grep("[M+",peakTableAnnoMSdat$isotopes,fixed = TRUE)
    tsv <- tsv[- iso1, ]
    } else {
      iso1 = NULL
      tsv = tsv }

    #if(adductRemove){ }

    }

    cat("\nA 'peakTarget' object with ", dim(tsv)[2] - 11, " samples", "\n")

    #cat("\nNumbers of samples found:", dim(tsv)[2] - 11)
    cat("\no Metabolic peaks found:", dim(tsv)[1])
    cat("\no Filtered isotope peaks ([M + n]):", length(iso1))
    cat("\no Time range: ",round(min(mat$rtmin),2),"-", round(max(mat$rtmax),2), sep = ""," min")
    cat("\no Mass range: ",round(min(mat$mzmin),4),"-", round(max(mat$mzmax),4), sep = ""," m/z\n")


    peakTable <- new("peakTarget")
    peakTable@targetPara <- targetPara
    peakTable@filterPeaks <- iso1
    peakTable@peakTable <- tsv
    peakTable@peakChroma <- MSdat
    return(peakTable)
}

