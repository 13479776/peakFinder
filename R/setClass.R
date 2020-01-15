
#########class###################

##' An S4 class to represent the parmater for peakTarget
##'
##' @export
##' @slot inputPara The peak picking and anotation parameters.
##' @slot filterPeaks The filtered peak IDs as isotope peaks or others.
##' @slot peakTable Generate a matrix of integrated peak intensity with rows for every group
##' and columns for every sample.
##' @slot peakChroma a XCMSnExp object contains the results of a G/LC-MS data
##' preprocessing that comprises chromatographic peak detection, alignment and
##' correspondence.
##' @return A object of peakTarget
##' @exportClass peakTarget
setClass("peakTarget", slots = list(

  inputPara ="list",
  annoPeaks = "ANY",
  filterPeaks = "ANY",
  peakTable = "data.frame",
  peakChroma="ANY"
))
#setGeneric("peakTarget", function(object) standardGeneric("peakTarget"))
setMethod("show", "peakTarget",
                    function(object) {
                      tsv <- object@peakTable
                      iso <- object@annoPeaks
                      iso1 <- object@filterPeaks
                      mat <- object@peakTable

                      cat("\nA 'peakTarget' object with ", dim(tsv)[2] - 11, " samples", "\n")

                      #cat("\nNumbers of samples found:", dim(tsv)[2] - 11)
                      cat("\n  o Metabolic peaks found:", dim(tsv)[1])
                      #cat("\n  o Annotated isotope peaks ([M + n]):", length(iso))
                      cat("\n  o Confidence peaks: ", length(iso1))
                      cat("\n  o Time range: ",round(min(mat$rtmin),2),"-", round(max(mat$rtmax),2), sep = ""," min")
                      cat("\n  o Mass range: ",round(min(mat$mzmin),4),"-", round(max(mat$mzmax),4), sep = ""," m/z\n")
                      }
          )
