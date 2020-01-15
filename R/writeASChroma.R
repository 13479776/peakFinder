#' @name writeASChroma
#' @title writeASChroma
#' @description  writeASChroma
#' @param x a object of featureChroma
#' @examples
#' \dontrun{
#' writeASChroma(xx)
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export
writeASChroma <- function(x){
mat <- matrix(1:length(x),dim(x)[1],dim(x)[2],byrow=T)
for(i in 1:dim(x)[1]) {
  cat(i,"...",sep = "")
  for(j in 1:dim(x)[2]){
  name <- paste("F",mat[i,j],".pdf",sep = "")
  pdf(name, 7, 7)
  plot(x[i,j])
  dev.off()
  }
}
}
