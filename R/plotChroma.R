#' @name plotChroma
#' @title plotChroma
#' @description  plotChroma
#' @param object a object of featureChroma
#' @param num peaks ID
#' @param expandRT RT
#' @param facet_wrap facet_wrap
#' @param leg.pos xx
#' @param SavGolay Smooth data with a Savitzky-Golay smoothing filter from signal 0.7-6 package.
#' @param SavGolay.W critical frequencies of the Savitzky-Golay smoothing filter filter. W must be a scalar for low-pass and high-pass filters, and W must be a two-element vector c(low, high) specifying the lower and upper bands. For digital filters, W must be between 0 and 1 where 1 is the Nyquist frequency.
#' @examples
#' \dontrun{
#' plotChroma(xx)
#' }
#' @keywords peaks, peakTarget, readPKparam
#' @export

plotChroma <- function(object, num, expandRT = 15,facet_wrap = TRUE,leg.pos = "none",SavGolay = FALSE, SavGolay.W = 1/3) {

  dat <- chromaList(object, num)
  if(facet_wrap)  {
    minrt <- min(dat$rtime) + expandRT
    maxrt <- max(dat$rtime) - expandRT
    if(maxrt- minrt <= 0){
      stop("expandRT was too large!")
    }

    if(SavGolay) {
    sgdat <- dat[!is.na(dat$intensity),]
    #sgdat$intensity <- sgolayfilt(sgdat$intensity)
    sgdat$intensity <- signal::filtfilt(signal::butter(5,SavGolay.W),sgdat$intensity)
    dat <- sgdat
    }

    ggplot2::ggplot(dat,ggplot2::aes(y = intensity)) +
      ggplot2::geom_line(data=dat[!is.na(dat$intensity),],ggplot2::aes(x = rtime, color = sampleID), size = 0.3,na.rm =FALSE,linejoin = "mitre",show.legend = FALSE)  +
      ggplot2::geom_area(mapping = ggplot2::aes(x = as.numeric(ifelse(rtime>=minrt & rtime< maxrt , rtime, NA))),fill = "grey",alpha= 0.3) +
      ggplot2::geom_vline(xintercept= minrt, linetype="dashed", size=0.5, color="darkgrey",alpha = 0.7) +
      ggplot2::geom_vline(xintercept= maxrt, linetype="dashed", size=0.5, color="darkgrey",alpha = 0.7) +
      ggplot2::facet_wrap(~ sampleID) +
      ggplot2::labs(
        #title=paste("Peak m/z",round(mz(pk)[1],3),"-",round(mz(pk)[2],3)),
        #subtitle = "Perc Returns for Personal Savings",
        y="Intensity",
        x = "Retention Time (s)",
        caption=dat$name[1]) +
      ggplot2::theme_bw()
  } else {
    minrt <- min(dat$rtime) + expandRT
    maxrt <- max(dat$rtime) - expandRT

    if(SavGolay) {
      sgdat <- dat[!is.na(dat$intensity),]
      #sgdat$intensity <- sgolayfilt(sgdat$intensity)
      sgdat$intensity <- signal::filtfilt(signal::butter(5,SavGolay.W),sgdat$intensity)
      dat <- sgdat
    }

    ggplot2::ggplot(dat, ggplot2::aes(rtime, intensity)) +
      ggplot2::geom_line(data=dat[!is.na(dat$intensity),],ggplot2::aes(color = sampleID), size = 0.6,na.rm =FALSE,linejoin = "mitre")  +
      #geom_area(mapping = aes(x = ifelse(rtime>=minrt & rtime< maxrt , rtime, NA)), fill = "grey",alpha= 0.3) +
      #geom_smooth(data=dat[!is.na(dat$intensity),],aes(color = sampleID), size = 0.25,
      #            na.rm =FALSE,linejoin = "mitre",method = "glm", se = FALSE,
      #            formula=y~ns(x,15),
      #            family= gaussian(link="log"))  +

      ggplot2::geom_vline(xintercept= minrt , linetype="dashed", size=1, color="grey70",alpha = 0.7) +
      ggplot2::geom_vline(xintercept= maxrt, linetype="dashed", size=1, color="grey70",alpha = 0.7) +
      ggplot2::labs(
        #title=paste("Peak m/z",round(mz(pk)[1],3),"-",round(mz(pk)[2],3)),
        #subtitle = "Perc Returns for Personal Savings",
        y="Intensity",
        x = "Retention Time (s)",
        caption=dat$name[1]) +
      #theme_bw() +
      ggplot2::theme(
            legend.position= leg.pos,
            panel.background = element_rect(fill = "white", colour = "grey70"),
            #panel.grid=element_line(color="grey50",size=1),
            #panel.grid.major=element_line(size=1,linetype =0.5,color="grey70"),
            #panel.grid.minor=element_line(size=0.5,linetype ="dashed",color="grey70"),
            panel.border = element_rect(size = 1,fill = NA, color = "black"),
            #panel.grid.major = element_line(size = 0.1, color = "grey70")
            #panel.grid.minor = element_line(size = 0.5, color = "grey70")
            #coord_fixed(ratio=1.5)
            plot.margin = margin(2, 1, 1.5, 1, "cm")
            )
  }
}
