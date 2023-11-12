#' Plot AR or ARMA Spectral Density
#' 
#' Constructs a plot of the AR spectral density function.
#' 
#' @usage PlotARSdf(phi = NULL, theta = NULL, units = "radial", logSdf = FALSE, InnovationVariance = 1, main = NULL, sub = NULL, lwd=3, col="blue", plotQ=TRUE, ...)
#' @param phi AR Coefficients.
#' @param theta MA Coefficients.
#' @param units default is "radial".
#' @param logSdf default is FALSE otherwise log sdf is plotted.
#' @param InnovationVariance innovation variance, default is 1.
#' @param main optional plot title.
#' @param sub optional subtitle.
#' @param lwd optional lwd for plot, default lwd=3.
#' @param col optional col for plot. Default "blue".
#' @param plotQ True, plot otherwise not.
#' @param ... optional arguments.
#' @details The spectral density function is symmetric and defined in (-pi, pi) 
#'   but plotted over (0, pi). If units are not "radial", it is plotted over (0, 0.5).
#' @returns Plot is produced using plot. Matrix with 2 columns containing the 
#'   frequencies and spectral density is returned invisibly.
#' @author A.I. McLeod and Y. Zhang.
#' @seealso [ARSdf()].
#' @examples
#' # AR(1)
#' PlotARSdf(0.8)
#' # MA(1)
#' PlotARSdf(theta=0.8)
#' # ARMA(1,1)
#' PlotARSdf(0.9,0.5)
#' # white noise
#' PlotARSdf()
#' 
#' @export
PlotARSdf <-
  function(phi=NULL, theta=NULL, units="radial",logSdf=FALSE,InnovationVariance=1, main=NULL, sub=NULL, lwd=3, col="blue", plotQ=TRUE, ...){
    sdf<-InnovationVariance
    if (!is.null(phi))  
      sdf<-ARSdf(phi)*sdf
    if (!is.null(theta))
      sdf<-sdf/ARSdf(theta)
    if (is.null(phi)&&is.null(theta))
      sdf<-sdf*ARSdf(0)
    f0<-(1/length(sdf))*(1:length(sdf))*0.5
    if (units=="radial")
      f <- f0*pi
    else
      f<-f0
    yl<-"Spectral Density"
    if (logSdf){
      sdf<-log(sdf)
      yl<-"Log Sdf"
    }
    if (plotQ)    
      plot(f,sdf,type="l",xlab="frequency",ylab=yl, main=main, sub=sub, lwd=lwd, col=col, ...)
    invisible(matrix(c(f,sdf),ncol=2))
  }
