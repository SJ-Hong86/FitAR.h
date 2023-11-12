#' Basic ACF Plotting
#' 
#' Produces theoretical correlation plot
#' 
#' @usage AcfPlot(g, LagZeroQ= TRUE, ylab=NULL, main=NULL, ...)
#' @param g         vector of autocorrelations at lags 1,..,length(g).
#' @param LagZeroQ  start plot at lag zero with g\[0\]=1.
#' @param ylab	    vertical axis label.
#' @param main	    plot title.
#' @param ...	      optional graphical parameters.
#' @returns No value. Plot is produced via plot function.
#' @author A.I. McLeod and Y. Zhang.
#' @seealso [acf()].
#' @examples 
#' # Simple example, plot acf for AR(1)
#' phi<-0.8
#' maxLag<-20
#' g<-phi^(1:maxLag)
#' AcfPlot(g)
#' AcfPlot(g, LagZeroQ=FALSE)
#' 
#' Plot the sample inverse partial autocorrelations
#' On the basis of this plot, Cleveland (1972) suggested an ARp(1,2,7)
#' for this data
#' "InverseAcf" <-
#' function(z, p=15){
#'     g<-TacvfMA(GetFitARpLS(z-mean(z),1:p)$phiHat, lag.max=p)
#'     g/g[1]
#' }
#' data(SeriesA)
#' AcfPlot(InverseAcf(SeriesA),LagZeroQ=FALSE)
#' 
#' @export
AcfPlot <-
  function(g, LagZeroQ=TRUE, ylab=NULL, main=NULL, ...){
    if (LagZeroQ) {
      plot( 0:length(g), c(1,g), type="h", ylim=c(-1,1), xlab="lag", ylab=ylab, main=main)
      lines( c(0, length(g)), c(0,0), col="magenta")
    }
    else {
      plot( 1:length(g), g, type="h", ylim=c(-1,1), xlab="lag", ylab=ylab, main=main)
      lines( c(0, length(g)), c(0,0), col="magenta")
    }
  }
