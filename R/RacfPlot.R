#' Residual Autocorrelation Plot
#' 
#' Residual autocorrelation plot for "FitAR" objects. This plot is useful for 
#' diagnostic checking models fit with the function `FitAR`.
#' 
#' @usage RacfPlot(obj, lag.max = 1000, SquaredQ=FALSE, ylab="")
#' @param obj output from FitAR.
#' @param lag.max maximum lag. Set to 1000 since minimum of this value and the value in the obj is used.
#' @param SquaredQ default is FALSE. For squared residual autocorrelations, set to TRUE.
#' @param ylab y-axis label.
#' @details
#' The standard deviations of the residual autocorrelations are obtained from 
#' McLeod (1978, eqn.16) or McLeod and Zhang (2006, eqn.16). Simultaneous 
#' confidence bounds are shown and constructed using the Bonferonni 
#' approximation as suggested by Hosking and Ravishanker (1993).
#' @returns Plot is produced as a side-effect. No output.
#' @note This function is normally invoked when plot.FitAR is used.
#' @author A.I. McLeod and Y. Zhang.
#' @references
#' Hosking, J.R.M. and Ravishanker, N. (1993) Approximate simultaneous 
#' significance intervals for residual autocorrelations of autoregressive-moving 
#' average time series models. Journal of Time Series Analysis 14, 19-26.
#' 
#' McLeod, A.I. (1978), On the distribution and applications of residual 
#' autocorrelations in Box-Jenkins modelling, Journal of the Royal Statistical 
#' Society B 40, 296-302.
#' 
#' McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation parameterization 
#' for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [plot.FitAR()], [FitAR()].
#' @examples
#' # fit subset AR and plot diagnostic check
#' data(SeriesA)
#' out<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
#' RacfPlot(out)
#' # note that plot produces LBQPlot and RacfPlot
#' plot(out)
#' # check squared residuals
#' RacfPlot(out, SquaredQ=TRUE)
#' 
#' @export
RacfPlot <-
  function(obj, lag.max=1000, SquaredQ=FALSE, ylab=" "){
    modN<-paste("Residual autocorrelation:",obj$ModelTitle,"\nSimultaneous 95% Interval")
    tb<-obj$RacfMatrix
    ra<-tb[,1]
    sdra<-tb[,2]
    MXL<-min(length(ra),lag.max)
    if (SquaredQ) {
      res<-residuals(obj)
      ra<-(acf(res^2, MXL, type="correlation",plot=FALSE)$acf)[-1]
      sdra<-rep(1/sqrt(length(res)),MXL)
    }
    else {
      sigma<-sqrt(obj$sigsqHat)
      ra<-ra[1:MXL]
      sdra<-sdra[1:MXL]
    }
    clim<-qnorm(1-0.05/(2*MXL))
    glim<-max(abs(ra),clim*max(sdra))*1.05
    lags<-1:MXL
    plot(lags, ra, type="h", ylim=c(-glim,glim), xlab="lag", ylab=ylab)
    lines(c(0, length(ra)), c(0,0), col="magenta")
    lines(lags, clim*sdra, col="blue")
    lines(lags, -clim*sdra, col="blue")
    title(main=modN) 
  }
