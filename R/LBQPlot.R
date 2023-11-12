#' Plot Ljung-Box Test P-value vs Lag
#' 
#' The Ljung-Box portmanteau p-value is plotted vs lag.
#' 
#' @usage LBQPlot(res, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
#' @param res residuals.
#' @param lag.max maximum lag.
#' @param StartLag starting lag.
#' @param k number of AR parameters fit.
#' @param SquaredQ default, SquaredQ = FALSE, regular autocorrelations. 
#'   If SquaredQ = TRUE use autocorrelations of squared residuals.
#' @returns Plot is produced as a side-effect. No output.
#' @note This function is normally invoked when plot.FitAR is used.
#' @author A.I. McLeod and Y. Zhang.
#' @references 
#' Ljung, G.M. and Box, G.E.P. (1978) On a measure of lack of fit 
#' in time series models. Biometrika 65, 297-303.
#' 
#' McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation parameterization 
#' for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [plot.FitAR()], [FitAR()].
#' @examples
#' # fit subset AR and plot diagnostic check
#' out<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
#' res<-resid(out)
#' LBQPlot(res)
#' # note that plot produces LBQPlot and RacfPlot
#' plot(out)
#' 
#' @export
LBQPlot <-
  function(res, lag.max=30, StartLag=k+1, k=0, SquaredQ=FALSE){
    stopifnot(k>=0, lag.max-StartLag>0, length(res)>lag.max)
    ans<-LjungBoxTest(res, k=k, lag.max=lag.max, StartLag=StartLag, SquaredQ=FALSE)
    m <- ans[,"m"]
    pv <- ans[,"pvalue"]
    plot(m, pv, xlab="lag", ylab="p-value", ylim=c(0,1), main="Ljung-Box Test")
    abline(h=0.05, col="red", lty=2)
  }
