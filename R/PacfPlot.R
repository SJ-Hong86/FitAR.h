#' Plot Partial Autocorrelations and Limits
#' 
#' The sample partial autocorrelations and their individual 95 percent confidence 
#' intervals are plotted under the assumption the model is contained in an AR(P), 
#' where P is a specified maximum order.
#' 
#' @usage PacfPlot(z, lag.max = 15, ...)
#' @param z time series.
#' @param lag.max maximum lag, P.
#' @param ... optional parameters passed through to plot.
#' @details The Burg algorithm is used to estimate the PACF.
#' @returns No value is returned. Graphical output is produced as side-effect.
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [ar.burg()], [pacf()].
#' @examples
#' # For the log(lynx) series and taking lag.max=15, the PacfPlot and
#' # the minimum BIC subset selection produce the same result.
#' z<-log(lynx)
#' PacfPlot(z)
#' SelectModel(z,lag.max=15,ARModel="ARz",Best=1,Criterion="BIC")
#' 
#' @export
PacfPlot <-
  function(z, lag.max=15,...){
    ylab<-"pacf"
    lags<-1:lag.max
    zeta<-ARToPacf((ar.burg(z,aic=FALSE,order.max=lag.max))$ar)
    glim=1
    plot(lags, zeta,  type="n", ylim=c(-glim,glim), xlab="lag", ylab=ylab, ...)
    lines(c(0, length(zeta)), c(0,0), col="magenta")
    sdZeta<-sqrt(diag(solve(InformationMatrixARz(zeta,lags)))/length(z))
    segments(lags,zeta-1.96*sdZeta, lags, zeta+1.96*sdZeta, lwd=3, col="blue")
    title(sub="95% confidence intervals for pacf") 
  }
