#' Plot Method for "FitAR" Object
#' 
#' Diagnostic plots: portmanteau p-values; residual autocorrelation plot; normal 
#' probability plot and Jarque-Bera test; spectral density function.
#' 
#' @usage plot(x, terse=TRUE, ...) ## S3 method for class 'FitAR'
#' @param x object of class "FitAR".
#' @param terse if TRUE, only one graph is produced, otherwise many diagnostic plots.
#' @param ... optional arguments.
#' @returns No value is returned. Plots are produced as side-effect.
#' @note When terse=FALSE, numerous graphs are produced. 
#'   Turn on recording to be able to page back and forth between them.
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [summary.FitAR()], [FitAR()], [JarqueBeraTest()], [RacfPlot()], [LBQPlot()].
#' @examples
#' obj<-FitAR(SeriesA, c(1,2,6,7))
#' plot(obj)
#' 
#' @export
plot.FitAR <-
  function(x, terse=TRUE, ...){
    res <- resid(x)
    n <- length(res)
    k <- nrow(x$covHat)
    lag.max <- min(n/4,100) 
    #Ljung-Box plot and residual acf plot
    layout(matrix(c(1,2,1,2),ncol=2))
    LBQPlot(res, k=k, lag.max=lag.max)
    RacfPlot(x)
    if (terse) return(invisible())
    #squared residuals
    #Ljung-Box plot and residual acf plot
    layout(matrix(c(1,2,1,2),ncol=2))
    LBQPlot(res, k=k, lag.max=lag.max, SquaredQ=TRUE)
    RacfPlot(x, SquaredQ=TRUE)
    #normal plot and box plot
    layout(matrix(c(1,2,1,2),ncol=2))
    ans<-JarqueBeraTest(res)
    sti<-paste("Jarque-Bera Test, p-value =",format.pval(ans$pvalue))
    #qqnorm(res, ylab="residuals",main=sti)
    print(qqmath(~res,ylab="Residual Quantiles",xlab="N(0,1) quantiles",main=sti))
    print(bwplot(~res, xlab="residuals"))
    #trace of original and bootstrap
    layout(matrix(c(1,2,1,2),ncol=2))
    y<-x$fits+x$res
    TimeSeriesPlot(y,main="Trace of original time series")
    zBoot<-Boot(x)
    TimeSeriesPlot(zBoot,main="Trace of parametric bootstrap")
    #monotonic spread plot
    if (length(x$phiHat)>0){
      fits<-as.vector(x$fits)
      fits<-fits-mean(fits)
      plot(fits,sqrt(abs(res)), xlab="fit", ylab="sqrt abs residual")
      title(main="Monotone spread plot")
      lines(lowess(fits,sqrt(abs(res)),f=1),lwd=2,col="blue")
    }
    #residual/fit plot
    print(rfs(x))
    #compare theoretical and sample autocorrelations
    #layout(matrix(c(1,2,1,2),ncol=2))
    layout(matrix(c(1,2,3,1,2,3),ncol=2))
    if (length(y)<50)
      LMX<-10
    else
      LMX<-40
    r<-acf(y, lag.max=LMX, type="correlation",plot=FALSE)$acf
    AcfPlot(r[-1], main="Sample ACF")
    g<-TacvfAR(x$phiHat, lag.max=LMX)
    r<-g[-1]/g[1]
    AcfPlot(r, main="Theoretical ACF")
    r<-acf(y, lag.max=LMX, type="correlation",plot=FALSE)$acf
    AcfPlot(r[-1], main="Sample ACF")
    #spectral density plot
    PlotARSdf(x$phiHat,InnovationVariance=var(resid(x)),logSdf=TRUE, main="fitted model sdf")
    invisible()
  }
