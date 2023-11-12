#' Parametric Time Series Bootstrap
#' 
#' @description
#' An AR(p) model is fit to the time series using the AIC and then it is simulated.
#' 
#' @usage Boot(obj, R=1, ...) # S3 method for class 'ts'
#' @param obj	 a time series, class "ts".
#' @param R	 number of bootstrap replications.
#' @param ...	 optional arguments.
#' @note Parametric and nonparametric time series bootstraps are discussed
#'   by Davison and Hinkley (1997, Ch.8.2). Nonparametric bootstrap
#'   for time series is available in the function tsboot in the library boot.
#' @returns A time series or vector.
#' @author A.I. McLeod and Y. Zhang.
#' @references Davison, A.C. and Hinkley, D.V. (1997), Bootstrap Methods
#'   and Their Application. Cambridge University Press.
#' @seealso [Boot.FitAR()].
#' @examples 
#' layout(matrix(c(1,2,1,2),ncol=2))
#' TimeSeriesPlot(SeriesA)
#' TimeSeriesPlot(Boot(SeriesA),main="Bootstrap of Series A")
#' 
#' @export
Boot.ts <-
  function(obj, R=1, ...){
    p<-SelectModel(obj, lag.max=ceiling(length(obj)/4), Best=1)
    ans<-FitAR(obj, 1:p)
    Boot.FitAR(ans, R=R)
  }
