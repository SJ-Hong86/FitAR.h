#' Spectral Density of Fitted ARIMA Model
#' 
#' Method for class "Arima" for sdfplot.
#' 
#' @usage sdfplot(obj, ...) ## S3 method for class 'Arima'
#' @param obj object class "Arima".
#' @param ... optional arguments.
#' @returns Plot is produced using plot. Matrix with 2 columns containing 
#'   the frequencies and spectral density is returned invisibly.
#' @author A.I. McLeod.
#' @seealso [ARSdf()], [sdfplot()], [sdfplot.FitAR()], [sdfplot.ts()].
#' @examples
#' sdfplot(SeriesA,c(1,0,1))
#' 
#' @export
sdfplot.Arima <-
  function(obj, ...){
    sti<-paste("ARMA(",obj$arma[1],",",obj$arma[2],")",sep="")
    if (obj$arma[5]>1){
      sti<-paste("S",sti, "(",obj$arma[3],",",obj$arma[4],")",obj$arma[5],sep="")
    }
    PlotARSdf(phi=obj$model$phi, theta=-(obj$model$theta), InnovationVariance=obj$sigma2, logSdf=TRUE, sub=sti, ...)
  }
