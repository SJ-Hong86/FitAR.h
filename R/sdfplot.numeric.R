#' Autoregressive Spectral Density Estimation for "numeric"
#' 
#' Method function for vectors, class "numeric".
#' 
#' @usage sdfplot(obj, ...) ## S3 method for class 'numeric'
#' @param obj object, class"numeric", a vector.
#' @param ... optional arguments.
#' @returns Plot is produced using plot. Matrix with 2 columns containing 
#'   the frequencies and spectral density is returned invisibly.
#' @author A.I. McLeod.
#' @seealso [sdfplot()].
#' @examples
#' sdfplot(lynx)
#' 
#' @export
sdfplot.numeric <-
  function(obj, ...){
    p<-SelectModel(obj, lag.max=ceiling(length(obj)/4), Best=1)
    y<-obj-mean(obj)
    ans<-FitAR(y, 1:p)
    PlotARSdf(ans$phiHat, InnovationVariance=ans$sigsqHat, logSdf=TRUE, ...)
  }
