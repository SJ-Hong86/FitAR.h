#' Autoregressive Spectral Density Estimation for "ar"
#' 
#' Method for class "ar" for sdfplot.
#' 
#' @usage sdfplot(obj, ...) ## S3 method for class 'ar'
#' @param obj class "ar" object, output from `ar`.
#' @param ... optional arguments.
#' @returns Plot is produced using plot. Matrix with 2 columns containing 
#'   the frequencies and spectral density is returned invisibly.
#' @author A.I. McLeod.
#' @seealso [ARSdf()], [sdfplot()], [sdfplot.FitAR()], [sdfplot.ts()].
#' @examples
#' # Fit AR(p) using AIC model selection and Burg estimates.
#' # Plot spectral density estimate.
#' ans<-ar(lynx, lag.max=20)
#' sdfplot(ans)
#' 
#' @export
sdfplot.ar <-
  function(obj, ...){
    sti<-paste("AR(",obj$order,")",sep="")
    PlotARSdf(obj$ar, InnovationVariance=obj$var.pred, logSdf=TRUE, sub=sti, ...)
  }
