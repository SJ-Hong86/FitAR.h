#' Autoregressive Spectral Density Estimation for "FitAR"
#' 
#' Methods function for sdfplot. Autoregressive spectral density function 
#' estimation using the result output from FitAR.
#' 
#' @usage sdfplot(obj, ...) ## S3 method for class 'FitAR'
#' @param obj object class "FitAR".
#' @param ... optional arguments.
#' @returns Plot is produced using plot. Matrix with 2 columns containing 
#'   the frequencies and spectral density is returned invisibly.
#' @author A.I. McLeod.
#' @seealso [sdfplot()], [FitAR()].
#' @examples
#' # Use AIC to select best subset model to fit to lynx data and
#' # plot spectral density function
#' pvec<-SelectModel(SeriesA, ARModel="ARp", lag.max=10, Best=1)
#' ans<-FitAR(SeriesA, pvec)
#' sdfplot(ans)
#' 
#' # plot sdf and put your own title
#' z<-c(SeriesA)
#' pvec<-SelectModel(z, ARModel="ARp", lag.max=10, Best=1)
#' ans<-FitAR(z, pvec)
#' sdfplot(ans)
#' title(main="Example SDF")
#' 
#' @export
sdfplot.FitAR <-
  function(obj, ...){
    mti<-obj$DataTitle
    sti<-obj$ModelTitle
    PlotARSdf(obj$phiHat,InnovationVariance=obj$sigsqHat,logSdf=TRUE, main=mti, sub=sti, ...)
  }

