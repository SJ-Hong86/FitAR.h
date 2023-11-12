#' Fitted Values from "FitAR" Object
#' 
#' Method function, extracts fitted values from FitAR object.
#' 
#' @usage fitted(object, ...) ## S3 method for class 'FitAR'
#' @param object object of class "FitAR".
#' @param ... optional arguments.
#' @returns Vector of fitted values.
#' @author A.I. McLeod and Y. zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()].
#' @examples 
#' out<-FitAR(SeriesA, c(1,2,6,7))
#' fitted(out)
#' 
#' @export
fitted.FitAR <-
  function(object, ...){
    object$fits
  }
