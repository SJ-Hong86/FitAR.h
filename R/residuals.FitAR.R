#' Extract Residuals from "FitAR" Object
#' 
#' Method function.
#' 
#' @usage residuals(object, ...) ## S3 method for class 'FitAR'
#' @param object object of class "FitAR".
#' @param ... optional arguments.
#' @returns Vector of residuals.
#' @author A.I. McLeod.
#' @seealso [FitAR()].
#' @examples
#' out<-FitAR(SeriesA, c(1,2,6,7))
#' resid(out)
#' 
#' @export
residuals.FitAR <-
  function(object, ...) {
    object$res
  }
