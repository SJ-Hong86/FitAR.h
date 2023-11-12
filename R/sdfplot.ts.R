#' Autoregressive Spectral Density Estimation for "ts" Object
#' 
#' Methods function for "ts".
#' 
#' @usage sdfplot(obj, ...) ## S3 method for class 'ts' 
#' @param obj object, class"ts".
#' @param ... optional arguments.
#' @returns Plot is produced using plot. Matrix with 2 columns containing 
#'   the frequencies and spectral density is returned invisibly.
#' @author A.I. McLeod.
#' @seealso [sdfplot()].
#' @examples
#' data(SeriesA)
#' sdfplot(SeriesA) 
#' 
#' @export
sdfplot.ts <-
  function(obj, ...){
    sdfplot.numeric(as.vector(obj))
  }

