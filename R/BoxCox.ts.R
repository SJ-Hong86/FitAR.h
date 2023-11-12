#' Box-Cox Analysis for a Time Series
#' 
#' The time series is converted to a vector and BoxCox.numeric is used. 
#' 
#' @usage BoxCox(object, interval = c(-1, 1), ...) ## S3 method for class 'ts'
#' @param object a vector of time series values.
#' @param interval interval to be searched.
#' @param ...	 optional arguments.
#' @details 
#' For \eqn{\lambda \ne 0}, the Box-Cox transformation is of x is \eqn{(x^{\lambda}-1)/\lambda)}.
#' If the minimum data value is <= 0, a small positive constant, equal 
#' to the negative of the minimum plus 0.25, is added to all the data values.
#' 
#' It is important not to transform the data when fitting it with AR since the 
#' optimal transformation would be found for the transformed data – not the 
#' original data. Normally this would not be a sensible thing to do.
#' @returns No value returned. Graphical output is produced as side-effect.
#'   The plot shows relative likelihood function as well as the MLE and a confidence interval.
#' @note The MASS package has a similar function `boxcox` 
#'   but this is implemented only for regression and analysis of variance.    
#' @author A.I. McLeod
#' @references Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. 
#'   Journal of Royal Statistical Society, Series B, vol. 26, pp. 211-246.
#' @seealso [BoxCox.FitAR()], [BoxCox.Arima()], [BoxCox.numeric()].
#' @examples 
#' ## Not run:  # takes a few seconds
#' BoxCox(sunspot.year)
#' ## End(Not run)
#' 
#' @export
BoxCox.ts <-
  function(object, interval=c(-1,1), ...){
    BoxCox.numeric(object)
  }
