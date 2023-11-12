#' Generic Box-Cox Analysis Function
#' 
#' The function is implemented as a generic function with methods for classes 
#' "FitAR", "Arima", "ts" and "numeric".
#' 
#' For \eqn{\lambda \ne 0}, the Box-Cox transformation is of x is \eqn{(x^{\lambda}-1)/\lambda)}.
#' If the minimum data value is <= 0, a small positive constant, equal 
#' to the negative of the minimum plus 0.25, is added to all the data values.
#' 
#' @usage BoxCox(object, ...)
#' @param object model object.
#' @param ...	 optional arguments. 
#' @returns No value returned. Graphical output is produced as side-effect.
#'   The plot shows relative likelihood function as well as the MLE and a confidence interval.
#' @note The MASS package has a similar function `boxcox` 
#'   but this is implemented only for regression and analysis of variance.    
#' @author A.I. McLeod and Y. Zhang.
#' @seealso [BoxCox.Arima()], [BoxCox.FitAR()], [BoxCox.ts()], [BoxCox.numeric()].
#' @examples 
#' ## Not run:  # takes a few seconds
#' BoxCox(lynx)
#' out<-FitAR(lynx, c(1,2,4,10,11), ARModel="ARp", MLEQ=FALSE)
#' BoxCox(out)
#' out<-FitAR(lynx, c(1,2,4,5,7,8,10,11,12))
#' BoxCox(out)
#' ## End(Not run)
#' 
#' @export
BoxCox <-
  function(object, ...){
    UseMethod("BoxCox")}
