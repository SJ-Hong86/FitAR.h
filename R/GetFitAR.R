#' MLE for AR, ARp and ARz
#' 
#' Obtains the exact MLE for AR(p) or subset AR models ARp or ARz. 
#' This function is used by FitAR. One might prefer to use GetFitAR 
#' for applications such as bootstrapping since it is faster than FitAR.
#' 
#' @usage GetFitAR(z, p, ARModel = "ARz", ...)
#' @param z time series.
#' @param p model order or subset lags.
#' @param ARModel either "ARp" or "ARz" corresponding to `GetFitARp` or `GetFitARz`.
#' @param ... optional arguments which are passed to `GetFitARp` or `GetFitARz`.
#' @details This is just a shell which simply invokes either GetFitARp or GetFitARz.
#' @returns
#' The returned values are as follow:
#' 
#' * `loglikelihood` value of maximized loglikelihood.
#' * `zetaHat` estimated zeta parameters.
#' * `phiHat` estimated phi parameters.
#' * `convergence` result from optim.
#' @author A.I. McLeod.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()].
#' @examples
#' # compare results from GetFitAR and FitAR
#' z<-log(lynx)
#' z<-z - mean(z)
#' GetFitAR(z, c(1,2,8))
#' out<-FitAR(log(lynx), c(1,2,8))
#' out
#' coef(out) 
#' 
#' @export
GetFitAR <-
  function(z, p, ARModel="ARz", ...){
    if (ARModel=="ARp")
      GetFitARpLS(z, p, ...)
    else #also for ARModel="AR"
      GetFitARz(z, p, ...)
  }
