#' Exact Loglikelihood for AR
#' 
#' The exact loglikelihood function, defined in eqn. (6) of McLeod & 
#' Zhang (2006) is computed. Requires O(n) flops, n = length(z).
#' 
#' @usage LoglikelihoodAR(phi, z, MeanValue = 0)
#' @param phi AR parameters.
#' @param z time series data, not assumed mean corrected.
#' @param MeanValue usually this is mean(z) but it could be another value for 
#'   example the MLE of the mean.
#' @details
#' Eqn (6) of McLeod and Zhang (2006) may be written 
#' \deqn{-(n/2)\mathrm{log}(\hat{\sigma}_a^2)-(1/2)\mathrm{log}(g_p)} 
#' where \eqn{\hat{\sigma}_a^2} is the residual variance and \eqn{g_p} is the
#' covariance determinant. 
#' @returns The value of the loglikelihood is returned.
#' @note
#' No check is done for stationary-causal process.
#' 
#' For MLE computation it is better to use FastLoglikelihoodAR since for repeated 
#' likelihood evaluations this requires only O(1) flops vs O(n) flops, where n = 
#' length(z).
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FastLoglikelihoodAR()].
#' @examples
#' # Fit a subset model to Series A and verify the loglikelihood
#' out<-FitAR(SeriesA, c(1,2,7))
#' out
#' # either using print.default(out) to see the components in out
#' # or applying LoglikelihoodAR () by first obtaining the phi parameters as out$phiHat.
#' LoglikelihoodAR(out$phiHat, SeriesA, MeanValue=mean(SeriesA))
#' 
#' @export
LoglikelihoodAR <-
  function(phi, z, MeanValue=0){
    if(length(phi)==0) phi=0
    phis<-c(1,-phi)
    y<-z-MeanValue
    n<-length(z)
    -log(DetAR(phi))/2 - (n/2)*log(sum(crossprod(phis,ChampernowneD(y,length(phis)-1,MeanZero=TRUE))*phis)/n)
  }

