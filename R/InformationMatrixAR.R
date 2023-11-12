#' Information Matrix for AR(p)
#' 
#' The Fisher large-sample information matrix per observation for the 
#' p coefficients in an AR(p) is computed.
#' 
#' @usage InformationMatrixAR(phi)
#' @param phi vector of length p corresponding to the AR(p) coefficients.
#' @details
#' The Fisher information matrix is computed as the covariance matrix of an AR(p) 
#' process with coefficients given in the argument `phi` and with unit innovation 
#' variance. The `TacvfAR` function is used to compute the necessary autocovariances. 
#' `FitAR` uses `InformationMatrixAR` to obtain estimates of the standard 
#' errors for the estimated parameters in the case of the full AR(p) model.
#' @returns a p-by-p Toeplitz matrix, p = length(phi).
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()], [InformationMatrixARp()], [TacvfAR()], [InformationMatrixARz()].
#' @examples
#' InformationMatrixAR(c(1.8,-0.6))
#' 
#' @export
InformationMatrixAR <-
  function(phi){
    toeplitz(TacvfAR(phi, length(phi)-1))
  }
