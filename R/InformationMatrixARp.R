#' Fisher Information Matrix Subset Case, ARp
#' 
#' The large-sample information matrix per observation is computed in a subset 
#' AR with the usual parameterization, that is, a subset of the AR coefficients.
#' 
#' @usage InformationMatrixARp(phi, lags)
#' @param phi vector of coefficients in the subset AR.
#' @param lags vector indicating lags present in phi.
#' @details
#' The subset information matrix is obtained simply by selecting the appropriate 
#' rows and columns from the full information matrix. This function is used by  
#' `FitARp` to obtain the estimated standard errors of the parameter estimates.
#' @returns a p-by-p Toeplitz matrix, p = length(phi).
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [InformationMatrixARp()], [FitARp()], [InformationMatrixARz()].
#' @examples
#' # variances of parameters in a subset ARp(1,2,6)
#' fi<-InformationMatrixARp(c(0.36,0.23,0.23),c(1,2,6))
#' sqrt(diag(solve(fi*197)))
#' 
#' @export
InformationMatrixARp <-
  function(phi, lags){
    p<-length(lags)
    PMAX<-max(lags)
    iar<-toeplitz(TacvfAR(phi, PMAX-1))
    iar[lags,lags]
  }
