#' Fisher Information Matrix Subset Case, ARz
#' 
#' Computes the large-sample Fisher information matrix per observation for 
#' the AR coefficients in a subset AR when parameterized by the partial autocorrelations.
#' 
#' @usage InformationMatrixARz(zeta, lags)
#' @param zeta vector of coefficients, ie. partial autocorrelations at lags 
#'   specified in the argument `lags`.
#' @param lags lags in subset model, same length as zeta argument.
#' @details
#' The details of the computation are given in McLeod and Zhang (2006, eqn 13). 
#' `FitAR` uses `InformationMatrixARz` to obtain estimates of the standard 
#' errors of the estimated parameters in the subset AR model when partial 
#' autocorrelation parameterization is used.
#' @returns a p-by-p Toeplitz matrix, p = length(zeta).
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()], [InformationMatrixAR()], [InformationMatrixARp()].
#' @examples
#' # Information matrix for ARz(1,4) with parameters 0.9 and 0.9.
#' InformationMatrixARz(c(0.9, 0.9), lags=c(1,4))
#' 
#' @export
InformationMatrixARz <-
  function(zeta,lags){
    pmax<-max(lags)
    z<-numeric(pmax)
    z[lags]<-zeta
    Iar<-InformationMatrixAR(PacfToAR(z))
    if (pmax == 1)
      return(Iar)
    J<-Jacobian(z)
    Iz<-t(J)%*%Iar%*%J
    q<-length(lags)
    Izeta <-matrix(rep(0,q*q),nrow=q, ncol=q)
    for (i in 1:q)
      for (j in 1:q)
        Izeta[i,j]=Iz[lags[i],lags[j]]
    Izeta
  }
