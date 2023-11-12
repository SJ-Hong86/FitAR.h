#' Covariance Matrix Residual Autocorrelations for ARz
#' 
#' The ARz subset model is defined by taking a subset of the partial 
#' autocorrelations (zeta parameters) in the AR(p) model. With this function one 
#' can obtain the standard deviations of the residual autocorrelations which can 
#' be used for diagnostic checking with `RacfPlot`.
#' 
#' @usage VarianceRacfARz(zeta, lags, MaxLag, n)
#' @param zeta zeta parameters (partial autocorrelations).
#' @param lags lags in model.
#' @param MaxLag covariance matrix for residual autocorrelations at 
#'   lags 1 ,..., m, where m = MaxLag is computes.
#' @param n length of time series.
#' @details
#' The covariance matrix for the residual autocorrelations is the subset ARz 
#' case is derived in McLeod and Zhang (2006, eqn. 16).
#' @returns The m-by-m covariance matrix of residual autocorrelations at 
#'   lags 1, ..., m, where m = MaxLag.
#' @author A.I. McLeod and Y. Zhang.
#' @references
#' McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation parameterization 
#' subset autoregression. Journal of Time Series Analysis, 27, 599-612. 
#' @seealso [VarianceRacfAR()], [VarianceRacfARp()], [RacfPlot()].
#' @examples
#' # The standard deviations of the first 5 residual autocorrelations 
#' # to a subset AR(1,2,6) model fitted to Series A is
#' v<-VarianceRacfARz(c(0.36,0.23,0.23),c(1,2,6), 5, 197)
#' sqrt(diag(v))
#' 
#' @export
VarianceRacfARz <-
  function(zeta, lags, MaxLag, n){
    PMAX<-max(lags)
    p<-max(lags)
    zta<-numeric(p)
    zta[lags]<-zeta
    phi<-PacfToAR(zta)
    psi<-ARToMA(phi, MaxLag-1)
    PMAX<-length(phi)
    X<-matrix(rep(c(-psi,0),PMAX)[1:(PMAX*MaxLag)],ncol=PMAX)
    X<-X*outer((1:MaxLag),(1:PMAX),">=")
    J<-Jacobian(zta)[,lags]
    X<-X%*%J
    id<-matrix(rep(c(1,rep(0,MaxLag)),MaxLag)[1:(MaxLag^2)],nrow=MaxLag,ncol=MaxLag)
    (id-X%*%solve(InformationMatrixARz(zeta,lags))%*%t(X))/n
  }
