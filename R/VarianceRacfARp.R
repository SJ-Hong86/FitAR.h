#' Covariance Matrix Residual Autocorrelations for ARp
#' 
#' The ARp subset model is defined by taking a subset of the parameters in the 
#' regular AR(p) model. With this function one can obtain the standard deviations 
#' of the residual autocorrelations which can be used for diagnostic checking 
#' with `RacfPlot`.
#' 
#' @usage VarianceRacfARp(phi, lags, MaxLag, n)
#' @param phi vector of AR coefficients.
#' @param lags lags in subset AR.
#' @param MaxLag covariance matrix for residual autocorrelations at 
#'   lags 1 ,..., m, where m = MaxLag is computes.
#' @param n length of time series.
#' @details
#' The covariance matrix for the residual autocorrelations is derived in McLeod 
#' (1978, eqn. 15) for the general ARMA case. McLeod (1978, eqn. 35) specializes 
#' this result to the subset case.
#' @returns The m-by-m covariance matrix of residual autocorrelations at 
#'   lags 1, ..., m, where m = MaxLag.
#' @author A.I. McLeod.
#' @references
#' McLeod, A.I. (1978), On the distribution and applications of residual 
#' autocorrelations in Box-Jenkins modelling, Journal of the Royal 
#' Statistical Society B, 40, 296–302.
#' @seealso [VarianceRacfAR()], [VarianceRacfARz()], [RacfPlot()].
#' @examples
#' # The standard deviations of the first 5 residual autocorrelations 
#' # to a subset AR(1,2,6) model fitted to Series A is
#' v<-VarianceRacfARp(c(0.36,0.23,0.23),c(1,2,6), 5, 197)
#' sqrt(diag(v))
#' 
#' @export
VarianceRacfARp <-
  function(phi,lags, MaxLag, n){
    psi<-ARToMA(phi, MaxLag-1)
    PMAX<-max(lags)
    X<-matrix(rep(c(psi,0),PMAX)[1:(PMAX*MaxLag)],ncol=PMAX)
    X<-X*outer((1:MaxLag),(1:PMAX),">=")
    X<-X[,lags]
    id<-matrix(rep(c(1,rep(0,MaxLag)),MaxLag)[1:(MaxLag^2)],nrow=MaxLag,ncol=MaxLag)
    (id-X%*%solve(InformationMatrixARp(phi, lags))%*%t(X))/n
  }
