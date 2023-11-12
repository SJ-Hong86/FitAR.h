#' Covariance Matrix Residual Autocorrelations for AR
#' 
#' Computes the variance-covariance matrix for the residual autocorrelations in an AR(p).
#' 
#' @usage VarianceRacfAR(phi, MaxLag, n)
#' @param phi vector of AR coefficients.
#' @param MaxLag covariance matrix for residual autocorrelations at 
#'   lags 1 ,..., m, where m = MaxLag is computes.
#' @param n length of time series.
#' @details
#' The covariance matrix for the residual autocorrelations is derived in McLeod 
#' (1978, eqn. 15) for the general ARMA case. With this function one can obtain 
#' the standard deviations of the residual autocorrelations which can be used for 
#' diagnostic checking with `RacfPlot`.
#' @returns The m-by-m covariance matrix of residual autocorrelations at 
#'   lags 1, ..., m, where m = MaxLag.
#' @note
#' The derivation assumes normality of the innovations, mle estimation of the 
#' parameters and a known mean-zero time series. It is easily seen that 
#' the same result still holds for IID innovations with mean zero and finite 
#' variance, any first-order efficient estimates of the parameters including 
#' the AR coefficients and mean.
#' @author A.I. McLeod.
#' @references
#' McLeod, A.I. (1978), On the distribution and applications of residual 
#' autocorrelations in Box-Jenkins modelling, Journal of the Royal 
#' Statistical Society B, 40, 296–302.
#' @seealso [VarianceRacfARp()], [VarianceRacfARz()], [RacfPlot()].
#' @examples
#' VarianceRacfAR(0.5,5,100)
#' 
#' @export
VarianceRacfAR <-
  function(phi, MaxLag, n){
    psi<-ARToMA(phi, MaxLag-1)
    PMAX<-length(phi)
    X<-matrix(rep(c(psi,0),PMAX)[1:(PMAX*MaxLag)],ncol=PMAX)
    X<-X*outer((1:MaxLag),(1:PMAX),">=")
    id<-matrix(rep(c(1,rep(0,MaxLag)),MaxLag)[1:(MaxLag^2)],nrow=MaxLag,ncol=MaxLag)
    (id-X%*%solve(InformationMatrixAR(phi))%*%t(X))/n
  }
