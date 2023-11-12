#' Internal Utility Function: BLUE Mean
#' 
#' This function is not normally used directly by the user. 
#' It is used in the exact mle for mean.
#' 
#' @usage Get1G(phi, n)
#' @param phi a vector of AR coefficients.
#' @param n length of series.
#' @returns A vector used in the mle computation of the mean.
#' @author A.I. McLeod.
#' @seealso [GetARMeanMLE()].
#' @examples
#' # Simulate an AR(2) and compute the exact mle for mean
#' set.seed(7771111)
#' n<-50
#' phi<-c(1.8,-0.9)
#' z<-SimulateGaussianAR(phi, n)
#' g1<-Get1G(phi, length(z))
#' sum(g1*z)/sum(g1)
#' # sample mean
#' mean(z)
#' # more directly with getArMu
#' GetARMeanMLE(z,phi)
#' 
#' @export
Get1G <-
  function(phi,n){
    p<-length(phi)
    x0<-sum(c(1,-phi))^2
    x<-x0-rowSums(GetB(phi))-GetKappa(phi)
    c(x,rep(x0,n-2*p),rev(x))
  }
