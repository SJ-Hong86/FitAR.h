#' Fast Computation of the Loglikelihood Function in AR
#' 
#' Computation of the loglikelihood is O(1) flops in repeated evaluations of 
#' the loglikelihood holding the data fixed and varying the parameters. 
#' This is useful in exact MLE estimation.
#' 
#' @usage FastLoglikelihoodAR(phi, n, CD)
#' @param phi AR coefficients.
#' @param n length of series.
#' @param CD Champernowne matrix.
#' @details 
#' The details of this computation are described in McLeod and Zhang (2006).
#' @returns Loglikelihood.
#' @author A.I. McLeod and Y. zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [ChampernowneD()], [LoglikelihoodAR()].
#' @examples 
#' # Compute the loglikelihood using the direct method as implemented
#' # in LoglikelihoodAR and using the fast method
#' phi<-PacfToAR(rep(0.5,10))
#' p<-length(phi)
#' z<-SeriesA-mean(SeriesA)
#' n<-length(z)
#' L1<-LoglikelihoodAR(phi, z)
#' cd<-ChampernowneD(z,p,MeanZero=TRUE)
#' L2<-FastLoglikelihoodAR(phi,n,cd)
#' out<-c(L1,L2)
#' names(out)<-c("direct","fast")
#' out
#' 
#' @export
FastLoglikelihoodAR <-
  function(phi, n, CD){
    phis<-c(1,-phi)
    LL <- -log(DetAR(phi))/2-(n/2)*log(sum(crossprod(phis,CD)*phis)/n)
    if (!is.finite(LL)){
      LL<--1e35 }
    LL
  }
