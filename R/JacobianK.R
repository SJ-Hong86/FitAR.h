#' Internal Utility Function
#' 
#' The matrix defined in eqn. (10) of McLeod and Zhang (2006). 
#' Used by the function Jacobian.
#' 
#' @usage JacobianK(zeta, k)
#' @param zeta partial autocorrelations.
#' @param k k-th Jacobian.
#' @returns Matrix.
#' @author A.I. McLeod.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [Jacobian()].
#' @examples
#' JacobianK(rep(0.8,4),3)
#' 
#' @export
JacobianK <-
  function(zeta, k){
    p<-length(zeta)
    pmk<-p-k
    ZeroMatrix<-matrix(rep(0, pmk*k),nrow=pmk,ncol=k)
    IdMatrix<-matrix(rep(c(1,rep(0,k)),k)[1:(k*k)],nrow=k,ncol=k)
    J<-matrix(rep(c(1,rep(0, pmk)),pmk)[1:(pmk^2)],nrow=pmk)
    A<-matrix(rep(0,pmk*k),nrow=pmk,ncol=k)
    phim1j<-c(rev(-PacfToAR(zeta[1:pmk])),1)
    A[,1]<-phim1j[1:pmk]
    for (j in 1:pmk) J[j,pmk+1-j] <- J[j,pmk+1-j]-zeta[pmk+1]
    rbind(cbind(J,A),t(rbind(ZeroMatrix,IdMatrix)))
  }
