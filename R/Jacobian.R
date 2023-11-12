#' Jacobian AR-coefficients to Partial Autocorrelations
#' 
#' This is more or less and internal routine used by InformationMatrixZeta 
#' but it is described in more details since it may be useful in other computations.
#' 
#' @usage Jacobian(zeta)
#' @param zeta partial autocorrelation parameters.
#' @details The computation is described in detail in McLeod and Zhang (2006, Section 2.2).
#' @returns square matrix of order length(zeta).
#' @author A.I. McLeod.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [InformationMatrixARz()].
#' @examples
#' #In McLeod and Zhang (2006, p.603) a symbolic example is given for the AR(4).
#' Jacobian(rep(0.8,4))
#' 
#' @export
Jacobian <-
  function(zeta){
    if(length(zeta)==1) 
      return(1)
    J<-JacobianK(zeta,1)
    if (length(zeta)==2) 
      return(J)
    for (j in 2:(length(zeta)-1)){
      J<-J%*%JacobianK(zeta,j)
    }
    J
  }
