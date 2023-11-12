#' Internal Utility Function
#' 
#' The user would not normally use this function. The function is needed for 
#' exact mle for mean. Used in Get1G which is called from GetARMeanMLE.
#' 
#' @usage GetB(phi)
#' @param phi a vector of AR coefficients.
#' 
#' @export
GetB <-
  function(phi){
    p<-length(phi)
    if (p == 1)
      a<-phi^2
    else {
      a<-p*as.vector(acf(phi,lag.max=p,type="covariance",demean=FALSE,plot=FALSE)$acf)
      for (i in 1:(p-1))
        a<-c(a,p*as.vector(acf(c(phi[-(1:i)],rep(0,i)),lag.max=p,type="covariance",demean=FALSE,plot=FALSE)$acf)[1:(p-i)])
    }
    FromSymmetricStorageUpper(a)
  }
