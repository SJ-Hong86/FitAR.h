#' Champernowne Matrix
#' 
#' Computes sufficient statistics for AR.
#' 
#' @usage ChampernowneD(z, p, MeanZero = FALSE)
#' @param z time series data.
#' @param p order of the AR.
#' @param MeanZero Assume mean is zero. Default is FALSE so the sample mean is 
#'   subtracted from the data first. Otherwise no sample mean correction is made.
#' @details 
#' This matrix is defined in McLeod & Zhang (2006).
#' @returns The matrix D defined following eqn. (3) of McLeod & Zhang (2006) is computed.
#' @note This function is used by GetFitAR. 
#'   It may be used to compute the exact loglikelihood for an AR.  
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation parameterization 
#'   for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [GetFitARz()], [FastLoglikelihoodAR()], [FitAR()].
#' @examples 
#' # compute the exact concentrated loglikelihood function, (McLeod & Zhang, 2006, eq.(6)),
#' # for AR(p) fitted by Yule-Walker to logged lynx data.
#' #
#' p<-8
#' CD<-ChampernowneD(log(lynx), p)
#' n<-length(lynx)
#' phi<-ar(log(lynx), order.max=p, aic=FALSE, method="yule-walker")$ar
#' LoglYW<-FastLoglikelihoodAR(phi,n,CD)
#' phi<-ar(log(lynx), order.max=p, aic=FALSE, method="burg")$ar
#' LoglBurg<-FastLoglikelihoodAR(phi,n,CD)
#' phi<-ar(log(lynx), order.max=p, aic=FALSE, method="ols")$ar
#' LoglOLS<-FastLoglikelihoodAR(phi,n,CD)
#' phi<-ar(log(lynx), order.max=p, aic=FALSE, method="mle")$ar
#' LoglMLE<-FastLoglikelihoodAR(phi,n,CD)
#' ans<-c(LoglYW,LoglBurg,LoglOLS,LoglMLE)
#' names(ans)<-c("YW","Burg","OLS","MLE")
#' ans
#' # compare the MLE result given by ar with that given by FitAR
#' FitAR(log(lynx),p)
#' 
#' @export
ChampernowneD <-
  function(z, p, MeanZero=FALSE){
    n<-length(z)
    if(MeanZero) y<-z
    else y<-z-mean(z)
    x0<-x<-y
    for (i in 1:p)
      x<-c(x,c(rep(0,i),y[1:(n-i)]))
    x<-matrix(x, nrow=p+1, ncol=length(x0), byrow=TRUE)
    C<-c(x%*%x0)
    A<-toeplitz(C)
    E<-matrix(0, nrow=p+1, ncol=p+1)
    for (j in 1:p)
      for (i in 1:j){
        E[i+1,j+1] <- E[i,j]+y[i]*y[j]+y[n+1-i]*y[n+1-j]
      }
    for (j in 1:(p+1))
      for (i in 1:j)
        E[j,i]=E[i,j]
    A-E
  }
