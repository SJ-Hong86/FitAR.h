#' Innovation Residuals in AR
#' 
#' @description
#' Obtains the residuals (estimated innovations).
#' The residuals for t=1,...,p are obtained using
#' the backforecasting algorithm of Box and Jenkins (1970).
#' 
#' @usage BackcastResidualsAR(y, phi, Q = 100, demean=TRUE)
#' @param y	 a time series or vector.
#' @param phi	 AR coefficients, lags 1,...,p.
#' @param Q	 for backcasting, the AR is approximated by an MA(Q).
#' @param demean	 subtract sample mean.
#' @details 
#' The backforecasting algorithm is described in detail in the book of Box and
#' Jenkins (1970). The idea is to compute the expected value of the innovation
#' assuming a high-order MA(q).
#' @returns Vector of residuals.
#' @note No check is done that the AR is causal-stationary.
#' @author A.I. McLeod and Y. Zhang.
#' @references Box and Jenkins (1970). Time Series Analysis: Forecasting and Control.
#' @seealso [InvertibleQ()], [FitAR()].
#' @examples 
#' # compare residuals obtained using backcasting with fitted parameters and
#' # the residuals extracted from output of FitAR. They are identical.
#' p<-11
#' out<-FitAR(log(lynx), p)
#' phi<-out$phiHat # fitted parameters
#' resphi<-BackcastResidualsAR(log(lynx), phi)
#' sum(abs(resphi-resid(out)))
#' 
#' @export
BackcastResidualsAR <-
  function(y,phi,Q=100,demean=TRUE){
    if (demean)
      z <- y-mean(y)
    else
      z <- y
    p <- length(phi)
    if (p==0)
      a<-z
    else {
      n<-length(z)
      p<-length(phi)
      nQ<-n+Q
      a<-e<-zF<-zR<-numeric(nQ)
      r<-p+1
      zR[1:n]<-rev(z)
      zF[Q+(1:n)]<-z
      for (i in r:n)
        e[i]<-zR[i]-crossprod(phi,zR[i-(1:p)])
      for (i in 1:Q)
        zR[n+i]<-crossprod(phi,zR[(n+i)-(1:p)])
      zF[1:Q]<-rev(zR[n+(1:Q)])
      for (i in r:nQ)
        a[i]<-zF[i]-crossprod(phi,zF[i-(1:p)])
      zF[Q+(1:n)]<-z
      a<-a[Q+(1:n)]
    }
    a
  }
