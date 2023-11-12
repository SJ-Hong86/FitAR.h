#' Exact MLE for AR(p) and Subset ARz – Short Version
#' 
#' Obtain the exact MLE for AR(p) or subset ARz model. This function is used 
#' by `FitAR` and `FitARz`. One might prefer to use GetFitARz for applications 
#' such as bootstrapping since it is faster.
#' 
#' @usage GetFitARz(z, pvec, MeanValue=0, ...)
#' @param z time series.
#' @param pvec lags included in AR model. If pvec = 0, white noise model assumed.
#' @param MeanValue by default it is assumed the mean of z is 0.
#' @param ... optional arguments passed through to optim.
#' @details The built-in function optim is used to obtain the MLE estimates 
#'   for an AR or subset AR. First "BFGS" is tried. This usually works fine. 
#'   In the rare cases where convergence is not obtained, "Nelder-Mead" is used. 
#'   A warning message is given if this happens.
#' @returns
#' a list with components:
#' 
#' * `loglikelihood` value of maximized loglikelihood.
#' * `zetaHat` estimated zeta parameters.
#' * `phiHat` estimated phi parameters.
#' * `convergence` result from optim.
#' * `pvec` lags of estimated AR coefficient.
#' * `algorithm` "BFGS" or "Nelder-Mead".
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()], [FitARz()], [GetFitARpMLE()], [RacfPlot()].
#' @examples
#' # compare results from GetFitARz and FitAR
#' z<-log(lynx)
#' z<-z - mean(z)
#' GetFitARz(z, c(1,2,8))
#' out<-FitAR(log(lynx), c(1,2,8), ARModel="ARz")
#' out
#' coef(out)
#' 
#' @export
GetFitARz <-
  function(z, pvec, MeanValue=0, ...){
    stopifnot(length(z)>0, length(z)>=2*max(pvec), length(pvec)>0,all(pvec>=0))
    pVec<-pvec[pvec>0]
    y <- z-MeanValue
    n<-length(y)
    if (length(pVec)==0 || length(pvec)==0) 
      return(list(loglikelihood=-(n/2)*log(sum(y^2)/n),zetaHat=NULL,phiHat=NULL,convergence=0,algorithm="cubic"))
    PMAX <- max(pVec)
    if (PMAX == 1){
      phiHat <- AR1Est(y)
      LogL <- LoglikelihoodAR(phiHat,z)
      return(list(loglikelihood=LogL, zetaHat=phiHat, phiHat=phiHat, convergence=0))
    }
    
    PEFF<-length(pVec)
    CD<-ChampernowneD(y,PMAX,MeanZero=TRUE)
    xpar<-numeric(PMAX)
    EntropyAR<-function(x){
      if (max(abs(x))>0.999)
        out<-1e35
      else {
        xpar[pVec]<-x
        out<- -FastLoglikelihoodAR(PacfToAR(xpar),n,CD)
      }
      out
    }
    xinit<-ARToPacf(ar.burg(y, aic=FALSE, order.max=PMAX, demean=FALSE)$ar)[pVec]
    #Sometimes there are problems with "L-BFGS-B" -- it frequently tests the endpoints which is bad news due
    #to numerical problems such as ARToPacf(PacfToAR(rep(0.99,20))) is not correct!
    #So it is better to use "BFGS" with a penalty instead.
    #ans<-optim(xinit,EntropyAR,method="L-BFGS-B", lower=rep(-0.9999,PEFF), upper=rep(0.9999,PEFF),control=list(trace=6),...)
    ans<-optim(xinit,EntropyAR,method="BFGS", control="trace", ...)
    alg<-1
    if(ans$convergence != 0) {
      alg<-2
      warning("Convergence problem. convergence=", ans$convergence)
      warning("Trying Nelder-Mead algorithm ...")
      ans<-optim(xinit,EntropyAR,method="Nelder-Mead", ...)
      if(ans$convergence != 0)
        warning("Still convergence problem, convergence= ", ans$convergence)
    }
    zetaHat<-ans$par
    zetas<-numeric(PMAX)
    zetas[pVec]<-zetaHat
    list(loglikelihood=-ans$value, zetaHat=ans$par, phiHat=PacfToAR(zetas),convergence=ans$convergence, algorithm=c("BFGS","Nelder-Mead")[alg],pvec=pvec)
  }
