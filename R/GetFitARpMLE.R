#' Exact MLE for subset ARp Models
#' 
#' Uses built-in function `arima` to fit subset ARp model, that is, 
#' the subset model is formed by constraining some coefficients to zero.
#' 
#' @usage GetFitARpMLE(z, pvec)
#' @param z time series.
#' @param pvec lags included in AR model. If pvec = 0, white noise model assumed.
#' @details Due to the optimization algorithms used by `arima`, this method 
#'   is not very reliable. The optimization may simply fail. Example 1 shows 
#'   it working but in Example 2 below it fails.
#' @returns
#' a list with components:
#' 
#' * `loglikelihood` the exact loglikelihood.
#' * `phiHat` estimated AR parameters.
#' * `constantTerm` constant term in the linear regression.
#' * `pvec` lags of estimated AR coefficient.
#' * `res` the least squares regression residuals.
#' * `InvertibleQ` True, if the estimated parameters are in the AR admissible region.
#' @author A.I. McLeod.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()], [FitARz()], [GetFitARz()], [FitARp()], [RacfPlot()].
#' @examples
#' # Example 1. MLE works
#' z<-log(lynx)
#' p<-c(1,2,4,7,10,11)
#' GetFitARpMLE(z, p)
#' #
#' # Example 2. MLE fails with error.
#' p<-c(1,2,9,12)
#' ## Not run: GetFitARpMLE(z, p)
#' 
#' @export
GetFitARpMLE <-
  function(z, pvec){
    stopifnot(length(z)>0)
    if ((length(pvec)==1 && pvec==0) || length(pvec)==0){
      phiHat<-numeric(0)
      Iq<-TRUE
      constantTerm<-mean(z)
      LL<-LoglikelihoodAR(phiHat,z, MeanValue=constantTerm)
      res<-z
    }
    else {
      P <- max(pvec)
      stopifnot(length(z)>P)
      ind <- rep(0, P+1)
      ind[pvec] <- NA
      ind[P+1] <- NA
      out<-arima(z, order=c(P, 0, 0), fixed=ind,transform.pars=FALSE)
      res<-resid(out)
      estimates<-coef(out)
      constantTerm<-as.vector(estimates[P+1])
      if (P>0) {
        phiHat<-estimates[1:P]
        Iq<-InvertibleQ(phiHat)
        if (Iq)
          LL<-LoglikelihoodAR(phiHat,z, MeanValue=constantTerm)
        else
          LL<- -1e30
      }
    }
    list(loglikelihood=LL, phiHat=phiHat, constantTerm=constantTerm,res=res,pvec=pvec,InvertibleQ=Iq)
  }
