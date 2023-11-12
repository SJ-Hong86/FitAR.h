#' LS for AR(p) and Subset ARp – Short Version
#' 
#' For ARp subset models, the least squares estimates are computed. 
#' The exact loglikelihood is then determined. The estimated parameters 
#' are checked to see if they are in the AR admissible region.
#' 
#' @usage GetFitARpLS(z, pvec)
#' @param z vector or ts object, the time series.
#' @param pvec lags included in subset AR. If pvec = 0, white noise assumed.
#' @details The R function `lsfit` is used.
#' @returns
#' a list with components:
#' 
#' * `loglikelihood` the exact loglikelihood.
#' * `phiHat` estimated AR parameters.
#' * `constantTerm` constant term in the linear regression.
#' * `pvec` lags of estimated AR coefficient.
#' * `res` the least squares regression residuals.
#' * `InvertibleQ` True, if the estimated parameters are in the AR admissible region.
#' * `yX` the y vector and X matrix used for the regression fitting.
#' @note This is a helper function for `FitARp` which is invoked by the main 
#'   package function `FitAR`. Normally the user would `FitAR` since this 
#'   function provides generic print, summary, resid and plot methods but 
#'   `GetFitARpLS` is sometimes useful in iterative computations like 
#'   bootstrapping since it is faster.
#' @author A.I. McLeod.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()], [FitARz()], [GetFitARz()], [FitARp()], [GetFitARpMLE()], [RacfPlot()].
#' @examples
#' # Fit subset AR using LS
#' # normally use FitAR
#' ans<-FitAR(SeriesA, c(1,2,7), ARModel="ARp", MLEQ=FALSE)
#' # could also use FitARp
#' ans<-FitARp(SeriesA, c(1,2,7))
#' # for some applications GetFitARpLS is simpler and faster
#' ansLS<-GetFitARpLS(SeriesA, c(1,2,7))
#' ansLS
#' 
#' @export
GetFitARpLS <-
  function(z, pvec){
    stopifnot(length(z)>0)
    if (length(pvec)==0 || pvec==0){
      phiHat<-numeric(0)
      constantTerm<-mean(z)
      res<-z-constantTerm
      yX<-matrix(z,ncol=1)
      Iq<-TRUE
      LL<-LoglikelihoodAR(0,z, MeanValue=constantTerm)
      covHat<-numeric(0)
    }
    else {
      PMAX<-max(pvec)
      stopifnot(PMAX < length(z))
      Xy <- embed(z, PMAX+1)
      y <- Xy[,1]
      if (length(pvec)==1 && pvec == 1)
        X <- matrix(Xy[,-1], ncol=1)
      else    
        X <- (Xy[,-1])[,pvec]
      yX<-cbind(y,X)
      ans <- lm(y~X)
      res <- resid(ans)
      betaHat <- as.vector(coef(ans)[-1])
      constantTerm <- as.vector(coef(ans)[1])
      phiHat<-numeric(PMAX)
      phiHat[pvec]<-betaHat
      covHat <- vcov(ans)
    }
    Iq<-InvertibleQ(phiHat)
    if (Iq)
      LL<-LoglikelihoodAR(phiHat,z, MeanValue=mean(z))
    else
      LL<- -1e30
    list(loglikelihood=LL, phiHat=phiHat, constantTerm=constantTerm, res=res, pvec=pvec,
         InvertibleQ=Iq,  yX=yX, covHat=covHat)
  }
