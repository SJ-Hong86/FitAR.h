#' Ljung-Box Test for Randomness
#' 
#' The Ljung-Box Portmanteau test for the goodness of fit of ARIMA models 
#' is implemented.
#' 
#' @usage LjungBoxTest(res, k=0, lag.max=30, StartLag=1, SquaredQ=FALSE)
#' @param res residuals.
#' @param k number of ARMA parameters, default k = 0.
#' @param lag.max maximum lag, default MaxLag = 30.
#' @param StartLag test is done for lags m=StartLag:MaxLag, default StartLag = 1.
#' @param SquaredQ if TRUE, use squared residuals for ARCH test, default Squared = FALSE.
#' @details
#' This test is described in detail in Wei (2006, p.153, eqn. 7.5.1). 
#' The df are given by h-k, where h is the lag, running from StartLag to lag.max, 
#' when h-k < 1, it is reset to 1. This is ok, since the test is conservative in this case.
#' 
#' A powerful test for ARCH and other nonlinearities is obtained by using 
#' squared values of the series to be tested (McLeod & Li, 1983). Note that if 
#' Squared=TRUE is used the data "res" is centered by sample mean correction 
#' before squaring.
#' @returns A matrix with columns labelled m, Qm, pvalue, where m is the lag 
#'   and Qm is the Ljung-Box Portmanteau statistic and pvalue its p-value.
#' @note This test may also be used to test a time series for randomness taking k = 0.   
#' @author A.I. McLeod.
#' @references
#' W.W.S. Wei (2006, 2nd Ed.), Time Series Analysis: Univariate and Multivariate Methods.
#' 
#' A.I. McLeod. & W.K. Li (1983), Diagnostic checking ARMA time series models 
#' using squared-residual autocorrelations, Journal of Time Series Analysis 4, 269–273.
#' @seealso [Box.test()].
#' @examples
#' # test goodness-of-fit of AR(2) model fit to log(lynx)
#' data(lynx)
#' z<-log(lynx)
#' ans<-FitAR(z, 1:2)
#' # notice that the test is also available as a component of the output of FitAR 
#' ans$LjungBox
#' # a plot of the test is produced
#' plot(ans)
#' # doing the test manually
#' res<-resid(ans)
#' LjungBoxTest(res, k=2, lag.max=20, StartLag=5)
#' 
#' # test for subset case
#' z<-log(lynx)
#' pvec<-SelectModel(z, ARModel="ARz", Criterion="BIC", lag.max=10, Best=1)
#' ans<-FitAR(z, pvec)
#' plot(ans)
#' res<-resid(ans)
#' LjungBoxTest(res, k=length(pvec), lag.max=20, StartLag=11)
#' # test for ARCH effect,
#' LjungBoxTest(res,SquaredQ=TRUE)
#' 
#' @export
LjungBoxTest <-
  function(res, k=0, lag.max=30, StartLag=1, SquaredQ=FALSE){
    stopifnot(k>=0, StartLag>=1, lag.max>=StartLag)
    n<-length(res)
    L0<-StartLag
    if (SquaredQ) {
      z<-(res-mean(res))^2
      kpar<-0
    }
    else {
      z<-res
      kpar<-k
    }
    ra<-(acf(z, lag.max=lag.max, plot=FALSE)$acf)[-1]
    lags<-L0:lag.max
    QQ<-n*(n+2)*cumsum((ra^2)/(n-(1:lag.max)))[lags]
    df <- ifelse(lags-kpar > 0, lags-kpar, 1)
    pv<-1-pchisq(QQ,df)
    QQ<-round(QQ,2)
    a<-matrix(c(lags,QQ,pv),ncol=3)
    dimnames(a)<-list(rep("",length(QQ)),c("m","Qm", "pvalue"))
    a
  }
