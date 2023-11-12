#' Partial Autocorrelations via Durbin-Levinson
#' 
#' Given autocovariances, the partial autocorrelations and/or autoregressive 
#' coefficients in an AR may be determined using the Durbin-Levinson algorithm. 
#' If the autocovariances are sample autocovariances, this is equivalent to using 
#' the Yule-Walker equations. But as noted below our function is more general    
#' than the built-in R functions.
#' 
#' @usage PacfDL(c, LinearPredictor = FALSE)
#' @param c autocovariances at lags 0,1,...,p = length(c)-1.
#' @param LinearPredictor if TRUE, AR coefficients are also determined using 
#'   the Yule-Walker method.
#' @details
#' The Durbin-Levinson algorithm is described in many books on time series and 
#' numerical methods, for example Percival and Walden (1993, eqn 403).
#' @returns
#' If LinearPredictor = FALSE, vector of length p = length(c)-1 containing the 
#' partial autocorrelations at lags 1,...,p. Otherwise a list with components:
#' 
#' * `Pacf` vector of partial autocorrelations.
#' * `ARCoefficients` vector of AR coefficients.
#' * `ResidualVariance` residual variance for AR(p).
#' @note
#' Stationarity is not tested.
#' 
#' Sample partial autocorrelations can also be computed with the acf function 
#' Yule-Walker estimates can be computed with the ar function. Our function 
#' `PacfDL` provides more flexibility since then input c may be any valid 
#' autocovariances not just the usual sample autocovariances. For example, 
#' we can determine the minimum mean square error one-step ahead linear 
#' predictor of order p for theoretical autocovariances from a fractional 
#' arma or other linear process.
#' @author A.I. McLeod and Y. Zhang.
#' @references Percival, D.B. and Walden, A.T. (1993). Spectral Analysis For 
#'   Physical Applications, Cambridge University Press.
#' @seealso [acf()], [ar()].
#' @examples
#' # first define a function to compute the Sample Autocovariances
#' sacvf<-function(z, lag.max){
#' c(acf(z, plot=FALSE, lag.max=lag.max)$acf)*(length(z)-1)/length(z)
#' }
#' # now compute PACF and also fit AR(7) to SeriesA
#' ck<-sacvf(SeriesA, 7)
#' PacfDL(ck)
#' PacfDL(ck, LinearPredictor = TRUE)
#' # compare with built-in functions
#' pacf(SeriesA, lag.max=7, plot=FALSE)
#' ar(SeriesA, lag.max=7, method="yw")
#' # fit an optimal linear predictor of order 10 to MA(1)
#' g<-TacvfMA(0.8,5)
#' PacfDL(g, LinearPredictor=TRUE)
#' 
#' # Compute the theoretical pacf for MA(1) and plot it
#' ck<-c(1,-0.4,rep(0,18))
#' AcfPlot(PacfDL(ck)$Pacf)
#' title(main="Pacf of MA(1), r(1)=-0.4")
#' 
#' @export
PacfDL <-
  function(c, LinearPredictor=FALSE){
    L<-length(c)-1
    d<-c[-1]
    if (L==0) {
      phik<-numeric(0)
      vk<-c[1]
      pi<-numeric(0)
    }
    else {
      phik<-c[2]/c[1]
      pi<-numeric(L)
      pi[1]<-phik
      vk <- c[1]*(1 - pi[1]^2)
    }
    if (L>1){
      for (k in 2:L) {
        vkm1 <- vk
        phikm1 <- phik
        a <- sum(c(1,-phikm1)*rev(d[1:k]))/vk
        phik <- c(phikm1-a*rev(phikm1),a)
        vk <- vkm1*(1-a^2)
        pi[k] <- a
      }
    }
    if (!LinearPredictor)
      list(Pacf=pi)
    else
      list(Pacf=pi, ARCoefficients=phik, ResidualVariance=vk)
  }
