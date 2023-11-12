#' Theoretical Autocovariances for Moving Average Process
#' 
#' The theoretical autocovariance function of a MA(q) with unit variance is computed. 
#' 
#' @usage TacvfMA(theta, lag.max = 20)
#' @param theta q parameters in MA(q).
#' @param lag.max number of lags required.
#' @details
#' The first q+1 values are determined using a matrix multiplication - 
#' avoiding a loop. The remaining values set to zero.
#' @returns Vector of length q+1 containing the autocovariances at lags 0,1,...,lag.max.
#' @note See Details in `TacvfAR` for why we prefer to use this algorithm instead of `ARMAacf`.
#' @author A.I. McLeod and Y. Zhang.
#' @references
#' McLeod, A.I. and Zhang, Y. (2006), Partial autocorrelation parameterization  
#' for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [ARMAacf()], [TacvfAR()].
#' @examples
#' TacvfMA(c(1.8,-0.9), 10)
#' 
#' @export
TacvfMA <-
  function(theta, lag.max = 20)
  {
    if(length(theta) == 0)
      if(lag.max >= 0) 
        return(c(1, numeric(lag.max+1)))
    else 
      stop("maxlag invalid")
    maxlagp1 <- lag.max+1
    g<-rep(0, maxlagp1)
    th<-c(-1,theta)
    qth <- length(th)
    x<-c(th, rep(0,qth))
    A <- matrix(0, qth, qth)
    B<-matrix(x[abs(col(A) + row(A)) - 1], qth, qth)
    g1<-c(B%*%th)
    if (length(g1)<maxlagp1)
      g[1:qth]<-g1
    else
      g<-g1[1:maxlagp1]
    g
  }
