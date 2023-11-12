#' Theoretical Autocovariance Function of AR
#' 
#' The theoretical autocovariance function of an AR(p) with unit variance is 
#' computed. This algorithm has many applications. In this package it is used 
#' for the computation of the information matrix, in simulating p initial 
#' starting values for AR simulations and in the computation of the exact 
#' mle for the mean. 
#' 
#' @usage TacvfAR(phi, lag.max = 20) 
#' @param phi vector of AR coefficients.
#' @param lag.max computes autocovariances lags 0,1,...,maxlag.
#' @details
#' The algorithm given by McLeod (1975) is used.
#' 
#' The built-in R function ARMAacf could also be used but it is quite complicated 
#' and apart from the source code, the precise algorithm used is not described. 
#' The only reference given for ARMAacf is the Brockwell and Davis (1991) but 
#' this text does not give any detailed exact algorithm for the general case. 
#' 
#' Another advantage of TacvfAR over ARMAacf is that it will be easier for to 
#' translate and implement this algorithm in other computing environments 
#' such as MatLab etc. since the code is entirely written in R.
#' @returns Vector of length = (lag.max+1) containing the autocovariances 
#'   at lags 0,...,lag.max is returned.
#' @author A.I. McLeod.
#' @references
#' McLeod, A.I. (1975), Derivation of the theoretical autocorrelation function 
#' of autoregressive moving-average time series. Applied Statistics, 24, 255-256.
#' @seealso [ARMAacf()], [InformationMatrixAR()], [GetARMeanMLE()], [SimulateGaussianAR()].
#' @examples
#' # calculate and plot the autocorrelations from an AR(2) model
#' # with parameter vector c(1.8,-0.9).
#' g<-TacvfAR(c(1.8,-0.9),20)
#' AcfPlot(g/g[1], LagZeroQ=FALSE)
#' 
#' @export
TacvfAR <-
  function(phi, lag.max = 20)
  {
    p<-length(phi)
    maxlagp1<-lag.max+1
    if(p == 0)
      if(lag.max >= 0) 
        return(c(1, numeric(lag.max)))
    else 
      stop("maxlag invalid")
    r <- p + 1
    b <- numeric(r)
    C <- 1
    phiStar <- numeric(3 * r)
    phiStar[r] <- -1
    phiStar[r + 1:p]<-phi
    a <- matrix(numeric(r^2), ncol = r)
    b[1] <- 1
    for(i in 1:r)
      for(j in 1:r)
        if(j == 1)
          a[i, j] <- phiStar[r + i - 1]  
    else 
      a[i, j] <- phiStar[r + i - j] + phiStar[r + i + j - 2]
    g <- solve(a,  -b)
    if(length(g) < maxlagp1) {
      g <- c(g, numeric(maxlagp1 - r))
      for(i in (r + 1):maxlagp1) 
        g[i] <- phi %*% g[i - 1:p]
    }
    g[1:maxlagp1]
  }
