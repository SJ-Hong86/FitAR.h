#' Autoregression Simulation
#' 
#' Simulate a mean-zero stationary Gaussian AR(p) time series. 
#' 
#' @usage SimulateGaussianAR(phi, n = 100, InnovationVariance = 1)
#' @param phi vector containing AR coefficients.
#' @param n length of time series.
#' @param InnovationVariance innovation variance.
#' @details
#' The p initial values are simulated using the appropriate multivariate 
#' distribution as was suggested in McLeod (1975). The R function `rnorm()` is used. 
#' @returns A vector of length n, the simulated series.
#' @note 
#' If the process is non-stationary, then random initial values are used  
#' determined by the first p innovations. 
#' @author A.I. McLeod.
#' @references 
#' McLeod, A.I. (1975), Derivation of the theoretical autocorrelation function of 
#' autoregressive moving-average time series, Applied Statistics 24, 255–256.
#'  
#' Percival, D.B. and Walden, A.T. (1993), Spectral Analysis for Physical Applications.
#' @seealso [Boot.FitAR()].
#' @examples
#' # Percival and Walden (1993, p.46) illustrated a time series with a
#' # very peaked spectrum with the AR(4) with coefficients
#' # c(2.7607,-3.8106,2.6535,-0.9238) with NID(0,1) innovations.
#' z<-SimulateGaussianAR(c(2.7607,-3.8106,2.6535,-0.9238),1000)
#' library(lattice)
#' TimeSeriesPlot(z)
#' 
#' @export
SimulateGaussianAR <-
  function(phi, n=100, InnovationVariance=1)
  {
    p<-length(phi)
    a<-rnorm(n, mean=0, sd=sqrt(InnovationVariance))
    if(p==0) return(a)
    if (p==1 && phi==1) return(cumsum(a)) #convenient for unit root test
    z<-numeric(n)
    g<-TacvfAR(phi,p-1)
    if (is.null(g)){ #is null only if non-causal
      warning("Simulating non-stationary stochastic difference equation")
      z[1:p]<-a[1:p]
    } 
    if (p>0)
      z[1:p]<-crossprod(a[1:p],chol(toeplitz(g)))
    for (i in (p+1):n) 
      z[i]=a[i]+sum(rev(phi)*z[(i-p):(i-1)])
    z   
  }
