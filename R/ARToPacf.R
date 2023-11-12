#' Reparametrize AR Coefficients In Terms of PACF
#' 
#' @description
#' Transform AR parameter coefficients into partial autocorrelation function (PACF).
#' 
#' No check for invertibility is done for maximum computational efficiency
#' since this function is used extensively in the numerical optimization
#' of the AR loglikelihood function in FitAR.
#' Use InvertibleQ to test for invertible AR coefficients.
#' 
#' @usage ARToPacf(phi)
#' @param phi vector of AR parameter coefficients.
#' @details 
#' For details see McLeod and Zhang (2006).
#' @returns Vector of length(phi) containing the parameters in the transformed PACF domain.
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation parameterization
#'   for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [InvertibleQ()], [PacfToAR()].
#' @examples 
#' somePACF<-c(0.5,0.6,0.7,0.8,-0.9,-0.8)
#' # PacfToAR() transforms PACF to AR parameter coefficients.
#' someAR<-PacfToAR(somePACF)
#' test<-ARToPacf(someAR)
#' #This should be very small
#' sum(abs(test-somePACF))
#' 
#' @export
ARToPacf <-
  function(phi){
    phik=phi
    L=length(phi)
    if(L==0) return(0)
    pi=numeric(L)
    for (k in 1:L){
      LL=L+1-k
      a <- phik[LL]
      pi[L+1-k] <- a
      phikp1 <- phik[-LL]
      if(is.na(a) || abs(a)==1)
        stop("transformation is not defined, partial correlation = 1")
      phik <- (phikp1+a*rev(phikp1))/(1-a^2)
    }
    pi
  }
