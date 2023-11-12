#' Transform from PACF Parameters to AR Coefficients
#' 
#' Transforms AR partical autocorrelation function (PACF) parameters to AR 
#' coefficients based on the Durbin-Levinson recursion.
#' 
#' @usage PacfToAR(zeta)
#' @param zeta vector of AR PACF parameters.
#' @details See Mcleod and Zhang (2006).
#' @returns Vector of AR coefficients.
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [InvertibleQ()], [PacfToAR()].
#' @examples
#' somePACF<-c(0.5,0.6,0.7,0.8,-0.9,-0.8)
#' someAR<-PacfToAR(somePACF)
#' test<-ARToPacf(someAR)
#' # this should be very small
#' sum(abs(test-somePACF))
#' 
#' @export
PacfToAR <-
  function(zeta){
    L=length(zeta)
    if (L==0) return(numeric(0))
    if (L==1) return(zeta)
    phik=zeta[1]
    for (k in 2:L){
      phikm1=phik
      phik=c(phikm1-zeta[k]*rev(phikm1),zeta[k])
    }
    phik
  }
