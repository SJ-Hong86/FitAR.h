#' Test if Invertible or Stationary-casual
#' 
#' Tests if the polynomial \deqn{1-\phi (1)B-\cdots-\phi (p)B^p} 
#' where p=length[phi] has all roots outside the unit circle. 
#' This is the invertibility condition for the polynomial.
#' 
#' @usage InvertibleQ(phi)
#' @param phi a vector of AR coefficients.
#' @details The PACF is computed for lags 1, ..., p using eqn. (1) in McLeod 
#'   and Zhang (2006). The invertibility condition is satisfied if and only 
#'   if all PACF values are less than 1 in absolute value.
#' @returns TRUE, if invertibility condition is satisfied. FALSE, if not invertible.
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [ARToPacf()].
#' @examples
#' # simple examples
#' InvertibleQ(0.5)
#' # find the area of the invertible region for AR(2).
#' # We assume that the parameters must be less than 2 in absolute value.
#' # From the well-known diagram in the book of Box and Jenkins (1970), 
#' # this area is exactly 4.
#' NSIM<-10^4
#' phi1<-runif(NSIM, min=-2, max=2)
#' phi2<-runif(NSIM, min=-2, max=2)
#' k<-sum(apply(matrix(c(phi1,phi2),ncol=2), MARGIN=1, FUN=InvertibleQ))
#' area<-16*k/NSIM
#' area
#' 
#' @export
InvertibleQ <-
  function(phi){
    identical(TRUE,try(all(abs(ARToPacf(phi))<1),silent=TRUE))
  }
