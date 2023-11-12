#' Coefficients In Infinite Moving Average Expansion
#' 
#' A stationary-causal AR(p) can be written as a general linear process (GLP).
#' This function obtains the moving-average expansion out to the L-th lag,
#' \eqn{z_t=a_t+\psi_1 a_{t-1}+\cdots+\psi_L a_{t-L}}.
#' 
#' @usage ARToMA(phi, lag.max)
#' @param phi	 AR Coefficient vector.
#' @param lag.max	 maximum lag.
#' @details 
#' The coefficients are computed recursively as indicated in Box and Jenkins (1970).
#' @returns Vector of length L+1 containing, \eqn{(1,\psi_1,\dots,\psi_L)}.
#' @author A.I. McLeod and Y. Zhang.
#' @references Box and Jenkins (1970), Time Series Analysis, Forecasting & Control.
#' @seealso [InvertibleQ()].
#' @examples 
#' ARToMA(0.5,20)
#' ARToMA(c(0.2,0.5), 15)
#' 
#' @export
ARToMA <-
  function(phi, lag.max)
  {
    p <- length(phi)
    x <- numeric(lag.max + 1)
    x <- 1
    for(i in 1:p) {
      x[i + 1] <- crossprod(phi[1:i], (rev(x))[1:i])
    }
    if(lag.max > p) {
      for(i in (p + 1):lag.max) {
        x[i + 1] <- crossprod(phi, (rev(x))[1:p])
      }
    }
    return(x)
  }
