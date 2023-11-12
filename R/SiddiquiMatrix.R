#' Covariance Matrix of MLE Parameters in an AR(p)
#' 
#' A direct method of computing the inverse of the covariance matrix of p 
#' successive observations in an AR(p) with unit innovation variance given by 
#' Siddiqui (1958) is implemented. This matrix, divided by n = length of series, 
#' is the covariance matrix for the MLE estimates in a regular AR(p). 
#' 
#' @usage SiddiquiMatrix(phi)
#' @param phi coefficients in a regular AR(p).
#' @returns Matrix, covariance matrix of MLE estimates.
#' @note 
#' No check on whether the parameters are in the stationary region is done. It has 
#' been shown a necessary and sufficient condition for the parameters to be in 
#' the stationary region is that this matrix should be positive-definite (Pagano, 
#' 1973). But computationally it is probably better to test for stationarity by using 
#' `ARToPacf` to transform to the PACF and then check that the absolute value of all 
#' partial autocorrelations are less than 1.
#' @author A.I. McLeod.
#' @references 
#' Siddiqui, M.M. (1958) On the inversion of the sample covariance matrix in a 
#' stationary autoregressive process. Annals of Mathematical Statistics 29, 585-588.
#' 
#' Pagano, M. (1973), When is an autoregressive scheme stationary? 
#' Communications in Statistics A 1, 533-544.
#' @seealso [FitAR()].
#' @examples
#' # compute the inverse directly and by Siddiqui's method and compare:
#' phi<-PacfToAR(rep(0.8,5))
#' A<-SiddiquiMatrix(phi)
#' B<-solve(toeplitz(TacvfAR(phi, lag.max=length(phi)-1)))
#' max(abs(A-B))
#' 
#' @export
SiddiquiMatrix <-
  function(phi)
  {
    p <- length(phi)
    phis <- c(-1, phi)
    A <- matrix(numeric(p^2), nrow = p, ncol = p)
    for(j in 1:p) 
      for(i in 1:p) 
        if(j > i) 
          A[i, j] <- A[j, i]
    else {
      k <- 1:min(i, j)
      A[i, j] <- sum(phis[1 + i - k] * phis[1 + j - k] - 
                       phis[1 + p + k - i] * phis[1 + p + k - j])
    }
    A
  }
