#' Exact MLE Mean-Zero AR(1)
#' 
#' @description
#' This function is used by GetFitAR in the AR(1) case.
#' It is a fast exact solution using the root of a cubic equation.
#' 
#' @usage AR1Est(z, MeanValue = 0)
#' @param z 	time series or vector.
#' @param MeanValue  known mean.
#' @details 
#' The exact MLE for mean-zero AR(1) satisfies a cubic equation.
#' The solution of this equation for the MLE given by Zhang (2002) is used.
#' This approach is more reliable as well as faster than
#' the usual approach to the exact MLE using a numerical optimization technique
#' which can occasionally have convergence problems.
#' @returns MLE for the parameter.
#' @author A.I. McLeod and Y. Zhang.
#' @references Zhang, Y. (2002). Topics in Autoregression,
#'   Ph.D. Thesis, University of Western Ontario.
#' @seealso [GetFitARz()]
#' @examples 
#' AR1Est(lynx-mean(lynx))
#' @export
AR1Est <-
  function(z, MeanValue=0){
    stopifnot(length(z)>1)
    n=length(z)
    m <- MeanValue
    mz <- rep(m,n)
    a <- sum((z-mz)^2)
    b <- sum((z[-1]-mz[-1])*(z[-n]-mz[-n]))
    c <- sum((z[c(-1,-n)]-mz[c(-1,-n)])^2)
    i <- complex(1,0,1)
    x <- ((-16)*b^3+18*a*b*c+24*b^3*n-27*a*b*c*n-9*b*c^2*n-12*b^3*n^2+9*a*b*c*n^2+27*b*c^2*n^2+2*b^3*n^3-18*b*c^2*n^3)
    y <- (-b^2*(-2+n)^2+3*c*(-1+n)*(-a-c*n))
    f <- complex(1,x^2+4*y^3,0)
    z <- (x+sqrt(f))^(1/3)
    g <- x^2+4*y^3
    z1 <- (x+(-g)^(1/2)*i)^(1/3)
    part1 <- (n-2)*b/(3*c*(n-1))
    part2 <- (1-sqrt(3)*i)*y/(3*2^(2/3)*c*(n-1)*z)
    part3 <- (1+sqrt(3)*i)*z/(6*2^(1/3)*c*(n-1))
    Re(part1+part2-part3) 
  }
