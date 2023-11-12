#' Binary representation of non-negative integer
#' 
#' A non-negative integer is represented as a binary number. 
#' The digits, 0 or 1, of this number are returned in a vector.
#' 
#' @usage toBinary(n, k = ceiling(logb(n+1,base=2)))
#' @param n a non-negative integers.
#' @param k number of digits to be returned.
#' @returns A vector of length k. The first element is the least significant digit.
#' @author A.I. McLeod.
#' @examples
#' toBinary(63)
#' toBinary(64)
#' # sometimes we want to pad result with 'leading' 0's
#' toBinary(63, k=20)
#' toBinary(64, k=20)
#' 
#' @export
toBinary <-
  function(n, k=ceiling(logb(n+1, base=2)))
  {
    if (n < 0)
      stop("n must be non-negative integer")
    if (!is.loaded("asBinary"))
      stop("asBinary not loaded")
    m <- numeric(k)
    z<-.C("asBinary", as.integer(n), as.integer(m), as.integer(k))[[2]]
    z   
  }
