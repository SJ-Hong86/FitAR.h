#' Normalized rho unit root test statistic
#' 
#' Utility function used by UnitRootTest.
#' 
#' @usage getRho(ans)
#' @param ans output from FitAR.
#' @returns Value of the test statistic.
#' @author A.I. McLeod.
#' @seealso [getT()], [UnitRootTest()].
#' @examples
#' z <- cumsum(rnorm(100))
#' ans <- FitAR(z, p=1)
#' getRho(ans)
#' 
#' @export
getRho <-
  function(ans){
    #rho test-statistic
    phiHat <- ans$phiHat
    res <- ans$res
    n <- length(res)
    n*(phiHat-1)
  }
