#' t-statistic for unit root test
#' 
#' Utility function used by UnitRootTest.
#' 
#' @usage getT(ans)
#' @param ans output from FitAR.
#' @returns Value of the test statistic.
#' @author A.I. McLeod.
#' @seealso [getRho()], [UnitRootTest()].
#' @examples
#' z <- cumsum(rnorm(100))
#' ans <- FitAR(z, p=1)
#' getT(ans)
#' 
#' @export
getT <-
  function(ans){
    #T statistic, mean correction case
    (ans$phiHat-1)/sqrt(ans$covHat[1,1])
  }
