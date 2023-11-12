#' Internal Utility Function
#' 
#' Used by Get1G.
#' 
#' @usage GetKappa(phi)
#' @param phi AR coefficients.
#' 
#' @export
GetKappa <-
  function(phi)
    rev(cumsum(rev(TacvfMA(phi, length(phi))[-1])))
