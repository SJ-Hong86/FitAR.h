#' Generic Bootstrap Function
#' 
#' @description
#' Generic function to bootstrap a fitted model.
#' 
#' @usage Boot(obj, R=1, ...)
#' @param obj	 fitted object.
#' @param R	 number of bootstrap replications.
#' @param ...	 optional arguments.
#' @details At present, the only function implemented is `Boot.FitAR`.
#' @returns Parametric bootstrap simulation.
#' @author A.I. McLeod and Y. Zhang.
#' @seealso [Boot.FitAR()].
#' @examples 
#' out<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
#' Boot(out)
#' 
#' @export
Boot <-
  function(obj, R=1, ...){
    UseMethod("Boot")
  }
