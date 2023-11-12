#' Display Estimated Parameters from Output of "FitAR"
#' 
#' Method function to display fitted parameters, their standard errors and 
#' Z-ratio for AR models fit with FitAR.
#' 
#' @usage coef(object, ...) ## S3 method for class 'FitAR'
#' @param object obj the output from FitAR.
#' @param ... optional parameters.
#' @returns A matrix is returned. The columns of the matrix are labeled MLE, sd and Z-ratio.
#'   The rows labels indicate the AR coefficients which were estimated followed by 
#'   mu, the estimate of mean.
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation parameterization 
#'   for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @examples 
#' # Fit subset AR to SeriesA
#' outA<-FitAR(SeriesA, c(1,2,7), ARModel="ARz")
#' coef(outA)
#' #
#' outALS<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
#' coef(outALS)
#' 
#' @export
coef.FitAR <-
  function(object, ...){
    if (object$SubsetQ)
      if (object$FitMethod=="MLE")
        BETA<-object$zetaHat
    else
      BETA<-(object$phiHat)[object$pvec]
    else
      BETA<-object$phiHat
    p<-length(BETA)
    BETA<-c(BETA,object$muHat)
    sdB<-sqrt(diag(object$covHat))
    sdB<-c(sdB,sqrt((object$sigsq)/(length(object$res)*sum(c(1,-object$phiHat))^2)))
    Z<-BETA/sdB
    if (object$SubsetQ)
      if (object$FitMethod=="MLE")
        rn<-c(paste("zeta(",object$pvec,")", sep=""),"mu")
    else
      rn<-c(paste("phi(",object$pvec,")", sep=""),"mu")
    
    else
      rn<-c(paste("phi(",1:p,")", sep=""),"mu")
    cn<-c("MLE","sd","Z-ratio")
    ans<-matrix(c(BETA,sdB, Z),ncol=3)
    dimnames(ans)<-list(rn, cn)
    ans
  }
