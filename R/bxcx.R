#' Box-Cox Transformation and its Inverse
#' 
#' Box-Cox or power transformation or its inverse. For \eqn{\lambda \ne 0}, the Box-Cox 
#' transformation of x is \eqn{(x^{\lambda}-1)/\lambda}, whereas the regular power 
#' transformation is simply \eqn{x^{\lambda}}. When \eqn{\lambda = 0}, it is log in both cases. 
#' The inverse of the Box-Cox and the power transform can also be obtained.
#' 
#' @usage bxcx(x, lambda, InverseQ = FALSE, type = "BoxCox")
#' @param x a vector or time series.
#' @param lambda power transformation parameter.
#' @param InverseQ if TRUE, the inverse transformation is done.
#' @param type either "BoxCox" or "power".
#' @returns A vector or time series of the transformed data.
#' @author A.I. McLeod.
#' @references Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
#'   Journal of Royal Statistical Society, Series B, vol. 26, pp. 211-246.
#' @seealso [BoxCox()].
#' @examples 
#' # lambda=0.5
#' z<-AirPassengers; lambda<-0.5
#' y<-bxcx(z, lambda)
#' z2<-bxcx(y, lambda, InverseQ=TRUE)
#' sum(abs(z2-z))
#' #
#' z<-AirPassengers; lambda<-0.0
#' y<-bxcx(z, lambda)
#' z2<-bxcx(y, lambda, InverseQ=TRUE)
#' sum(abs(z2-z))
#' ## End(Not run)
#' 
#' @export
bxcx <-
  function(x, lambda, InverseQ=FALSE, type="BoxCox"){
    if (type=="BoxCox"){    
      if (!InverseQ) {
        if (min(x) <= 0){
          cat("min data value <= 0, 0.25-min(x) added to data", fill=TRUE)
          x<- x + 0.25-min(x)
        }    
        if (abs(lambda)<1e-6)
          log(x)
        else
          (x^lambda-1)/lambda
      }
      else
        if (abs(lambda)<1e-6)
          exp(x)
      else {
        y <- lambda*x+1
        if (min(y) <= 0) 
          cat("Warning: modified inverse Box-Cox transformation used", fill=TRUE)
        abs(y)^(1/lambda)
      }
    }
    else {    #simple power transformation
      if (!InverseQ) {
        if (min(x) <= 0){
          cat("min data value <= 0, 0.25-min(x) added to data", fill=TRUE)
          x<- x + 0.25-min(x)
        }    
        if (abs(lambda)<1e-6)
          log(x)
        else
          x^lambda
      }
      else
        if (abs(lambda)<1e-6)
          exp(x)
      else {
        y <- x
        if (min(y) <= 0) 
          cat("Warning: modified inverse Box-Cox transformation used", fill=TRUE)
        abs(y)^(1/lambda)
      }
    }
  }
