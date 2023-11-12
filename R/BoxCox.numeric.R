#' Box-Cox Analysis for a Time Series
#' 
#' An AR(p) model is selected using AIC and then the best Box-Cox transformation 
#' is determined. Requires package FitAR.
#' 
#' @usage BoxCox(object, interval = c(-1, 1), IIDQ = FALSE, ...) ## S3 method for class 'numeric'
#' @param object a vector of time series values.
#' @param interval interval to be searched.
#' @param IIDQ If true, IID is assumed, ie. p=0. If FALSE, AR(p) is fit with 
#'   p determined using AIC.
#' @param ...	 optional arguments. 
#' @details 
#' For \eqn{\lambda \ne 0}, the Box-Cox transformation is of x is \eqn{(x^{\lambda}-1)/\lambda)}.
#' 
#' If the minimum data value is <= 0, a small positive constant, equal 
#' to the negative of the minimum plus 0.25, is added to all the data values. 
#' If length(object) < 20, no AR model is used, that is, p=0.
#' 
#' @returns No value returned. Graphical output is produced as side-effect.
#'   The plot shows relative likelihood function as well as the MLE and a confidence interval.
#' @note The MASS package has a similar function `boxcox` 
#'   but this is implemented only for regression and analysis of variance.    
#' @author A.I. McLeod and Y. Zhang.
#' @references Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
#'   Journal of Royal Statistical Society, Series B, vol. 26, pp. 211-246.
#' @seealso [BoxCox.FitAR()], [BoxCox.Arima()], [BoxCox.ts()].
#' @examples 
#' ## Not run:  # takes a few seconds
#' # annual sunspot series
#' BoxCox(sunspot.year, IIDQ=FALSE)
#' #
#' # non-time series example, lengths of rivers
#' BoxCox(rivers)
#' ## End(Not run)
#' 
#' @export
BoxCox.numeric <-
  function(object, interval=c(-1,1), IIDQ = FALSE, ...){
    n<-length(object)
    if (n <= 20 || IIDQ)
      p<-0
    else 
      p<-SelectModel(object, lag.max=min(30,round(length(object)/4)), Best=1)
    z<-as.vector(object)
    if (min(z) <= 0){
      cat(" minimum data value <= 0 so -min+0.25 added to all values", fill=TRUE)
      z<- z + 0.25 - min(z)
    }
    #now get Jacobian
    J <- sum(log(z))
    #definite loglikelihood function
    LogL<-function(lam){
      y<-bxcx(z,lam)
      y<-y-mean(y)
      out<-GetFitAR(y,p)
      (out$loglikelihood) + (lam-1)*J
    }
    #optimize
    ans<-optimize(LogL,interval=interval,maximum=TRUE)
    lamHat<-ans$maximum
    LLlamHat<-ans$objective
    #determine right and left ends for plotting RL
    RL<-1
    rightLam<-lamHat
    while (RL>0.01){
      rightLam<-rightLam+0.1
      RL<-exp(LogL(rightLam)-LLlamHat)
    }
    RL<-1
    leftLam<-lamHat
    while (RL>0.01){
      leftLam<-leftLam-0.1
      RL<-exp(LogL(leftLam)-LLlamHat)
    }
    lams<-seq(leftLam,rightLam,length=21)
    LL<-numeric(length(lams))
    for (i in 1:length(lams))
      LL[i]<-exp(LogL(lams[i])-LLlamHat)
    ans<-spline(lams,LL)
    dataTI<-attr(data,"title")
    plot(ans$x, ans$y, type="l", xlab=expression(lambda), ylab=expression("R("*lambda*")"), main="Relative Likelihood Analysis\n95% Confidence Interval",
         sub=dataTI)
    abline(h=0.1465, col="blue", lwd=2)
    text((lamHat+rightLam)/1.85,0.8,bquote(hat(lambda)==.(round(lamHat,3))))
  }
