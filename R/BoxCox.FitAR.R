#' Box-Cox Analysis for "FitAR" Objects
#' 
#' This is a methods function to do a Box-Cox analysis for models fit using FitAR.
#' 
#' @usage BoxCox(object, interval = c(-1, 1), type = "BoxCox", InitLambda = "none", ...) ## S3 method for class 'FitAR'
#' @param object output from FitAR.
#' @param interval interval to be searched for the optimal transformation.
#' @param type ignored unless, InitLambda!="none". Type of transformation,
#'   default is "BoxCox". Otherwise a simple power transformation.
#' @param InitLambda default "none". Otherwise a numerical value
#'   giving the transformation parameter.
#' @param ...	 optional arguments passed to optimize. 
#' @details 
#' If no transformation is used on the data, then the original data is used. 
#' But if a transformation has already been used, we need to inverse 
#' transform the data to recover the untransformed data.
#' 
#' For \eqn{\lambda \ne 0}, the Box-Cox transformation is of x is \eqn{(x^{\lambda}-1)/\lambda)}.
#' If the minimum data value is <= 0, a small positive constant, equal 
#' to the negative of the minimum plus 0.25, is added to all the data values.
#' 
#' @returns No value returned. Graphical output is produced as side-effect.
#'   The plot shows relative likelihood function as well as the MLE and a confidence interval.
#' @note The MASS package has a similar function `boxcox` 
#'   but this is implemented only for regression and analysis of variance.    
#' @author A.I. McLeod.
#' @references Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
#'   Journal of Royal Statistical Society, Series B, vol. 26, pp. 211-246.
#' @seealso [BoxCox()], [BoxCox.Arima()].
#' @examples 
#' ## Not run:  # takes a few seconds
#' # lynx time series. ARp subset model.
#' out<-FitAR(lynx, c(1,2,4,10,11), ARModel="ARp")
#' BoxCox(out)
#' #
#' p<-SelectModel(lynx, ARModel="ARz", lag.max=25, Best=1)
#' out<-FitAR(lynx, p)
#' BoxCox(out)
#' ## End(Not run)
#' 
#' @export
BoxCox.FitAR <-
  function(object, interval=c(-1,1), type="BoxCox", InitLambda="none", ...){
    if (InitLambda=="none")
      initLQ <- FALSE
    else
      initLQ <- TRUE
    if (object$ARModel=="ARp") 
      ycall<-"GetFitARpLS(z=y, p=pBXCX)"
    if (object$ARModel=="ARz") 
      ycall<-"GetFitARz(z=y, p=pBXCX)"
    zBXCX<-object$z
    pBXCX <- object$pvec
    if (initLQ)
      zBXCX<-bxcx(zBXCX, InitLambda, type=type, InverseQ=TRUE)
    if (min(zBXCX) <= 0){
      cat(" minimum data value <= 0 so -min+0.25 added to all values", fill=TRUE)
      data<- data + 0.25 - min(data)
    }
    #now get Jacobian
    J <- sum(log(zBXCX))
    #definite loglikelihood function
    LogL<-function(lam){
      y<-bxcx(zBXCX,lam)
      y<-y-mean(y)
      out<-eval(parse(text=ycall))
      (out$loglikelihood) + (lam-1)*J
    }
    #optimize
    ans<-optimize(LogL,interval=interval,maximum=TRUE)
    lamHat<-ans$maximum
    LLlamHat<-ans$objective
    #determine right and left ends for plotting RL
    RL<-1
    rightLam<-lamHat
    MaxIt<-10
    iter<-0
    while (RL>0.01 && iter<MaxIt){
      rightLam<-rightLam+0.1
      RL<-exp(LogL(rightLam)-LLlamHat)
      iter<-iter+1
    }
    RL<-1
    leftLam<-lamHat
    MaxIt<-10
    iter<-0
    while (RL>0.01 && iter<MaxIt){
      leftLam<-leftLam-0.1
      RL<-exp(LogL(leftLam)-LLlamHat)
      iter<-iter+1
    }
    lams<-seq(leftLam,rightLam,length=21)
    LL<-numeric(length(lams))
    for (i in 1:length(lams)){
      LL[i]<-exp(LogL(lams[i])-LLlamHat)
    }
    ans<-spline(lams,LL)
    dataTI<-attr(data,"title")
    plot(ans$x, ans$y, type="l", xlab=expression(lambda), ylab=expression("R("*lambda*")"), main="Relative Likelihood Analysis\n95% Confidence Interval",
         sub=dataTI)
    abline(h=0.1465, col="blue", lwd=2)
    text((lamHat+rightLam)/1.85,0.8,bquote(hat(lambda)==.(round(lamHat,3))))
  }
