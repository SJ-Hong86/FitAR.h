#' Box-Cox Analysis for "Arima" Objects
#' 
#' Implements Box-Cox analysis for "Arima" class objects, the output from `arima`,
#' a R built-in function. Variance change in time series is an important topic.
#' In some cases using a Box-Cox transformation will provide a much simpler
#' analysis than the much more complex ARMA-GARCH approach. 
#' See US Tobacco series example given below for an example.
#' 
#' @usage BoxCox(object, interval = c(-1, 1), type = "BoxCox", InitLambda = "none", ...)
#' @param object output from arima, a R built-in function.
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
#' The log of the Jacobian is \eqn{(\lambda-1)\sum_{t=D+1}^{n}\mathrm{log}(z_t)}, where 
#' \eqn{\lambda} is the transformation, n=length(z), z is the vector of data and 
#' D=d+ds*s, where d is the degree of regular differencing, ds is the degree 
#' of seasonal differencing and s is the seasonal period.
#' The correct expression for the loglikelihood function was first 
#' given in Hipel and McLeod (1977, eqn. 10). Using the wrong expression 
#' for the Jacobian has a disasterous effect in many situations. For example with 
#' the international airline passenger time series, the MLE for lambda would 
#' be about 1.958 instead of close to zero.
#' 
#' If the minimum data value is <= 0, a small positive constant, equal to the 
#' negative of the minimum plus 0.25, is added to all the data values.
#' @returns No value returned. Graphical output is produced as side-effect.
#'   The plot shows relative likelihood function as well as the MLE and a confidence interval.
#' @note The MASS package has a similar function `boxcox` 
#'   but this is implemented only for regression and analysis of variance.    
#' @author A.I. McLeod and Y. Zhang.
#' @references Hipel, K.W. and McLeod, A.I. (1977). Advances in Box-Jenkins Modelling. 
#'   Part 1, Model Construction. Water Resources Research 13, 567-575.
#' @seealso [arima()], [BoxCox()], [BoxCox.FitAR()].
#' @examples 
#' ## Not run:  # not run to save time!
#' # Tobacco Production
#' plot(USTobacco)
#' USTobacco.arima<-arima(USTobacco,order=c(0,1,1))
#' BoxCox(USTobacco.arima)
#' #
#' air.arima<-arima(AirPassengers, c(0,1,1), seasonal=list(order=c(0,1,1), period=12))
#' BoxCox(air.arima)
#' #
#' # In this example, we fit a model to the square-root of the sunspots and
#' # back transform in BoxCox.
#' sqrtsun.arima<-arima(sqrt(sunspot.year),c(2,0,0))
#' BoxCox(sqrtsun.arima, InitLambda=0.5, type="power")
#' #
#' # Back transform with AirPassengers
#' Garima<-arima(log(AirPassengers), c(0,1,1), seasonal=list(order=c(0,1,1),period=12))
#' BoxCox(Garima, InitLambda=0)
#' ## End(Not run)
#' 
#' @export
BoxCox.Arima <-
  function(object, interval=c(-1,1), type="BoxCox", InitLambda="none", ...){
    if (InitLambda=="none")
      initLQ <- FALSE
    else
      initLQ <- TRUE
    xcall<-deparse(object$call,width.cutoff=500)
    #replace x = ? with x=y
    ycall<-sub("x = .*, o", "x=y, o", xcall, perl=TRUE)
    #get data
    test<-gsub("arima(x =","",xcall,fixed=TRUE)
    data<-gsub(",.*)","",test)
    data<-eval(parse(text=data))
    if (initLQ)
      data<-bxcx(data, InitLambda, type=type, InverseQ=TRUE)
    if (min(data) <= 0){
      cat(" minimum data value <= 0 so -min+0.25 added to all values", fill=TRUE)
      data<- data + 0.25 - min(data)
    }
    #parse object to determine d, ds, s
    Gcall<-deparse(object$call,width.cutoff=500)
    modelorder<-sub("^.*, order = ", "", Gcall, perl=TRUE)
    modelorder<-sub("\\, sea.*", "", modelorder, perl=TRUE )
    modelorder<-sub("))",")",modelorder)
    modelorder<-eval(parse(text=modelorder))
    d<-modelorder[2]
    ds<-s<-0
    ind<-grep("seasonal", Gcall)
    if (length(ind) > 0){
      test<-sub("^.*list","list",Gcall)
      FirstRB<-c(regexpr(")",test))
      SecondRB<-c(regexpr(")",substring(test,FirstRB+1,nchar(test))))
      test<-substring(test,1,FirstRB+SecondRB)
      seaord<-eval(parse(text=test))
      ds<-seaord[[1]][2]
      s<-seaord$period
    }
    D<-d+s*ds
    #now get Jacobian
    J <- sum(log(data[(D+1):length(data)]))
    # For comparison, this is WRONG: J <- sum(log(data))
    #definite loglikelihood function
    LogL<-function(lam){
      y<-bxcx(data,lam)
      out.arima<-eval(parse(text=ycall))
      (out.arima$loglik) + (lam-1)*J
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
    for (i in 1:length(lams)){
      LL[i]<-exp(LogL(lams[i])-LLlamHat)
    }
    ans<-spline(lams,LL)
    dataTI<-attr(data,"title")
    plot(ans$x, ans$y, type="l", xlab=expression(lambda), ylab=expression("R("*lambda*")"), main="Relative Likelihood Analysis\n95% Confidence Interval",
         sub=dataTI)
    abline(h=0.1465, col="blue", lwd=2)
    text((lamHat+rightLam)/1.8,0.8,bquote(hat(lambda)==.(round(lamHat,3))))
  }
