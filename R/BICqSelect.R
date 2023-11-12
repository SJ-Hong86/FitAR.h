#' Select best model using BICq
#' 
#' @description
#' Given the loglikelihoods for a set of models arranged in ascending order of
#' size, the best model is selected using the BICq criterion for a specified size.
#' 
#' @usage BICqLL(logL, n, level = 0.99, mSize = 1:length(logL), mComplex = function(k) k)
#' @param logL	 vector of loglikelihoods.
#' @param n	 sample size.
#' @param level	 probability.
#' @param mSize	 model sizes.
#' @param mComplex	 a complexity function.
#' @details 
#' See reference.
#' @returns
#' * `khat` dataframe with columns: k, a.1, a.2 q.1, q.2, level=level, where k is the
#'   optimal model, (a.1,a.2) is the interval for alpha in the GIC, (q.1, q.2) is
#'   the interval for q and level is the probability. Each row corresponds to
#'   an entry in 'level'.
#'   
#' * `table` This table indicates which models can be selected for some values of
#'   alpha or q.    
#' @note AIC corresponds to setting level=0.84. BIC corresponds to setting 
#'   level=pchisq(log(n), 1). So for n=100, 1000; BIC=0.96, 0.97.
#' @author Changjiang Xu and A. Ian McLeod.
#' @references Changjiang Xu and A. Ian McLeod (2010). Bayesian information criterion
#'   Bernoulli prior. Submitted for publication.
#' @seealso [SelectModel()].
#' @examples
#' # Example 1.
#' # AR(p) Order selection for 'lynx' series
#' z <- log(lynx)
#' n <- length(z)
#' lag.max <- 20
#' zta<-ARToPacf(ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
#' LagsEntering<-1:lag.max
#' LLapprox <- (-n)*log(cumprod(1-zta[LagsEntering]^2))
#' ans<-BICqLL(logL=LLapprox, n=n, level=c(0.9, 0.95, 0.99))
#' ans$khat
#' ans$table
#' # if we just want the best model for level=0.99 then,
#' (BICqSelect(logL=LLapprox, n=n, level=0.99)$khat)[[1]]
#' # aic for comparison
#' aic<-(-2*LLapprox)+2*LagsEntering
#' which.min(aic)
#' plot(LagsEntering, aic)
#' 
#' @export
BICqSelect <-
  function(logL, n, level=0.99, mSize=1:length(logL), mComplex=function(k) k)
  {
    #AICa = -2logL(k) + a*C(k)
    #BICq = -2logL(k) + C(k)*[log(n)-2log{q/(1-q)}]
    #C(k) is the model complexity or degree of freedom, usually C(k)=k.
    #logL: log-likelihood
    #n: sample size
    #level: confidence level for controlling overfitting
    #multiple level can be used, such as level=c(0.95, 0.99).
    #mSize: the set of model sizes {k1, k2, ..., k_P}, usually mSize=1:P.
    #Model complex: 
    mC <- mComplex(mSize) 
    #The number of candidate models:
    P<- length(logL)
    stopifnot(P>2)
    stopifnot(length(unique(mSize))==P) #Model sizes are not unique.
    if (!all(mSize[2:P]-mSize[1:(P-1)]>=0) && P>1) stop("Model sizes are not ordered.")
    #The tuning parameter a for AICa:
    aChosen<- qchisq(level,1) 
    ## Ranges of a for AICa 
    a12<- matrix(rep(NA,P*2),ncol=2)
    a12[1,]<- 2*c(max((logL[2:P]-logL[1])/(mC[2:P]-mC[1])),Inf)
    a12[P,]<- 2*c(0, min((logL[1:(P-1)]-logL[P])/(mC[1:(P-1)]-mC[P])))
    for (k in 2:(P-1)){
      i1<-1:(k-1) 
      i2<-(k+1):P
      a12[k,]<- 2*c(max((logL[i2]-logL[k])/(mC[i2]-mC[k])), 
                    min((logL[i1]-logL[k])/(mC[i1]-mC[k])))
    }
    ## Ranges of q for BICq
    q12 <- 1/(1+exp(a12/2)/sqrt(n))
    q12 <- q12[,2:1]
    ## Select the best model
    ks<- rep(NA, length(aChosen))
    as<- matrix(rep(NA, 2*length(aChosen)),ncol=2)
    qs<- matrix(rep(NA, 2*length(aChosen)),ncol=2)
    for (i in 1:length(aChosen)){
      indx<- which(a12[,1]<=aChosen[i] & a12[,2]>=aChosen[i])
      ks[i]<- mSize[indx]
      as[i,]<- a12[indx,]
      qs[i,]<- q12[indx,]
    }
    kaqChosen<- data.frame(k=ks, a=as, q=qs, level=level)
    aqTable<- data.frame(k=mSize, a=a12, q=q12, ms=a12[,1]<=a12[,2])
    colnames(aqTable)[6]="a1<=a2"
    list(kHat=kaqChosen, table=aqTable)
  }

