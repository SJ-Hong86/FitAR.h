#' Select lags for Best Subset ARp Model
#' 
#' The subset ARp model is the usual subset model, for example see Tong (1977). 
#' This function is used by SelectModel for model identification for ARp models.
#' 
#' @usage GetLeapsAR(z, lag.max = 15, Criterion = "UBIC", Best = 3, Candidates=5, t="default", ExactQ=FALSE)
#' @param z ts object or vector containing time series.
#' @param lag.max maximum order of the AR.
#' @param Criterion default UBIC, other choices are "AIC", "BIC", "EBIC", "BICq", "GIC".
#' @param Best the number of based selected. Ignore with "GIC".
#' @param Candidates number of models initially selected using the approximate criterion.
#' @param t tuning parameter, EBIC, BICq, GIC.
#' @param ExactQ exhaustive numeration using exact likelihood. Still under under development. NOT AVAILABLE IN THIS VERSION.
#' @details
#' The R function leaps in the R package leaps is used to compute the subset 
#' regression model with the smallest residual sum of squares containing 1, ..., 
#' lag.max parameters. The mean is always included, so the only parameters 
#' considered are the phi coefficients. After the best models containing 1, ..., 
#' lag.max parameters are selected the models are individually refit to determine 
#' the exact likelihood function for each selected model. Based on this likelihood 
#' the UBIC/BIC/AIC is computed and then the best models are selected. The UBIC 
#' criterion was developed by Chen and Chen (2007). The EBIC using a tuning 
#' parameter, G, where 0 <= G <= 1. The BICq takes a tuning parameter, Q, where 
#' 0 < Q < 1. The GIC takes a tuning parameter, p, where 0<p<0.25.
#' @returns
#' When 'Criterion' is one of UBIC, AIC, BIC, EBIC, BICq, a list with components:
#' 
#' * `p` lags present in model.
#' * `UBIC` approximate UBIC (Chen & Chen, 2007), if Criterion=="UBIC".
#' * `AIC` approximate AIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="AIC".
#' * `BIC` approximate BIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="BIC".
#' * `EBIC` approximate EBIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="EBIC".
#' * `BICq` approximate BICq, if Criterion=="BICq".
#' * `GIC` approximate GIC, if Criterion=="GIC".
#' @note 
#' AIC and BIC values produced are not comparable to AIC and BIC produced by 
#' SelectModel for ARz models. However comparable AIC/BIC values are produced 
#' when the selected models are fit by FitAR.
#' 
#' Requires leaps package. Since the least-squares is used, the number of 
#' observations depends on 'lag.max'. Hence different subsets may be chosen 
#' depending on the 'lag.max. See example below.
#' @author A.I. McLeod.
#' @references Tong, H. (1977) Some comments on the Canadian lynx data. 
#'   Journal of the Royal Statistical Society A 140, 432-436.
#' @seealso [SelectModel()], [GetFitARpLS()], [leaps()].
#' @examples
#' # Example 1: Simple Example
#' # for the log(lynx) Tong (1977) selected an ARp(1,2,4,10,11)
#' # using the AIC and a subset selection algorithm. Our more exact
#' # approach shows that the ARp(1,2,3,4,10,11) has slightly lower
#' # AIC (using exact likelihood evaluation).
#' z<-log(lynx)
#' GetLeapsAR(z, lag.max=11)
#' GetLeapsAR(z, lag.max=11, Criterion="BIC")
#' 
#' # Example 2: Subset autoregression depends on lag.max!
#' # Because least-squares is used, P=lag.max observations are deleted.
#' # This phenomenon does not happen with "ARz" subset models
#' # ARp models depend on lag.max
#' GetLeapsAR(z, lag.max=15, Criterion="BIC")
#' GetLeapsAR(z, lag.max=20, Criterion="BIC")
#' 
#' # Example 3: Comparing GIC with BIC, AIC, UBIC and BICq
#' z <- log(lynx)
#' GetLeapsAR(z, lag.max=15, Criterion="BIC", Best=1)
#' GetLeapsAR(z, lag.max=15, Criterion="AIC", Best=1)
#' GetLeapsAR(z, lag.max=15, Criterion="UBIC", Best=1)
#' GetLeapsAR(z, lag.max=15, Criterion="BICq", Best=1, t=0.25)
#' GetLeapsAR(z, lag.max=15, Best=1, Criterion="GIC", t=0.01)
#' ans<-GetLeapsAR(z, lag.max=15, Best=3, Criterion="GIC", t=0.001)
#' plot(ans)
#' 
#' @export
GetLeapsAR <-
  function(z, lag.max=15, Criterion="UBIC", Best=3, Candidates=5, t="default", ExactQ=FALSE){
    stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0)
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    stopifnot(is.wholenumber(lag.max))
    method<-Criterion
    if (is.na(pmatch(method,c("UBIC","AIC","BIC","EBIC","BICq","GIC"))))
      method<-"UBIC"
    if (ExactQ){
      stop("Sorry this option is not available yet!")
      M<-2^lag.max
      logL <- numeric(M)
      logL[1] <- GetFitARpMLE(z, pvec=0)$loglikelihood
      for (i in 1:(M-1)){
        ind<-as.logical(rev(toBinary(i, lag.max)))
        pvec <- (1:lag.max)[ind]
        out<-GetFitARpMLE(z, pvec=pvec)
        logL[i+1] <- out$loglikelihood
      }
    }
    #set tuning parameter
    P<-0.01
    Q<-0.25
    G<-1
    if (method=="EBIC"  && t!="default")  G <- t
    if (method=="QBIC"  && t!="default")  Q <- t
    if (method=="GIC"   && t!="default")  P <- t
    if (P>=0.25 || P<=0)
      stop("error: GIC tuning parameter invalid")
    #level <- 1 - P
    if (Q<=0 || Q>=1)
      stop("error: BICq tuning parameter invalid")
    #GIC/BICq treated as a special case
    pvec <- 1:lag.max
    LagRange<-1:lag.max
    n <- length(z)-lag.max
    ind <- (lag.max+1):length(z)
    y<-z[ind]
    X<-matrix(rep(0,n*lag.max), nrow=n, ncol=lag.max)
    for (i in 1:lag.max)
      X[,i] <- z[ind-pvec[i]]
    outLeaps <- leaps(y=y,x=X,nbest=1,method="r2", strictly.compatible=FALSE)
    k <- outLeaps$size
    #
    #this defines an approximate likelihood approach
    TotSS <- sum((y-mean(y))^2)
    RSS <- TotSS*(1-outLeaps$r2)
    LogL <- (-n/2)*log(RSS/n)
    if (method=="AIC")
      ic<- -2*LogL + 2*k
    if (method=="BIC")
      ic<- -2*LogL + log(n)*k
    if (method=="UBIC")
      ic<- -2*LogL + log(n)*k + 2*lchoose(lag.max+1, k)
    if (method=="EBIC")
      ic<- -2*LogL + log(n)*(1+LagRange)+2*G*lchoose(lag.max, k)
    if (method=="BICq")
      ic<- -2*LogL + log(n)*(1+LagRange)-2*(LagRange*log(Q)+(lag.max+1-k)*log(1-Q))
    if (method=="GIC")
      ic<- -2*LogL + k*qchisq(p=(1+sqrt(1-4*P))/2, df=1) 
    indBest<-order(ic)
    #extra step needed because leaps does not include null model
    #Very important: refit with exact MLE
    LogL<-numeric(Candidates+1)
    #LogL[1]<-GetFitARpLS(z, 0)$loglikelihood #null model included here
    LogL <-GetFitAR(z, 0)$loglikelihood #null model included here
    for (i in 1:Candidates)
      LogL[i+1]<-GetFitAR(z,pvec[outLeaps$which[indBest[i],]])$loglikelihood
    #LogL[i+1]<-GetFitARpLS(z,pvec[outLeaps$which[indBest[i],]])$loglikelihood
    k<-c(1,k[indBest[1:Candidates]])
    if (method=="AIC")
      ic<- -2*LogL + 2*k
    if (method=="BIC")
      ic<- -2*LogL + log(n)*k
    if (method=="EBIC")
      ic<- -2*LogL + log(n)*k + 2*G*lchoose(lag.max+1, k)
    if (method=="UBIC")
      ic<- -2*LogL + log(n)*k + 2*lchoose(lag.max+1, k)
    if (method=="BICq")
      ic<- -2*LogL + k*(log(n) - 2*log(Q/(1-Q)))
    if (method=="GIC")
      ic<- -2*LogL + k*qchisq(p=(1+sqrt(1-4*P))/2, df=1)
    icBest<-order(ic)[1:Best] #best models exact
    m<-as.list(numeric(Best))
    for (i in 1:Best){
      ind<-icBest[i]
      if (indBest[ind] ==1)
        p<-0
      else
        p<-pvec[outLeaps$which[indBest[ind-1],]]
      if (method == "AIC")
        m[[i]] <-list(p=p, AIC=ic[ind])
      if (method == "BIC")
        m[[i]] <-list(p=p, BIC=ic[ind])
      if (method == "UBIC")
        m[[i]] <-list(p=p, UBIC=ic[ind])
      if (method == "EBIC")
        m[[i]] <-list(p=p, EBIC=ic[ind], g=G)
      if (method == "BICq")
        m[[i]] <-list(p=p, BICq=ic[ind], q=Q)
      if (method == "GIC")
        m[[i]] <-list(p=p, GIC=ic[ind], p=P)
    }
    class(m)<-"Selectmodel"
    attr(m,"model")<-"ARp"
    if (Best>1)
      ans<-m
    else {
      ans <- m[[1]]$p
      if (length(ans)==0)
        ans<-0
    }
    ans
  }
