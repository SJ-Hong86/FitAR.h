#' Select Best AR, ARz or ARp Model
#' 
#' The AIC/BIC/UBIC/EBIC/BICq criterion is used to select the best fitting AR 
#' or subset AR model. When Best > 1, the result may be plotted using `plot`.
#' 
#' @usage SelectModel(z, lag.max = 15, ARModel = c("AR", "ARz", "ARp"), Criterion = "default", Best = 3, Candidates = 5, t="default")
#' @param z time series data.
#' @param lag.max maximum order of autoregression.
#' @param ARModel "AR" for full AR(p) or "ARp"/"ARz" corresponding to subset models.
#' @param Criterion default is "BIC" for order selection and "BICq" for subset selection. 
#'   Options: "AIC", "BIC", "UBIC", "EBIC", "BICq" and "GIC".
#' @param Best final number of models to be selected.
#' @param Candidates number of models initially selected using the approximate criterion.
#' @param t tuning parameter, EBIC, BICq, GIC.
#' @details 
#' McLeod and Zhang (2006) outline an approximate AIC/BIC selection algorithm. 
#' This algorithm is a refinement of that method. The refinement consists of 
#' automatically look for the best k candidates, where k = `Candidates`. Then the 
#' exact likelihood is evaluated for all k candidates. Out of these k candidates, the 
#' best q = `Best` are then selected. This two-step procedure is needed because 
#' if k is too low, the approximate AIC/BIC rankings may not agree with the exact 
#' rankings. This strategy is used for model selection for AR, ARz and ARp models. 
#' A plot method is available to graph the output. The UBIC and EBIC developed 
#' by Chen and Chen (2007) are an extension of the BIC criterion for subset 
#' selection. In the non-subset case UBIC is equivalent to BIC. The EBIC using a 
#' tuning parameter, G, where 0 <= G <= 1. The BICq takes a tuning parameter, 
#' Q, where 0 < Q < 1. The GIC takes a tuning parameter, p, where 0<p<0.25.
#' @returns
#' When `Best`=1, a vector is returned indicated the lag or lags included in the 
#' model. The null model is indicated by returning 0 for the lag. An object with 
#' class "Selectmodel" is returned when `Best` > 1. If ARModel = "AR", a matrix is 
#' return whose first column shows p and second AIC or BIC. Otherwise for 
#' ubset selection, the result is a list with q components, where q = `BEST`. 
#' When Criterion = "UBIC", the components in this list are:
#' 
#' * `p` lags present, a 0 indicates the null model.
#' * `UBIC` exact UBIC.
#' 
#' similarly for the AIC/BIC case. 
#' The components are arranged in order of the criterion used. 
#' When ARModel = "ARp" or "ARz", an attribute "model" indicating "ARp" 
#' "ARz" is included.
#' @note
#' Setting `Candidates` too low can result in anomalous results. For example if 
#' `Candidates` = 1, we find that the top ranking model may depend on how large 
#' `Best` is set. This phenomenon is due to the fact that among the best AIC/BIC 
#' models there is sometimes very little difference in their AIC/BIC scores. 
#' Since the initial ranking of the Candidates is done using the approximate 
#' likelihood, the final ranking using the exact likelihood may change.
#' 
#' For white noise, the best model is the null model, containing no lags. 
#' This is indicating by setting the model order, \eqn{p=0}.
#' @author A.I. McLeod and Y. Zhang.
#' @references
#' McLeod, A.I. and Zhang, Y. (2006). Partial Autocorrelation Parameterization 
#' for Subset Autoregression. Journal of Time Series Analysis, 27, 599-612.
#' 
#' Chen, J. and Chen, Z. (2008). Extended Bayesian Information Criteria for 
#' Model Selection with Large Model Space. Biometrika.
#' @seealso [plot.Selectmodel()], [PacfPlot()], [PacfPlot()], [FitAR()].
#' @examples 
#' # Example 1: Find an ARp subset model for lynx data using BIC
#' z<-log(lynx)
#' out<-SelectModel(z, ARModel="ARp", Criterion="BIC", Best=5)
#' plot(out)
#' 
#' # Example 2: Find an ARz subset model for lynx data using BIC
#' out<-SelectModel(z, ARModel="ARz", Criterion="BIC", Best=5)
#' plot(out)
#' 
#' # Example 3: Select an AR(p) model
#' out<-SelectModel(z, ARModel="AR", Criterion="BIC", Best=5)
#' out
#' plot(out)
#' out<-SelectModel(z, ARModel="AR", Criterion="BIC", Best=1)
#' 
#' # Example 4: Fit subset models to lynx series
#' z<-log(lynx)
#' # requires library leaps. Should be automatically when FitAR package is loaded.
#' # first fit ARp
#' pvec <- SelectModel(z, lag.max=11, ARModel="ARp", Criterion="AIC", Best=1)
#' ans1 <- FitAR(z, pvec, ARModel="ARp", MLEQ=FALSE)
#' # now fit ARz
#' pvec <- SelectModel(z, lag.max=11, ARModel="ARz", Criterion="AIC", Best=1)
#' ans2<-FitAR(z, pvec, ARModel="ARz")
#' # compare
#' summary(ans1)
#' summary(ans2)
#' # Use UBIC
#' pvec <- SelectModel(z, ARModel="ARp",lag.max=11,Best=1)
#' ans3<-FitAR(z, pvec, ARModel="ARp")
#' pvec <- SelectModel(z, ARModel="ARz",lag.max=11,Best=1)
#' ans4<-FitAR(z, pvec, ARModel="ARz")
#' # compare
#' summary(ans3)
#' summary(ans4)
#' 
#' # Example 5: lynx data subset AR models
#' # The AIC and BIC choose the same models as the GIC with t=0.1 and t=0.01 respectively.
#' # An even more parsimonious model is chosen with t=0.001
#' SelectModel(z, lag.max=15, ARModel="ARp", Criterion="GIC", Best=1, Candidates=5, t=0.1)
#' SelectModel(z, lag.max=15, ARModel="ARp", Criterion="GIC", Best=1, Candidates=5, t=0.01)
#' SelectModel(z, lag.max=15, ARModel="ARp", Criterion="GIC", Best=1, Candidates=5, t=0.001)
#' ans<-SelectModel(z, lag.max=15, ARModel="ARp", Criterion="GIC", Best=3, Candidates=5, t=0.001)
#' plot(ans)
#' 
#' @export
SelectModel <-
  function(z, lag.max=15, ARModel=c("AR","ARz","ARp"), Criterion="default", Best=3, Candidates=5, t="default"){
    #AR - order selection problem
    #ARp - AR subset
    #ARz - AR subset, partials
    stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0, Candidates>0)
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    stopifnot(is.wholenumber(lag.max))
    BestCandidates<-Candidates
    IsValidCriterionQ <- Criterion %in% c("default", "AIC", "BIC", "UBIC", "EBIC", "BICq", "GIC")
    if (!IsValidCriterionQ)
      stop("Criterion = ", Criterion, " not known.")
    ARModel <- match.arg(ARModel)
    if (Best > BestCandidates)
      BestCandidates<-Best
    if (ARModel=="ARp") #subset ARp
      return(GetLeapsAR(z, lag.max=lag.max, Criterion=Criterion, Best=Best, Candidates=Candidates, t=t))
    if (ARModel=="ARz")
      SubsetQ <- TRUE
    else
      SubsetQ <- FALSE
    method<-Criterion
    if (Criterion == "default")
      if (SubsetQ)
        method <- "BICq"
    else
      method <- "BIC"
    if (!SubsetQ && Criterion=="UBIC")
      method <- "BIC"
    #set tuning parameter
    P<-0.01
    Q<-0.25
    G<-1
    if (method=="EBIC"  && t!="default")  G <- t
    if (method=="QBIC"  && t!="default")  Q <- t
    if (method=="GIC"   && t!="default")  P <- t
    if (P>=0.25 || P<=0)
      stop("error: GIC tuning parameter invalid")
    if (Q<=0 || Q>=1)
      stop("error: BICq tuning parameter invalid")
    #approximate likelihood, "AR" or "ARz"
    zta<-ARToPacf(ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
    n<-length(z)
    LagRange<-1:lag.max
    if (method=="UBIC"){
      mColNames<-list(c("p", "UBIC-Exact", "UBIC-Approx"))
      PENALTY1 <- log(n) + 2*lchoose(lag.max, 1)
      penalty<-log(n)*(1+LagRange)+2*lchoose(lag.max, LagRange)
    }
    if (method=="EBIC"){
      mColNames<-list(c("p", "UBIC-Exact", "UBIC-Approx"))
      PENALTY1 <- log(n) + 2*G*lchoose(lag.max, 1)
      penalty<-log(n)*(1+LagRange)+2*G*lchoose(lag.max, LagRange)
    }
    if (method=="BICq"){
      mColNames<-list(c("p", "BICq-Exact", "BICq-Approx"))
      PENALTY1 <- log(n) - 2*log(Q/(1-Q))
      penalty<-(1+LagRange)*PENALTY1
    }
    if (method=="BIC"){
      PENALTY1 <- log(n)
      mColNames<-list(c("p", "BIC-Exact", "BIC-Approx"))
      penalty<-(1+LagRange)*PENALTY1
    }
    if (method=="AIC"){
      PENALTY1 <- 2
      mColNames<-list(c("p", "AIC-Exact", "AIC-Approx"))
      penalty<-(1+LagRange)*PENALTY1
    }
    if (method=="GIC"){
      mColNames<-list(c("p", "GIC-Exact", "GIC-Approx"))
      PENALTY1 <- qchisq(p=(1+sqrt(1-4*P))/2, df=1)
      penalty<-(1+LagRange)*PENALTY1
    }
    if (SubsetQ)
      LagsEntering<-order(abs(zta),decreasing=TRUE)
    else  
      LagsEntering<-1:lag.max
    LLapprox <- n*log(cumprod(1-zta[LagsEntering]^2))
    AnIC <- LLapprox + penalty
    #
    IndCandidates<-order(AnIC)[1:BestCandidates]
    AnICexact<-numeric(BestCandidates+1)
    if (SubsetQ){ #subset. AR model subset selection.
      m<-as.list(numeric(BestCandidates+1))
      for (isub in 1:BestCandidates){
        ModelLags<-sort(LagsEntering[1:IndCandidates[isub]])
        LL<-GetFitAR(z-mean(z), ModelLags)$loglikelihood
        k<-length(ModelLags)+1 #mean is included and k>=2 here
        if (method=="UBIC") {
          UBIC <- -2*LL + log(n)*k + 2*lchoose(lag.max+1, k) 
          AnICexact[isub]<-UBIC
          m[[isub]] <- list(p=ModelLags, UBIC=UBIC)
        }
        if (method=="EBIC") {
          EBIC <- -2*LL + log(n)*k + 2*G*lchoose(lag.max+1, k) 
          AnICexact[isub]<-EBIC
          m[[isub]] <- list(p=ModelLags, EBIC=EBIC)
        }
        if (method=="BICq") {
          BICq <- -2*LL + log(n)*k -2*(k*log(Q)+(lag.max+1-k)*log(1-Q))
          AnICexact[isub]<-BICq
          m[[isub]] <- list(p=ModelLags, BICq=BICq)
        }
        if (method=="AIC"){
          AIC <- -2*LL+2*k
          AnICexact[isub]<-AIC
          m[[isub]] <- list(p=ModelLags, AIC=AIC)
        }
        if (method=="BIC") {
          BIC <- -2*LL+log(n)*k
          AnICexact[isub]<-BIC
          m[[isub]] <- list(p=ModelLags, BIC=BIC)
        }
        if (method=="GIC") {
          GIC <- -2*LL+k*qchisq(p=(1+sqrt(1-4*P))/2, df=1)
          AnICexact[isub]<-GIC
          m[[isub]] <- list(p=ModelLags, GIC=GIC)
        }
      }
      #null model, note: k=1
      LL<-GetFitAR(z-mean(z), 0)$loglikelihood
      if (method=="UBIC") {#parameters=1, just mean
        UBIC <- -2*LL + PENALTY1
        AnICexact[BestCandidates+1]<-UBIC
        m[[BestCandidates+1]] <- list(p=0, UBIC=UBIC)
      }
      if (method=="EBIC") {
        EBIC <- -2*LL + PENALTY1
        AnICexact[BestCandidates+1]<-EBIC
        m[[BestCandidates+1]] <- list(p=0, EBIC=EBIC)
      }
      if (method=="BICq") {
        BICq <- -2*LL + PENALTY1 
        AnICexact[BestCandidates+1]<-BICq
        m[[BestCandidates+1]] <- list(p=0, BICq=BICq)
      }
      if (method=="AIC"){
        AIC <- -2*LL+PENALTY1
        AnICexact[BestCandidates+1] <- AIC
        m[[BestCandidates+1]]<-list(p=0,AIC=AIC)
      }
      if (method=="BIC") {
        BIC <- -2*LL+PENALTY1
        AnICexact[BestCandidates+1] <- BIC
        m[[BestCandidates+1]]<-list(p=0,BIC=BIC)
      }
      if (method=="GIC") {
        GIC <- -2*LL+PENALTY1
        AnICexact[isub]<-GIC
        m[[isub]] <- list(p=ModelLags, GIC=GIC)
      }
      #final model select based on exact likelihood
      i<-order(AnICexact)
      m<-m[i]
      m<-m[1:Best] 
      attr(m, "model")<-ARModel              
    }
    else  { #non-subset. AR model order selection.
      AnICexact<-numeric(BestCandidates)
      AnICApprox<-numeric(BestCandidates)
      for (i in 1:BestCandidates){
        p<-LagsEntering[IndCandidates[i]]-1
        AnICApprox[i]<-AnIC[p+1]
        ans<-GetFitAR(z-mean(z), 0:p)
        LL<-ans$loglikelihood
        #mean included in all models, k=p+1
        if (method=="AIC")
          penalty<-2*(p+1)
        if (method=="BIC")
          penalty<-log(n)*(p+1)
        if (method=="BICq")
          penalty<-log(n)*(1+p)-2*((p+1)*log(Q/(1-Q)))
        if (method=="EBIC") #is equivalent to BIC really!
          penalty<-log(n) + 2*G*lchoose(p+1, 1)
        if (method=="UBIC")
          penalty<-log(n) + 2*lchoose(p+1, 1)
        if (method=="GIC")
          penalty<-(1+p)*PENALTY1
        if (method=="BIC")
          penalty<-log(n)*(p+1)      
        AnICexact[i]<- -2*LL+penalty
      }
      m<-c(LagsEntering[IndCandidates]-1,AnICexact,AnICApprox)
      m<-matrix(m,ncol=3)
      if (Best==1) {
        m<-m[order(AnICexact),drop=FALSE]
        m<-m[1:Best,drop=FALSE]
      }
      else {
        m<-m[order(AnICexact),]
        m<-m[1:Best,]
      }
      if (Best > 1)
        dimnames(m)<-c(list(1:Best), mColNames)
    }
    if (SubsetQ) class(m)<-"Selectmodel"
    if (Best > 1)
      m
    else 
      if (is.list(m))
        m[[1]]$p
    else
      as.vector(m[1])
  }
