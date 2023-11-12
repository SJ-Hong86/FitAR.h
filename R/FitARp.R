#' Fit subset ARp Models
#' 
#' The subset ARp is defined as an AR(p) in which some of the ar-coefficients 
#' are constrained to zero. This is the usual type of subset AR. In contrast 
#' the ARz model constrains some of the partial autocorrelation coefficients to zero.
#' 
#' @usage FitARp(z, p, lag.max = "default", MLEQ = FALSE)
#' @param z time series, vector or ts object.
#' @param p p specifies the model. If length(p) is 1, an AR(p) is assumed and 
#'   if p has length greater than 1, a subset ARp is assumed. For example, 
#'   to fit a subset model with lags 1 and 4 present, set p to c(1,4) or 
#'   equivalently c(1,0,0,4). To fit a subset model with just lag 4, 
#'   you must use p=c(0,0,0,4) since p=4 will fit a full AR(4).
#' @param lag.max the residual autocorrelations are tabulated for lags 
#'   1, ..., lag.max. Also lag.max is used for the Ljung-Box portmanteau test.
#' @param MLEQ TRUE, use MLE. FALSE, use LS.
#' @details 
#' Subset ARp model is fit using exact MLE. The built-in arima function is used 
#' for MLE. When MLEQ=FALSE, LS is used. LS is has been widely used in past for 
#' subset ARp fiting.
#' @returns
#' A list with class name "FitAR" and components:
#' 
#' * `loglikelihood` value of the loglikelihood.
#' * `phiHat` coefficients in AR(p) – including 0's.
#' * `sigsqHat` innovation variance estimate.
#' * `muHat` estimate of the mean.
#' * `covHat` covariance matrix of the coefficient estimates.
#' * `zetaHat` transformed parameters, length(zetaHat) = \# coefficients estimated.
#' * `RacfMatrix` residual autocorrelations and sd for lags 1, ..., lag.max.
#' * `LjungBox` table of Ljung-Box portmanteau test statistics.
#' * `SubsetQ` parameters in AR(p) – including 0's.
#' * `res` innovation residuals, same length as z.
#' * `fits` fitted values, same length as z.
#' * `pvec` lags used in AR model.
#' * `demean` TRUE if mean estimated otherwise assumed zero.
#' * `FitMethod` "MLE" or "LS".
#' * `IterationCount` number of iterations in mean mle estimation.
#' * `convergence` value returned by optim – should be 0.
#' * `MLEMeanQ` TRUE if mle for mean algorithm used.
#' * `ARModel` "ARp" if FitARp used, otherwise "ARz".
#' * `tsp` tsp(z).
#' * `call` result from match.call() showing how the function was called.
#' * `ModelTitle` description of model.
#' * `DataTitle` returns attr(z,"title").
#' * `z` time series data input.
#' @author A.I. McLeod.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()], [FitARz()], [GetFitARz()], [GetFitARpMLE()], [RacfPlot()].
#' @examples 
#' # First example: Fit to AR(4) 
#' set.seed(3323)
#' phi<-c(2.7607,-3.8106,2.6535,-0.9238)
#' z<-SimulateGaussianAR(phi,1000)
#' # MLE using arima
#' ans1<-FitARp(z,4,MLEQ=TRUE)
#' ans1
#' coef(ans1)
#' # OLS
#' ans2<-FitARp(z,4,MLEQ=FALSE)
#' ans2
#' coef(ans2)
#' 
#' ## Not run: # save time building package
#' # Second Example: Fit subset ARp model
#' z<-log(lynx)
#' # MLE 
#' FitARp(z, c(1,2,4,7,10,11),MLEQ=TRUE)
#' # LS
#' FitARp(z, c(1,2,4,7,10,11),MLEQ=FALSE)
#' 
#' # Third Example: Use UBIC model selection to fit subset models
#' z<-log(lynx)
#' p<-SelectModel(z,ARModel="ARp")[[1]]$p
#' # MLE # error returned by arima
#' # ans1<-FitARp(z, p, MLEQ=TRUE)
#' # ans1
#' # LS
#' ans2<-FitARp(z, p, MLEQ=FALSE)
#' ans2
#' 
#' @export
FitARp <-
  function(z,p,lag.max="default",MLEQ=FALSE)
  {
    stopifnot(length(z)>0, length(z)>length(p), length(p)>0)
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    stopifnot(is.wholenumber(p), p>0)
    n<-length(z)
    if (lag.max=="default")
      MaxLag <- min(300, ceiling(length(z)/5))
    else
      MaxLag <- lag.max
    pvec <- sort(p)
    pvec<-pvec[pvec>0]
    if (length(pvec)==0)
      pvec<-0
    PMAX<-max(pvec)
    if (length(p)==1 && PMAX!=0)  
      pvec <- 1:p
    SubsetQ <- length(pvec)<PMAX
    if (PMAX == 0) SubsetQ<-FALSE
    mz <- mean(z) 
    y <- z 
    #get parameter estimates
    if (MLEQ){
      ans<-GetFitARpMLE(y,pvec)
      FitMethod<-"MLE"
      MeanMLEQ<-TRUE    
    }
    else {
      ans<-GetFitARpLS(y,pvec)
      FitMethod<-"LS"
      MeanMLEQ<-FALSE
    }
    phiHat<-ans$phiHat
    res<-BackcastResidualsAR(y-mz, phiHat)
    fits<-y-res
    sigsq<-sum(res^2)/n
    racf<-(acf(res, plot=FALSE, lag.max=MaxLag)$acf)[-1]
    #covariance matrix via inverse Fisher information matrix
    #sd of racf
    if (SubsetQ){
      varNames<-paste("phi(",pvec,")",sep="")
      covHat<-solve(InformationMatrixARp(phiHat,pvec))/n
      dimnames(covHat)<-list(varNames,varNames)
      sdRacf<-sqrt(diag(VarianceRacfARp(phiHat,pvec,MaxLag,n)))
    }
    else {
      if (PMAX>0) {
        varNames<-paste("phi(",1:PMAX,")",sep="")
        covHat<-SiddiquiMatrix(phiHat)/n
        dimnames(covHat)<-list(varNames,varNames)
        sdRacf<-sqrt(diag(VarianceRacfAR(phiHat,MaxLag,n)))
      }
      else {
        varNames<-character(0)
        covHat<-numeric(0)
        sdRacf<-rep(1/sqrt(n),MaxLag)
      }
    }
    if (SubsetQ) {
      ModelTitle<-deparse(as.numeric(pvec),width.cutoff=180)
      ModelTitle<-paste("ARp",substr(ModelTitle,2,nchar(ModelTitle)),sep="")
      ModelTitle<-gsub(" ", "", ModelTitle)
    }
    else 
      ModelTitle<-paste("AR(",p,")",sep="")
    #
    LBQ<-LjungBoxTest(res, lag.max=MaxLag, k=length(pvec))
    RacfMatrix<-matrix(c(racf,sdRacf),ncol=2)
    dimnames(RacfMatrix)<-list(1:MaxLag, c("ra", "Sd(ra)"))
    zetaHat<-ARToPacf(phiHat)
    #
    if (!MLEQ) { #for LS, report usual LS covariance matrix
      covHat<-ans$covHat
      covHat <- covHat[-1,-1,drop=FALSE] #remove intercept
      dimnames(covHat)<-list(varNames,varNames)
    }
    ans<-list(loglikelihood=ans$loglikelihood,phiHat=phiHat,sigsqHat=sigsq,muHat=mz,covHat=covHat,zetaHat=zetaHat,
              RacfMatrix=RacfMatrix,LjungBoxQ=LBQ,res=res,fits=fits,SubsetQ=SubsetQ,pvec=pvec,FitMethod=FitMethod,
              MeanMLE=MeanMLEQ, ARModel="ARp", tsp=tsp(z),
              call=match.call(),DataTitle=attr(z,"title"),ModelTitle=ModelTitle,z=z)
    class(ans)<-"FitAR"
    ans
  }
