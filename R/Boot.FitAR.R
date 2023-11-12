#' Simulate a Fitted AR
#' 
#' @description
#' Simulate a realization from a fitted AR model. This is useful in the parametric
#' bootstrap. Generic function for "Boot" method.
#' 
#' @usage Boot(obj, R=1, ...) # S3 method for class 'FitAR'
#' @param obj	 the output from FitAR.
#' @param R	 number of bootstrap replications.
#' @param ...	 optional arguments.
#' @returns A simulated time series with the same length as the original
#'   fitted time series produced when R=1. When R>1, a matrix with R columns
#'   is produced with each column a separate bootstrap realization.
#' @author A.I. McLeod and Y. Zhang.
#' @seealso [Boot()], [SimulateGaussianAR()].
#' @examples 
#' # Plot log(lynx) time series and simulation
#' ans <- FitAR(log(lynx), 8)
#' z<-Boot.FitAR(ans)
#' par(mfrow=c(2,1))
#' TimeSeriesPlot(log(lynx))
#' title(main="log(lynx) time series")
#' TimeSeriesPlot(z)
#' title(main="Simulated AR(8), fitted to log lynx")
#' 
#' # par(mfrow=c(1,1)
#' # Use bootstrap to compute standard errors of parameters
#' # takes about 18 seconds on a 3.6 GHz PC
#' 
#' ## Not run: 
#' ptm <- proc.time() # user time
#' R<-100  # number of bootstrap iterations
#' p<-c(1,2,4,7,10,11)
#' ans<-FitAR(log(lynx),p)
#' out<-Boot(ans, R)
#' fn<-function(z) GetFitARz(z,p)$zetaHat
#' sdBoot<-sqrt(diag(var(t(apply(out,fn,MARGIN=2)))))
#' sdLargeSample<-coef(ans)[,2][1:6]
#' sd<-matrix(c(sdBoot,sdLargeSample),ncol=2)
#' dimnames(sd)<-list(names(sdLargeSample),c("Bootstrap","LargeSample"))
#' ptm<-(proc.time()-ptm)[1]
#' sd
#' ## End(Not run)
#' 
#' @export
Boot.FitAR <-
  function(obj, R=1, ...){
    n<-length(obj$res)
    phiHat<-obj$phiHat
    sigsq<-obj$sigsq
    zm<-obj$muHat
    if (R==1){
      z<-zm+SimulateGaussianAR(phiHat,n,InnovationVariance=sigsq)
      if (is.null(obj$tsp))
        z
      else
        ts(z, start=obj$tsp[1], frequency=obj$tsp[3])
    }
    else {
      x<-matrix(0, nrow=n, ncol=R)
      for (i in 1:R)
        x[,i]<-zm+SimulateGaussianAR(phiHat,n,InnovationVariance=sigsq)
      x
    }
  }
