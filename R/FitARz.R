#' Subset ARz Model Fitting
#' 
#' The subset ARz model, defined by constraining partial autocorrelations to 
#' zero, is fitted using exact MLE. When length(p)=1, an AR(p) is fit by MLE.
#' 
#' @usage FitARz(z, p, demean = TRUE, MeanMLEQ = FALSE, lag.max = "default")
#' @param z time series, vector or ts object.
#' @param p p specifies the model. If length(p) is 1, an AR(p) is assumed and 
#'   if p has length greater than 1, a subset ARz is assumed. For example, 
#'   to fit a subset model with lags 1 and 4 present, set p to c(1,4) or 
#'   equivalently c(1,0,0,4). To fit a subset model with just lag 4, 
#'   you must use p=c(0,0,0,4) since p=4 will fit a full AR(4).
#' @param demean TRUE, mean estimated. FALSE, mean is zero.
#' @param MeanMLEQ use exact MLE for mean parameter.
#' @param lag.max the residual autocorrelations are tabulated for lags 
#'   1, ..., lag.max. Also lag.max is used for the Ljung-Box portmanteau test.
#' @details 
#' The model and its properties are discussed in McLeod and Zhang (2006) 
#' and McLeod and Zhang (2008).
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
#' @note Normally one would use the `FitAR` function which then calls this 
#'   function for the ARz case.
#' @author A.I. McLeod.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [FitAR()], [FitARp()], [GetFitARz()], [GetFitARpMLE()], [RacfPlot()].
#' @examples 
#' # First example: Fit exact MLE to AR(4) 
#' set.seed(3323)
#' phi<-c(2.7607,-3.8106,2.6535,-0.9238)
#' z<-SimulateGaussianAR(phi,1000)
#' ans<-FitARz(z,4,MeanMLEQ=TRUE)
#' ans
#' coef(ans)
#' 
#' ## Not run: # save time building package
#' # Second Example: compare with sample mean result
#' ans<-FitARz(z,4)
#' coef(ans)
#' 
#' # Third Example: Fit subset ARz
#' z<-log(lynx)
#' FitARz(z, c(1,2,4,7,10,11))
#' # now obain exact MLE for Mean as well
#' FitARz(z, c(1,2,4,7,10,11), MeanMLE=TRUE)
#' 
#' # Fourth Example: Fit subset ARz
#' somePACF<-c(0.5,0,0,0,-0.9)
#' someAR<-PacfToAR(somePACF)
#' z<-SimulateGaussianAR(someAR,1000)
#' ans=FitARz(z, c(1,5),MeanMLEQ=TRUE)
#' coef(ans)
#' GetFitARz(z,c(1,5))#assuming a known zero mean
#' 
#' ## End(Not run)
#' 
#' @export
FitARz <-
  function (z, p, demean = TRUE, MeanMLEQ = FALSE, lag.max = "default") 
  {
    stopifnot(length(z) > 0, length(z) > max(p), length(p) > 
                0)
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    stopifnot(is.wholenumber(p), p>0)
    ztsp <- tsp(z)
    if (lag.max == "default") 
      MaxLag <- min(300, ceiling(length(z)/5))
    else    MaxLag = lag.max
    MaxIter <- 10
    n <- length(z)
    pvec <- sort(p)
    pvec <- pvec[pvec > 0]
    if (length(pvec) == 0) 
      pvec <- 0
    if (length(p) == 1 && pvec != 0) 
      pvec <- 1:p
    PMAX <- max(pvec)
    SubQ <- length(pvec) < PMAX
    indMeanQ <- demean || MeanMLEQ
    if (indMeanQ) 
      mz <- mean(z)
    else mz <- 0
    y <- z - mz
    ans <- GetFitARz(y, pvec)
    LL <- ans$loglikelihood
    etol <- 1
    mu <- iter <- 0
    if (MeanMLEQ && PMAX != 0) 
      while (etol > 1e-06 && iter < MaxIter) {
        LLPrev <- LL
        iter <- iter + 1
        mu <- GetARMeanMLE(y, ans$phiHat)
        ans <- GetFitAR(y - mu, pvec)
        LL <- ans$loglikelihood
        etol <- abs(LL - LLPrev)/LLPrev
        if (ans$convergence != 0) 
          stop("GetARFit returned convergence = ", ans$convergence)
      }
    muHat <- mu + mz
    zetaHat <- ans$zetaHat
    phiHat <- ans$phiHat
    if (PMAX != 0) 
      res <- BackcastResidualsAR(y, phiHat, Q = 100, demean = FALSE)
    else res <- y
    fits <- y - res
    sigsq <- sum(res^2)/n
    racf <- (acf(res, plot = FALSE, lag.max = MaxLag)$acf)[-1]
    if (SubQ) {
      varNames <- paste("zeta(", pvec, ")", sep = "")
      covHat <- solve(InformationMatrixARz(zetaHat, pvec))/n
      dimnames(covHat) <- list(varNames, varNames)
      sdRacf <- sqrt(diag(VarianceRacfARz(zetaHat, pvec, MaxLag, 
                                          n)))
    }
    else {
      if (PMAX > 0) {
        varNames <- paste("phi(", 1:PMAX, ")", sep = "")
        covHat <- SiddiquiMatrix(phiHat)/n
        dimnames(covHat) <- list(varNames, varNames)
        sdRacf <- sqrt(diag(VarianceRacfAR(phiHat, MaxLag, 
                                           n)))
      }
      else {
        varNames <- character(0)
        covHat <- numeric(0)
        sdRacf <- rep(1/sqrt(n), MaxLag)
      }
    }
    RacfMatrix <- matrix(c(racf, sdRacf), ncol = 2)
    dimnames(RacfMatrix) <- list(1:MaxLag, c("ra", "Sd(ra)"))
    LBQ <- LjungBoxTest(res, lag.max = MaxLag, k = length(zetaHat))
    if (SubQ) {
      m <- length(pvec)
      if (m < 13) {
        pVEC <- deparse(as.numeric(pvec), width.cutoff = 180)
        pVEC <- substr(pVEC, 2, nchar(pVEC))
      }
      else {
        pVECa <- deparse(as.numeric(pvec[1:4]), width.cutoff = 180)
        pVECa <- substr(pVECa, 2, nchar(pVECa)-1)
        pVECb <- deparse(as.numeric(pvec[(m-2):m]), width.cutoff = 180)
        pVECb <- substr(pVECb, 3, nchar(pVECb)) 
        pVEC <- paste(pVECa, ",...,", pVECb, ", m=",m)
      } 
      ModelTitle <- paste("ARz", pVEC, sep = "")
      ModelTitle <- gsub(" ", "", ModelTitle)
    }
    else ModelTitle <- paste("AR(", p, ")", sep = "")
    ans <- list(loglikelihood = ans$loglikelihood, phiHat = phiHat, 
                sigsqHat = sigsq, muHat = muHat, covHat = covHat, zetaHat = zetaHat, 
                RacfMatrix = RacfMatrix, LjungBoxQ = LBQ, res = res, 
                fits = fits + mz, SubsetQ = SubQ, pvec = pvec, demean = demean, 
                FitMethod = "MLE", iterationCount = iter, convergence = ans$convergence, 
                MeanMLE = MeanMLEQ, tsp = ztsp, call = match.call(), 
                ARModel = "ARz", DataTitle = attr(z, "title"), ModelTitle = ModelTitle, 
                z = z)
    class(ans) <- "FitAR"
    ans
  }
