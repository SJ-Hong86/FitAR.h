#' Exact MLE for Mean in AR(p)
#' 
#' Details of this algorithm are given in McLeod and Zhang (2007).
#' 
#' @usage GetARMeanMLE(z, phi)
#' @param z vector of length n containing the time series.
#' @param phi vector of AR coefficients.
#' @returns Estimate of mean.
#' @author A.I. McLeod and Y. Zhang.
#' @references McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation 
#'   parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
#' @seealso [mean()].
#' @examples
#' # Simulate a time series with mean zero and compute the exact 
#' # mle for mean and compare with sample average.
#' ## Not run: # save time building package!
#' set.seed(3323)
#' phi<-c(2.7607,-3.8106,2.6535,-0.9238)
#' z<-SimulateGaussianAR(phi,1000)
#' ans1<-mean(z)
#' ans2<-GetARMeanMLE(z,phi)
#' # define a direct MLE function
#' "DirectGetMeanMLE" <-
#' function(z, phi){
#'     GInv<-solve(toeplitz(TacvfAR(phi, length(z)-1)))
#'     g1<-colSums(GInv)
#'     sum(g1*z)/sum(g1)
#' }
#' ans3<-DirectGetMeanMLE(z,phi)
#' ans<-c(ans1,ans2,ans3)
#' names(ans)<-c("mean", "GetARMeanMLE","DirectGetMeanMLE")
#' ans
#' 
#' ## End(Not run)
#' 
#' @export
GetARMeanMLE <-
  function(z, phi){
    stopifnot (length(z)>=2*length(phi))
    g1<-Get1G(phi, length(z))
    sum(g1*z)/sum(g1)
  }
