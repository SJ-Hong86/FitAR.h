#' Jarque-Bera Normality Test
#' 
#' A powerful omnibus test for normality.
#' 
#' @usage JarqueBeraTest(z)
#' @param z vector of data.
#' @details This test is derived as a Lagrange multiplier test for normal 
#'   distribution in the family of Pearson distributions (Jarque and Bera, 1987).
#' @returns
#' The returned values are as follow:
#' 
#' * `LM` value of the LM statistic.
#' * `pvalue` p-value.
#' @author A.I. McLeod.
#' @references Jarque, C.M. and Bera, A.K. (1987). A Test for Normality of 
#'   Observations and Regression Residuals. International Statistical Review 55, 163-172.
#' @examples
#' # some normal data
#' z<-rnorm(100)
#' JarqueBeraTest(z)
#' # some skewed data
#' z<-rexp(100)
#' JarqueBeraTest(z)
#' # some thick tailed data
#' z<-rt(100,5)
#' JarqueBeraTest(z)
#' 
#' @export
JarqueBeraTest <-
  function(z){
    x<-z-mean(z)
    n<-length(x)
    m2<-sum(x^2)/n
    m3<-sum(x^3)/n
    m4<-sum(x^4)/n
    g3<-m3/(m2^(3/2))
    g4<-m4/(m2^2)
    LM<-n*((g3^2)/6+((g4-3)^2)/24)
    pv<-1-pchisq(LM,2)
    list(LM=LM,pvalue=pv)
  }
