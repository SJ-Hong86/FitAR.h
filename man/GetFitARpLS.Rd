% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetFitARpLS.R
\name{GetFitARpLS}
\alias{GetFitARpLS}
\title{LS for AR(p) and Subset ARp – Short Version}
\usage{
GetFitARpLS(z, pvec)
}
\arguments{
\item{z}{vector or ts object, the time series.}

\item{pvec}{lags included in subset AR. If pvec = 0, white noise assumed.}
}
\value{
a list with components:
\itemize{
\item \code{loglikelihood} the exact loglikelihood.
\item \code{phiHat} estimated AR parameters.
\item \code{constantTerm} constant term in the linear regression.
\item \code{pvec} lags of estimated AR coefficient.
\item \code{res} the least squares regression residuals.
\item \code{InvertibleQ} True, if the estimated parameters are in the AR admissible region.
\item \code{yX} the y vector and X matrix used for the regression fitting.
}
}
\description{
For ARp subset models, the least squares estimates are computed.
The exact loglikelihood is then determined. The estimated parameters
are checked to see if they are in the AR admissible region.
}
\details{
The R function \code{lsfit} is used.
}
\note{
This is a helper function for \code{FitARp} which is invoked by the main
package function \code{FitAR}. Normally the user would \code{FitAR} since this
function provides generic print, summary, resid and plot methods but
\code{GetFitARpLS} is sometimes useful in iterative computations like
bootstrapping since it is faster.
}
\examples{
# Fit subset AR using LS
# normally use FitAR
ans<-FitAR(SeriesA, c(1,2,7), ARModel="ARp", MLEQ=FALSE)
# could also use FitARp
ans<-FitARp(SeriesA, c(1,2,7))
# for some applications GetFitARpLS is simpler and faster
ansLS<-GetFitARpLS(SeriesA, c(1,2,7))
ansLS

}
\references{
McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation
parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
}
\seealso{
\code{\link[=FitAR]{FitAR()}}, \code{\link[=FitARz]{FitARz()}}, \code{\link[=GetFitARz]{GetFitARz()}}, \code{\link[=FitARp]{FitARp()}}, \code{\link[=GetFitARpMLE]{GetFitARpMLE()}}, \code{\link[=RacfPlot]{RacfPlot()}}.
}
\author{
A.I. McLeod.
}
