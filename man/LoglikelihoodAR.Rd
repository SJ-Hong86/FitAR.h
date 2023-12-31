% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoglikelihoodAR.R
\name{LoglikelihoodAR}
\alias{LoglikelihoodAR}
\title{Exact Loglikelihood for AR}
\usage{
LoglikelihoodAR(phi, z, MeanValue = 0)
}
\arguments{
\item{phi}{AR parameters.}

\item{z}{time series data, not assumed mean corrected.}

\item{MeanValue}{usually this is mean(z) but it could be another value for
example the MLE of the mean.}
}
\value{
The value of the loglikelihood is returned.
}
\description{
The exact loglikelihood function, defined in eqn. (6) of McLeod &
Zhang (2006) is computed. Requires O(n) flops, n = length(z).
}
\details{
Eqn (6) of McLeod and Zhang (2006) may be written
\deqn{-(n/2)\mathrm{log}(\hat{\sigma}_a^2)-(1/2)\mathrm{log}(g_p)}
where \eqn{\hat{\sigma}_a^2} is the residual variance and \eqn{g_p} is the
covariance determinant.
}
\note{
No check is done for stationary-causal process.

For MLE computation it is better to use FastLoglikelihoodAR since for repeated
likelihood evaluations this requires only O(1) flops vs O(n) flops, where n =
length(z).
}
\examples{
# Fit a subset model to Series A and verify the loglikelihood
out<-FitAR(SeriesA, c(1,2,7))
out
# either using print.default(out) to see the components in out
# or applying LoglikelihoodAR () by first obtaining the phi parameters as out$phiHat.
LoglikelihoodAR(out$phiHat, SeriesA, MeanValue=mean(SeriesA))

}
\references{
McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation
parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
}
\seealso{
\code{\link[=FastLoglikelihoodAR]{FastLoglikelihoodAR()}}.
}
\author{
A.I. McLeod and Y. Zhang.
}
