% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InformationMatrixARp.R
\name{InformationMatrixARp}
\alias{InformationMatrixARp}
\title{Fisher Information Matrix Subset Case, ARp}
\usage{
InformationMatrixARp(phi, lags)
}
\arguments{
\item{phi}{vector of coefficients in the subset AR.}

\item{lags}{vector indicating lags present in phi.}
}
\value{
a p-by-p Toeplitz matrix, p = length(phi).
}
\description{
The large-sample information matrix per observation is computed in a subset
AR with the usual parameterization, that is, a subset of the AR coefficients.
}
\details{
The subset information matrix is obtained simply by selecting the appropriate
rows and columns from the full information matrix. This function is used by
\code{FitARp} to obtain the estimated standard errors of the parameter estimates.
}
\examples{
# variances of parameters in a subset ARp(1,2,6)
fi<-InformationMatrixARp(c(0.36,0.23,0.23),c(1,2,6))
sqrt(diag(solve(fi*197)))

}
\references{
McLeod, A.I. and Zhang, Y. (2006). Partial autocorrelation
parameterization for subset autoregression. Journal of Time Series Analysis, 27, 599-612.
}
\seealso{
\code{\link[=InformationMatrixARp]{InformationMatrixARp()}}, \code{\link[=FitARp]{FitARp()}}, \code{\link[=InformationMatrixARz]{InformationMatrixARz()}}.
}
\author{
A.I. McLeod and Y. Zhang.
}
