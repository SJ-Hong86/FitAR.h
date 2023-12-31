% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdfplot.ts.R
\name{sdfplot.ts}
\alias{sdfplot.ts}
\title{Autoregressive Spectral Density Estimation for "ts" Object}
\usage{
sdfplot(obj, ...) ## S3 method for class 'ts'
}
\arguments{
\item{obj}{object, class"ts".}

\item{...}{optional arguments.}
}
\value{
Plot is produced using plot. Matrix with 2 columns containing
the frequencies and spectral density is returned invisibly.
}
\description{
Methods function for "ts".
}
\examples{
data(SeriesA)
sdfplot(SeriesA) 

}
\seealso{
\code{\link[=sdfplot]{sdfplot()}}.
}
\author{
A.I. McLeod.
}
