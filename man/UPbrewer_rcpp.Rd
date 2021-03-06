% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{UPbrewer_rcpp}
\alias{UPbrewer_rcpp}
\title{Brewer sampling}
\usage{
UPbrewer_rcpp(pik)
}
\arguments{
\item{pik}{Numeric vector of  the first order inclusion probabilities.}
}
\value{
Numeric vector of the same length as \code{pik}. Value of 1
indicates that unit is selected into sample, 0 means that unit is
not in the sample.
}
\description{
Uses the Brewer's method to select a sample of units (unequal probabilities,
without replacement, fixed sample size).
We follow the Algorithm 6.10 from Tille (2006), pp. 112-114, however details
of implementation differs a little from that of the \code{sampling} package.
}
\examples{
UPbrewer_rcpp((0:5) / 5)
UPbrewer_rcpp(c(10:0, 5) / 10)

}
\references{
Tille, Y. (2006) Sampling algorithms, Springer, New York.
}
