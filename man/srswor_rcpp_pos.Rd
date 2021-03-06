% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{srswor_rcpp_pos}
\alias{srswor_rcpp_pos}
\title{Simple Random Sampling Without Replacement}
\usage{
srswor_rcpp_pos(N, n)
}
\arguments{
\item{N}{Integer. Population size - number of elements to draw from.}

\item{n}{Integer. Sample size - number of elements to draw.}
}
\value{
Integer vector of length \code{n} of indices of selected units.
Indices have values from the range 0 to N-1 (C++ convention)
}
\description{
Function draws a simple random sample without replacement. This is a fast
replacement for R's \code{sample(N, n, replace = FALSE)}.
}
\examples{
srswor_rcpp_pos(10, 5)

}
