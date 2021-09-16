#' Sverchkov and Pfeffermann \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Sverchkov and Pfeffermann
#'   (2004, p. 85) bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param pi_sample Numeric vector. Values of inclusion probabilities for sample elements.
#' @param N Integer. Population size.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#' @param scheme Function. Implementation of a sampling scheme (takes a vector of inclusion probabilities as an argument and returns a sample: a vector of 0's and 1's of the same length as argument).
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)` columns. Each row is a bootstrap sample, and the number in i-th column denotes the number of times i-th unit from the original sample occurs in the bootstrap sample.
#'
#' @references Sverchkov M. and Pfeffermann D. (2004). Prediction of Finite Population Totals Based on the Sample Distribution, \emph{Survey Methodology}, 30(1):79-92.
#'
#' @examples
#' boot_SP_pips(1:3, c(0.005, 0.01, 0.015), 100, 10, UPbrewer_rcpp)
#' boot_SP_pips(1:10, (1:10) * 2e-5, 1000, 4, UPbrewer_rcpp)
#'
#' @export
boot_SP_pips <- function(x_sample, pi_sample, N, n_boots, scheme) {
  if (!is.numeric(x_sample) || !all(x_sample > 0)) {
      stop(paste("The first argument of boot_SP function should be a vector of positive numbers, but we got:\n",
                 paste(x_sample, collapse = " ")))
  }
  if (!is.numeric(pi_sample) || !all(pi_sample >= 0 & pi_sample <= 1)) {
      stop(paste("The second argument of boot_SP function should be a vector of probabilitites, but we got:\n",
                 paste(pi_sample, collapse = " ")))
  }
  if (!is.numeric(N) || !(N > 0) ||
      !(N == round(N))) {
    stop(paste("The third argument of boot_SP function should be a single positive integer, but we got:\n",
               paste(N, collapse = " ")))
  }
  if (!is.numeric(n_boots) || !(n_boots > 0) || !(n_boots == round(n_boots))) {
    stop(paste("The fourth argument of boot_SP_pips function should be a single positive integer, but we got:\n",
                paste(n_boots, collapse = " ")))
  }
  scheme <- match.fun(scheme)
  n <- length(x_sample)                           # sample size

  # ni = a number of repetitions in pseudo-population
  # don't try to use "sample", because we need the counts
  # ni <- sampling::UPmultinomial(N / sum(inv_pi_sample) * inv_pi_sample)
  ni <- stats::rmultinom(1, N, 1/pi_sample)[,1]   # that's what sampling::UPmultinomial really is

  x_b    <- rep.int(x_sample, ni)                 # x values in pseudo-population
  pr_b   <- inclusion_probabilities_rcpp(x_b, n)  # inclusion probabilities for pseudo-population
  M      <- matrix(0, nrow = n_boots, ncol = n)   # matrix for results
  groups <- as.factor(rep(1:n, ni))               # indicators of units from the original sample

  # The bootstrap replications are here:
  # - sample from the pseudo-population
  # - sum occurrences of units
  # - put this into matrix row-wise
  M[, ni > 0] <- matrix(replicate(n_boots,
                                  unlist(lapply(split(scheme(pr_b), groups), sum))),
                        nrow = n_boots, byrow = TRUE)
  # return this matrix
  M
}
