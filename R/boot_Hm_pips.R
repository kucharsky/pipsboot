#' Modified Holmberg \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Holmberg (1998, p. 380)
#' bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling modified such that
#' the generation of the pseudo-population occurs in each bootstrap iteration.
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the population.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#' @param scheme Function. Implementation of a sampling scheme (takes a vector
#'   of inclusion probabilities as an argument and returns a sample: a vector of
#'   0's and 1's of the same length as argument).
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#'   columns. Each row is a bootstrap sample, and the number in i-th column
#'   denotes the number of times i-th unit from the original sample occurs in
#'   the bootstrap sample.
#'
#' @references Holmberg, A.(1998). A bootstrap approach to probability proportional
#' to size sampling. \emph{Proceedings of Section on Survey Research Methods}.
#' American Statistical Association, Washington, 378-383.
#'
#' @examples
#' boot_Hm_pips(1:3, 100, 10, UPbrewer_rcpp)
#' boot_Hm_pips(1:10, 1000, 4, UPbrewer_rcpp)
#'
#' @export
boot_Hm_pips <- function(x_sample, x_pop_total, n_boots, scheme) {
  if (!is.numeric(x_sample) || !all(x_sample > 0)) {
    stop(paste("The first argument of boot_HBM_pips function should be a vector of positive numbers, but we got:\n",
               paste(x_sample, collapse = " ")))
  }
  if (!is.numeric(x_pop_total) || !length(x_pop_total == 1)) {
    stop(paste("The second argument of boot_HBM_pips function should be a single number, but we got:\n",
               paste(x_pop_total, collapse = " ")))
  }
  if (!is.numeric(n_boots) || !(n_boots > 0) || !(n_boots == round(n_boots))) {
    stop(paste("The third argument of boot_HBM_pips function should be a single positive integer, but we got:\n",
               paste(n_boots, collapse = " ")))
  }
  scheme <- match.fun(scheme)
  n <- length(x_sample)                        # sample size
  inv_prs <- x_pop_total / (x_sample * n)      # a number of repetitions in pseudo-population
  f_inv_prs <- floor(inv_prs)
  rk <- inv_prs - f_inv_prs
  M <- matrix(0, nrow = n_boots, ncol = n)     # matrix for results
  # Loop below replaces
  # matrix(replicate(n_boots, boot_H_pips(x_sample, x_pop_total, 1, scheme)),
  #        nrow = n_boots, byrow = TRUE)
  # to avoid some overhead caused by calling boot_H_pips multiple times
  for(i in 1:n_boots) {
    ni <- f_inv_prs + stats::rbinom(n, 1, rk)
    x_b <- rep.int(x_sample, ni)                        # x values in pseudo-population
    pr_b <- inclusion_probabilities_rcpp(x_b, n)        # inclusion probabilities for pseudo-population
    groups <- as.factor(rep.int(1:n, ni))               # indicators of units from the original sample
    M[i, ni > 0] <- unlist(lapply(split(scheme(pr_b), groups), sum))
  }
  M
}
