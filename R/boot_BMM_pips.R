#' Barbiero, Manzi and Mecatti \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Barbiero, Manzi and Mecatti
#'   (2015, pp. 610-611) bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS
#'   sampling via calibrated pseudo-population.
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the population.
#' @param g_weights Numeric vector. Values of weights for sample elements.
#' Products of `g_weights` and inverses of first order inclusion probabilities
#' are numbers of repetitions of sample elements in the pseudo-population.
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
#' @references Barbiero A., Manzi G. and Mecatti F. (2015) Bootstrapping
#' probability-proportional-to-size samples via calibrated empirical population.
#' Journal of Statistical Computation and Simulation, 85(3), 608-620.
#'
#' @examples
#' N <- 200   # population size
#' x_population <- rgamma(N, 1)
#' n <- 10    # sample size
#' pi <- inclusion_probabilities_rcpp(x_population, n)
#' set.seed(123)
#' sample <- UPbrewer_rcpp(pi)
#' x_s <- x_population[sample == 1]
#' pi_s <- pi[sample == 1]
#' # g-weigths
#' g <- sampling::calib(cbind(x_s, rep(1, n)),
#'                      d = 1 / pi_s,
#'                      c(sum(x_population), N),
#'                      method = "truncated",
#'                      bounds = c(low = 0, upp = 10))
#' boot_BMM_pips(x_s, sum(x_population), g, 5, UPbrewer_rcpp)
#'
#' @export
boot_BMM_pips <- function(x_sample, x_pop_total, g_weights, n_boots, scheme) {
  if (!is.numeric(x_sample) || !all(x_sample > 0)) {
    stop(paste("The first argument of boot_BMM_pips function should be a vector of positive numbers, but we got:\n",
               paste(x_sample, collapse = " ")))
  }
  if (!is.numeric(x_pop_total) || !length(x_pop_total == 1)) {
    stop(paste("The second argument of boot_BMM_pips function should be a single number, but we got:\n",
               paste(x_pop_total, collapse = " ")))
  }
  if (!is.numeric(g_weights)) {
    stop(paste("The third argument of boot_BMM_pips function should be a numeric vector, but we got:\n",
                 paste(g_weights, collapse = " ")))
  }
  if (!is.numeric(n_boots) || !(n_boots > 0) || !(n_boots == round(n_boots))) {
    stop(paste("The fourth argument of boot_BMM_pips function should be a single positive integer, but we got:\n",
               paste(n_boots, collapse = " ")))
  }
  scheme <- match.fun(scheme)
  n <- length(x_sample)                             # sample size
  ni_A <- g_weights * x_pop_total / (x_sample * n)  # a number of repetitions in pseudo-population
  ni <- round(ni_A) * as.numeric(ni_A >= 0)         # to remove negative values
  x_b <- rep.int(x_sample, ni)                      # x values in pseudo-population
  pr_b <- inclusion_probabilities_rcpp(x_b, n)      # inclusion probabilities for pseudo-population
  M <- matrix(0, nrow = n_boots, ncol = n)          # matrix for results
  groups <- as.factor(rep.int(1:n, ni))             # indicators of units from the original sample
  # The bootstrap replications are here:
  # - sample from the pseudo-population
  # - sum occurrences of units
  # - put this into matrix row-wise
  M[, ni > 0] <- matrix(replicate(n_boots,
                                  unlist(lapply(split(scheme(pr_b), groups), sum))
                                 ),
                        nrow = n_boots, byrow = TRUE)
  # return this matrix
  M
}
