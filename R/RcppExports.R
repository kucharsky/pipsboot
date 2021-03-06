# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Brewer sampling
#'
#' Uses the Brewer's method to select a sample of units (unequal probabilities,
#' without replacement, fixed sample size).
#' We follow the Algorithm 6.10 from Tille (2006), pp. 112-114, however details
#' of implementation differs a little from that of the \code{sampling} package.
#'
#' @param pik Numeric vector of  the first order inclusion probabilities.
#'
#' @return Numeric vector of the same length as \code{pik}. Value of 1
#'         indicates that unit is selected into sample, 0 means that unit is
#'         not in the sample.
#'
#' @references Tille, Y. (2006) Sampling algorithms, Springer, New York.
#'
#' @examples
#' UPbrewer_rcpp((0:5) / 5)
#' UPbrewer_rcpp(c(10:0, 5) / 10)
#'
#' @export
UPbrewer_rcpp <- function(pik) {
    .Call(`_pipsboot_UPbrewer_rcpp`, pik)
}

#' Poisson sampling
#'
#' Draws a Poisson sample (without replacement, random sample size) using a
#' given vector of first-order inclusion probabilities. Rewriten in Rcpp based
#' on the `\link[sampling]{UPpoisson}` function.
#'
#' @param pik  Numeric vector. First-order inclusion probabilities.
#'
#' @return Integer vector of size \code{length(pik)} (population size) with
#' elements 1 (if the unit is selected in the sample) or 0 (if it is not).
#'
#' @examples
#' UPpoisson_rcpp((1:10) / 55 * 5)
#' UPpoisson_rcpp(inclusion_probabilities_rcpp(runif(10), 5))
#'
#' @export
UPpoisson_rcpp <- function(pik) {
    .Call(`_pipsboot_UPpoisson_rcpp`, pik)
}

#' Barbiero and Mecatti 0.5 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Barbiero and Mecatti
#'   (2010, p. 62) "0.5 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS algorithm".
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for
#'   sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the
#'   population.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#' @param scheme Function. Implementation of a sampling scheme (takes a vector
#'   of inclusion probabilities as an argument and returns a sample: a vector of
#'   0's and 1's of the same length as argument).
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#'   columns. Each row is a bootstrap sample, and the number in i-th column
#'   denotes the number of times i-th unit from the original sample occurs in
#'   the bootstrap sample.
#' @export
#' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
#'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
#'   in: Complex Data Modeling and Computationally Intensive Statistical
#'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
#' @examples
#' boot_05_pips_rcpp(1:3, 100, 10, UPbrewer_rcpp)
#' boot_05_pips_rcpp(1:10, 1000, 5, UPbrewer_rcpp)
boot_05_pips_rcpp <- function(x_sample, x_pop_total, n_boots, scheme) {
    .Call(`_pipsboot_boot_05_pips_rcpp`, x_sample, x_pop_total, n_boots, scheme)
}

#' Barbiero and Mecatti 0.5 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Barbiero and Mecatti
#'   (2010, p. 62) "0.5 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS algorithm".
#'   Uses Brewers' design for sampling.
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for
#'   sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the
#'   population.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#'   columns. Each row is a bootstrap sample, and the number in i-th column
#'   denotes the number of times i-th unit from the original sample occurs in
#'   the bootstrap sample.
#' @export
#' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
#'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
#'   in: Complex Data Modeling and Computationally Intensive Statistical
#'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
#' @examples
#' boot_05_pips_rcpp(1:3, 100, 10, UPbrewer_rcpp)
#' boot_05_pips_rcpp(1:10, 1000, 5, UPbrewer_rcpp)
boot_05_pips_rcpp_brewer <- function(x_sample, x_pop_total, n_boots) {
    .Call(`_pipsboot_boot_05_pips_rcpp_brewer`, x_sample, x_pop_total, n_boots)
}

#' Antal and Tille (2011) bootstrap for simple random sampling without replacement
#'
#' This function implements method of resampling described in the algorithm 3 of
#' Antal & Tille (2011), p. 539.
#'
#' @param N Integer. Original population size.
#' @param n Integer. Sample size.
#' @param n_boots Integer. Number of bootstrap samples to be generated.
#'
#' @return Integer matrix with \code{n_boots} rows and \code{n} columns.
#' Each row is a bootstrap sample, with entries in columns denoting how many
#' times a given unit is resampled.
#'
#' @references Erika Antal and Yves Tille (2011) A Direct Bootstrap Method for
#' Complex Sampling Designs From a Finite Population, Journal of the American
#' Statistical Association, 106:494, 534-543, DOI: 10.1198/jasa.2011.tm09767
#'
#' @examples
#' boot_AT2011_alg3_rcpp(10, 5, 20)
#'
#' @export
boot_AT2011_alg3_rcpp <- function(N, n, n_boots) {
    .Call(`_pipsboot_boot_AT2011_alg3_rcpp`, N, n, n_boots)
}

#' Antal and Tille (2011) bootstrap for complex sampling designs
#'
#' Function implements a method of resampling from a sample selected with unequal
#' probabilities without replacement as described in algorithms 4 and 5 in
#' Antal & Tille (2011), pp. 539-540.
#'
#' @param pi_sample Numeric vector. Inclusion probabilities for elements in the original sample.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#' @param scheme Function. Implementation of a sampling scheme (takes a vector
#'   of inclusion probabilities as an argument and returns a sample: a vector of
#'   0's and 1's of the same length as the argument) used to select the original sample.
#'
#' @return Integer Matrix with \code{n_boots} rows and \code{length(pi_sample)}
#' columns. Each row is a bootstrap sample with entries denoting how many
#' times given unit is resampled.
#'
#' @references Erika Antal and Yves Tille (2011) A Direct Bootstrap Method for
#' Complex Sampling Designs From a Finite Population, Journal of the American
#' Statistical Association, 106:494, 534-543, DOI: 10.1198/jasa.2011.tm09767
#'
#' @examples
#' pi_sample <- inclusion_probabilities_rcpp(1:10, 6)
#' boot_AT2011_alg45_rcpp(pi_sample, 5, UPbrewer_rcpp)
#'
#' @export
boot_AT2011_alg45_rcpp <- function(pi_sample, n_boots, scheme) {
    .Call(`_pipsboot_boot_AT2011_alg45_rcpp`, pi_sample, n_boots, scheme)
}

#' Antal and Tille (2014) bootstrap for complex sampling designs
#'
#' This function implements resampling algorithm described in Antal and Tille
#' (2014), pp. 1355-1356 for sampling with unequal probabilities without replacement.
#'
#' @param pi_sample Numeric vector of inclusion probabilities of the first kind
#' (of elements in the original sample).
#' @param n_boots Integer. Number of bootstrap samples to create.
#'
#' @return Integer matrix with `n_boots` rows and `length(phik)` columns.
#' Each row is a bootstrap sample, with entries in columns denoting how many
#' times given unit is resampled.
#'
#' @references Erika Antal and Yves Tille (2014) A new resampling method for
#' sampling designs without replacement: the doubled half bootstrap.
#' Comput Stat 29, 1345-1363, https://doi.org/10.1007/s00180-014-0495-0.
#'
#' @export
boot_AT2014_rcpp <- function(pi_sample, n_boots) {
    .Call(`_pipsboot_boot_AT2014_rcpp`, pi_sample, n_boots)
}

#' Quatember \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Quatember (2014, p. 95)
#' bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling.
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param x_pop_total Numeric. Sum of the auxiliary variable in the population.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#'   columns. Each row is a bootstrap sample, and the number in i-th column
#'   denotes how many times a given unit is resampled.
#'
#' @references Quatember, A. (2014). The Finite Population Bootstrap - From the
#' Maximum Likelihood to the Horvitz-Thompson Approach. Austrian Journal of
#' Statistics, 43(2), 93-102. https://doi.org/10.17713/ajs.v43i2.10
#'
#' @examples
#' boot_Q_pips_rcpp(1:3, 100, 10)
#' boot_Q_pips_rcpp(1:10, 1000, 4)
#'
#' @export
boot_Q_pips_rcpp <- function(x_sample, x_pop_total, n_boots) {
    .Call(`_pipsboot_boot_Q_pips_rcpp`, x_sample, x_pop_total, n_boots)
}

#' Generalized Quatember \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using generalized Quatember (2014, p. 95)
#' bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling where additional weights
#' are taken into account to compute the number of replicates in the pseudo-population.
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the population.
#' @param g_weights Numeric vector of the length `length(x_sample)`. Positive
#' weights for sample elements. Should sum up to the sample size to ensure that
#' products of 'x_sample' and calibration weights (defined as products of
#' 'g_weights` and inverses of first order inclusion probabilities) are equal
#' to 'x_pop_total'. If this condition is not met, they are rescaled.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#'   columns. Each row is a bootstrap sample, with entries in columns denoting
#'   how many times given unit is resampled.
#'
#' @references Quatember, A. (2014). The Finite Population Bootstrap - From the
#' Maximum Likelihood to the Horvitz-Thompson Approach. Austrian Journal of
#' Statistics, 43(2), 93-102. https://doi.org/10.17713/ajs.v43i2.10
#'
#' @references Zadlo T. (2021) On the generalisation of Quatember's
#' bootstrap, Statistics in Transition new series, 22(1), 163-178,
#' https://doi.org/10.21307/stattrans-2021-009.
#'
#' @examples
#' N <- 200 # population size
#' x_population <- rgamma(N, 1)
#' n <- 10 # sample size
#' pi <- inclusion_probabilities_rcpp(x_population, n)
#' set.seed(123)
#' sample <- UPbrewer_rcpp(pi)
#' x_s <- x_population[sample == 1]
#' pi_s <- pi[sample == 1]
#' # g-weigths
#' g <- sampling::calib(Xs = cbind(x_s, rep(1, n)),
#'                      d = 1 / pi_s,
#'                      total = c(sum(x_population), N),
#'                      method = "truncated",
#'                      bounds = c(low = 0, upp = 10))
#' boot_Qg_pips_rcpp(x_s, sum(x_population), g, 5)
#'
#' @export
boot_Qg_pips_rcpp <- function(x_sample, x_pop_total, g_weights, n_boots) {
    .Call(`_pipsboot_boot_Qg_pips_rcpp`, x_sample, x_pop_total, g_weights, n_boots)
}

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
#' boot_SP_pips_rcpp(1:3, c(0.005, 0.01, 0.015), 100, 10, UPbrewer_rcpp)
#' boot_SP_pips_rcpp(1:10, (1:10) * 2e-5, 1000, 4, UPbrewer_rcpp)
#'
#' @export
boot_SP_pips_rcpp <- function(x_sample, pi_sample, N, n_boots, scheme) {
    .Call(`_pipsboot_boot_SP_pips_rcpp`, x_sample, pi_sample, N, n_boots, scheme)
}

#' Zadlo \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples mimicking Brewer's \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS
#'   sampling scheme (Tille (2006), 112-114) for pseudo-population where numbers
#'   of replications of sample elements do not have to be integers and the
#'   pseudo-population is not physically constructed.
#'
#' @param n_replications Numeric vector. Positive values (not necessarily
#' integers) of number of replications of sample elements in the pseudo-populations
#' e.g. inverses of first order inclusion probabilities for sample elements or
#' calibration weights for sample elements.
#' @param pi_sample Numeric vector. Values of first order inclusion probabilities for sample elements.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)` columns.
#' Each row is a bootstrap sample, and the number in i-th column denotes how many
#' times a given unit is resampled.
#'
#' @references Tille, Y. (2006) Sampling algorithms, Springer, New York.
#'
#' @examples
#' N <- 200 # population size
#' x_population <- rgamma(N, 1)
#' n <- 10  # sample size
#' pi <- inclusion_probabilities_rcpp(x_population, n)
#' set.seed(123)
#' sample <- UPbrewer_rcpp(pi)
#' x_s <- x_population[sample == 1]
#' pi_s <- pi[sample == 1]
#' # g-weigths
#' g <- sampling::calib(matrix(c(x_s, rep(1, n)), ncol = 2),
#'                      d = 1 / pi_s, c(sum(x_population), N), method = "linear")
#' # numbers of replications of sample elements in the pseudo-population
#' # equal inverses of first order inclusion probabilities:
#' boot_Z_pips_rcpp(1 / pi_s, pi_s, 5)
#' # numbers of replications of sample elements in the pseudo-population equal to calibration weights:
#' boot_Z_pips_rcpp(g / pi_s, pi_s, 5)
#'
#' @export
boot_Z_pips_rcpp <- function(n_replications, pi_sample, n_boots) {
    .Call(`_pipsboot_boot_Z_pips_rcpp`, n_replications, pi_sample, n_boots)
}

#' Barbiero and Mecatti x-balanced 1 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Barbiero and Mecatti
#'   (2010, p. 62-63) x-balanced 1 bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the population.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#' @param scheme Function. Implementation of a sampling scheme (takes a vector
#'   of inclusion probabilities as an argument and returns a sample: a vector of
#'   0's and 1's of the same length as argument).
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#' columns. Each row is a bootstrap sample, and the number in i-th column
#' denotes how many times a given unit is resampled.
#'
#' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
#'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
#'   in: Complex Data Modeling and Computationally Intensive Statistical
#'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
#'
#' @examples
#' boot_xbal1_pips_rcpp(1:3, 100, 10, UPbrewer_rcpp)
#' boot_xbal1_pips_rcpp(1:10, 1000, 4, UPbrewer_rcpp)
#'
#' @export
boot_xbal1_pips_rcpp <- function(x_sample, x_pop_total, n_boots, scheme) {
    .Call(`_pipsboot_boot_xbal1_pips_rcpp`, x_sample, x_pop_total, n_boots, scheme)
}

#' Barbiero and Mecatti x-balanced 1 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Barbiero and Mecatti
#'   (2010, p. 62-63) x-balanced 1 bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the population.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#' columns. Each row is a bootstrap sample, and the number in i-th column
#' denotes how many times a given unit is resampled.
#'
#' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
#'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
#'   in: Complex Data Modeling and Computationally Intensive Statistical
#'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
#'
#' @examples
#' boot_xbal1_pips_rcpp_brewer(1:3, 100, 10)
#' boot_xbal1_pips_rcpp_brewer(1:10, 1000, 4)
#'
#' @export
boot_xbal1_pips_rcpp_brewer <- function(x_sample, x_pop_total, n_boots) {
    .Call(`_pipsboot_boot_xbal1_pips_rcpp_brewer`, x_sample, x_pop_total, n_boots)
}

#' Barbiero and Mecatti x-balanced 2 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Barbiero and Mecatti
#'   (2010, p. 63) x-balanced 2 bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling
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
#'   denotes how many times a given unit is resampled.
#'
#' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
#'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
#'   in: Complex Data Modeling and Computationally Intensive Statistical
#'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
#'
#' @examples
#' boot_xbal2_pips_rcpp( 1:3, 100, 10, UPbrewer_rcpp)
#' boot_xbal2_pips_rcpp(1:10, 1000, 4, UPbrewer_rcpp)
#'
#' @export
boot_xbal2_pips_rcpp <- function(x_sample, x_pop_total, n_boots, scheme) {
    .Call(`_pipsboot_boot_xbal2_pips_rcpp`, x_sample, x_pop_total, n_boots, scheme)
}

#' Barbiero and Mecatti x-balanced 2 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
#'
#' The function generates bootstrap samples using Barbiero and Mecatti
#'   (2010, p. 63) x-balanced 2 bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling
#'
#' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
#' @param x_pop_total Number. Sum of the auxiliary variable in the population.
#' @param n_boots Integer. A number of bootstrap samples to be generated.
#'
#' @return Integer matrix with `n_boots` rows and `length(x_sample)`
#'   columns. Each row is a bootstrap sample, and the number in i-th column
#'   denotes how many times a given unit is resampled.
#'
#' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
#'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
#'   in: Complex Data Modeling and Computationally Intensive Statistical
#'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
#'
#' @references Tille, Y. (2006) Sampling algorithms, Springer, New York.
#'
#' @examples
#' boot_xbal2_pips_rcpp_brewer(1:3, 100, 10)
#' boot_xbal2_pips_rcpp_brewer(1:10, 1000, 4)
#'
#' @export
boot_xbal2_pips_rcpp_brewer <- function(x_sample, x_pop_total, n_boots) {
    .Call(`_pipsboot_boot_xbal2_pips_rcpp_brewer`, x_sample, x_pop_total, n_boots)
}

#' Doubled half sampling
#'
#' Draws a sample using "doubled half sampling" described in
#' Antal and Tille (2014), pp. 1351-1352. Samples are drawn with replacement.
#'
#' @param n Integer. Number of elements to sample from.
#'
#' @return Integer vector of size \code{n}. Values indicates how many times
#'         units are selected into sample, 0 means that the unit is not
#'         in the sample.
#'
#' @references Antal, E. and Tille, Y. (2014) A new resampling method for
#' sampling designs without replacement: the doubled half bootstrap.
#' Comput Stat 29, 1345-1363, https://doi.org/10.1007/s00180-014-0495-0
#'
#' @examples
#' halfdoubled_rcpp(5)
#' halfdoubled_rcpp(20)
#'
#' @export
halfdoubled_rcpp <- function(n) {
    .Call(`_pipsboot_halfdoubled_rcpp`, n)
}

#' Inclusion probabilities
#'
#' Calculates first order inclusion probabilities for a probability
#' proportional-to-size sampling design based on positive values of an
#' auxiliary variable. We use the algorithm described
#' in Tille (2006), pp. 18-19.
#'
#' @param values Numeric vector. Positive values of an auxiliary variable.
#' @param size   Integer. Sample size.
#'
#' @return Vector of size \code{length(values)} of inclusion probabilities.
#'
#' @references Tille, Y. (2006) Sampling algorithms, Springer, New York.
#'
#' @examples
#' inclusion_probabilities_rcpp(1:10, 5)
#' inclusion_probabilities_rcpp(1:10, 6)
#' inclusion_probabilities_rcpp(c(1:8, 50, 100), 5)
#'
#' @export
inclusion_probabilities_rcpp <- function(values, size) {
    .Call(`_pipsboot_inclusion_probabilities_rcpp`, values, size)
}

#' Inverse Hypergeometric random numbers
#'
#' Generates random number from a special case of the inverse hypergeometric distribution.
#'
#' Generates random number from a discrete distribution with the probability mass function:
#' \deqn{\mathbb{P}(X = j) = \dfrac{\dbinom{N + n - j - 2}{n - j}}{\dbinom{N + n - 1}{n}} \qquad j = 0,\dots,n,}
#' being a special case of the inverse hypergeometric distribution.
#'
#' This distribution may be described as follows: we have `n` red balls and `N-1` black balls.
#' We draw balls without replacement, and `X` is the number of red balls drawn before first black ball.
#'
#' The expectation of this distribution is \eqn{\mathbb{E}(X) = n/N}
#' and variance \eqn{\mathrm{Var}(X) = \dfrac{n(N-1)(N+n)}{N^2(N+1)}}.
#' Corner cases are covered in the following way:
#' if \eqn{n = 0} then \eqn{X \equiv 0} and if \eqn{N = 1} then \eqn{X \equiv n}.
#'
#' @param N Positive integer
#' @param n Vector of nonnegative integers.
#'
#' @return Integer.
#'
#' @examples
#' inverse_hypergeometric_rcpp_scalar(10, 20)
#'
#' @export
inverse_hypergeometric_rcpp_scalar <- function(N, n) {
    .Call(`_pipsboot_inverse_hypergeometric_rcpp_scalar`, N, n)
}

#' Bootstrap replications of sample data
#'
#' The function computes a matrix with bootstrap replications of sample data.
#'
#' @param bootstrap_samples Integer matrix with `B` rows and `n` columns,
#'   where `B` is the number of bootstrap replications and `n` is the sample
#'   size. Each row is a bootstrap sample, with entries in columns denoting how
#'   many times given unit is resampled.
#' @param x_sample Numeric vector. Values of the variable for sampled elements.
#'
#' @return Matrix with `B` rows and `n` columns. Each row is a bootstrap sample of values of `x_sample`.
#'
#' @examples
#' x_sample <- 1:10
#' set.seed(111)
#' bootstrap_samples <- boot_H_pips(x_sample, 1000, 4, UPbrewer_rcpp)
#' sampledata_bootrep_rcpp(bootstrap_samples, x_sample)
#'
#' @export
sampledata_bootrep_rcpp <- function(bootstrap_samples, x_sample) {
    .Call(`_pipsboot_sampledata_bootrep_rcpp`, bootstrap_samples, x_sample)
}

#' One-one resampling design
#'
#' Implements sampling design described in
#' Antal and Tille (2011, page 537, Algorithm 2).
#'
#' @param n Integer. Sample size.
#'
#' @return Integer vector of size `n. The value of `i`-th element tells how
#' many times `i`-th unit appears in the sample.
#'
#' @references Erika Antal and Yves Tille (2011),
#' A Direct Bootstrap Method for Complex Sampling Designs From a Finite
#' Population, Journal of the American Statistical Association 106:494,
#' 534-543, DOI: 10.1198/jasa.2011.tm09767
#'
#' @details
#' Fixed sample size resampling design with replacement. Citing
#' Antal and Tille (2011) p. 537: "an ad hoc construction used to randomly
#' select `n` units from a sample of size `n` in such a way that the
#' expectation and the variance of \eqn{S^*_k} are equal to 1 (...). This
#' sampling design is a mixture between a simple random sampling with
#' replacement and a simple random sampling with  over-replacement".
#'
#' @examples
#' sampling_oneone_rcpp(10)
#'
#' @export
sampling_oneone_rcpp <- function(n) {
    .Call(`_pipsboot_sampling_oneone_rcpp`, n)
}

#' Sampling with over-replacement
#'
#' Implements sampling design described in Antal and Tille (2011),
#' and called there "Simple random sampling with over-replacement".
#'
#' @param N Integer. Size of the population.
#' @param n Integer. Number of elements in the sample.
#'
#' @return Integer vector of length `N`.
#' The value of `i`-th coordinate tells how many times `i`-th unit appears
#' in the sample.
#'
#' @references E. Antal and Y. Tille (2011)
#' Simple random sampling with over-replacement,
#' Journal of Statistical Planning and Inference 141, 597-601.
#'
#' @details
#' A sample of fixed size `n` is drawn with replacement from population of
#' size `N`. All samples have the same probability of being selected:
#' \eqn{\left(\dbinom{N + n - 1}{n}\right)^{-1}}.
#'
#' After Antal and Tille (2011) we describe the derivation of this design.
#' Consider a sequence of `N` independent geometric random variables
#' \eqn{X_1, \dots, X_N}. Conditioning on \eqn{X_1 + \dots + X_N = n}, we
#' obtain a distribution such that
#' \deqn{P(X_1 = x_1, \dots, X_N = x_n) = 1 / \dbinom{N+n-1}{n}} for every
#' vector of nonnegative numbers
#' \eqn{(x_1,\dots, x_N)} with \eqn{x_1 + \dots + x_N = n}. Marginals of this
#' distribution are inverse hypergeometric. We have
#' \eqn{\mathbb{E}(X_k) = n/N},
#' \eqn{\mathrm{var}(X_k) = \dfrac{n(N-1)(N+n)}{N^2(N+1)}} and
#' \eqn{\mathrm{cov}(X_k, X_l) = -\mathrm{var}(X_k)/(N-1)} if \eqn{k\not=l}.
#'
#' @examples
#' sampling_with_overreplacement_rcpp_scalar(10, 8)
#'
#' @export
sampling_with_overreplacement_rcpp_scalar <- function(N, n) {
    .Call(`_pipsboot_sampling_with_overreplacement_rcpp_scalar`, N, n)
}

#' Simple Random Sampling Without Replacement
#'
#' Function draws a simple random sample without replacement. This is an Rcpp
#' replacement for \code{\link[sampling]{srswor}}, however with different
#' implementation.
#'
#' @param N Integer. Population size - number of elements to draw from.
#' @param n Integer. Sample size - number of elements to draw.
#'
#' @return Integer vector of size \code{N} with exactly \code{n} elements
#' equal to 1 (units selected into the sample) and others equal to 0 (units not
#' selected).
#'
#' @examples
#' srswor_rcpp_01(10, 5)
#'
#' @export
srswor_rcpp_01 <- function(N, n) {
    .Call(`_pipsboot_srswor_rcpp_01`, N, n)
}

#' Simple Random Sampling Without Replacement
#'
#' Function draws a simple random sample without replacement. This is a fast
#' replacement for R's \code{sample(N, n, replace = FALSE)}.
#'
#' @param N Integer. Population size - number of elements to draw from.
#' @param n Integer. Sample size - number of elements to draw.
#'
#' @return Integer vector of length \code{n} of indices of selected units.
#' Indices have values from the range 0 to N-1 (C++ convention)
#'
#' @examples
#' srswor_rcpp_pos(10, 5)
#'
#' @export
srswor_rcpp_pos <- function(N, n) {
    .Call(`_pipsboot_srswor_rcpp_pos`, N, n)
}

#' Simple Random Sampling With Replacement
#'
#' Samples `n` integers from the range `0` to `N-1`.
#' A call `srswr_rcpp_zero(N, n)` is statistically equivalent to R's
#' `sample(N, n, replace = TRUE) - 1` however the results of every individual
#' call are different due to differences in the implementation.
#' `srswr_rcpp_zero` is up to 5 times faster.
#'
#' @param N Integer. Number of elements to draw from.
#' @param n Integer. Number of elements in the sample.
#'
#' @return Vector of length `n` of integers from the range `0` to `N-1`.
#'
#' @examples
#' srswr_rcpp_zero(10, 20)
#'
#' @export
srswr_rcpp_zero <- function(N, n) {
    .Call(`_pipsboot_srswr_rcpp_zero`, N, n)
}

