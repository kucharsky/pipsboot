#' Precision of HT estimator of the population total
#'
#' The function computes bootstrap estimators of the variance and the standard
#' error of Horvitz-Thompson estimator of the population total as well as lower
#' and upper confidence limits.
#'
#' @param y_sample Numeric vector. Values of the variable of interest for sample elements.
#' @param pi_sample Numeric vector. Values of inclusion probabilities for sample elements.
#' @param bootstrap_samples Integer matrix with `B` rows and `n` columns, where
#'   `B` is the number of bootstrap replications and `n` is the sample size.
#'   Each row is a bootstrap sample, and the number in i-th column denotes the
#'   number of times i-th unit from the original sample occurs in the bootstrap sample.
#' @param confidence_level The assumed confidence level.
#'
#' @return A list with the following elements:
#'   - LCL_boot: Lower Confidence Limit based on the quantile of bootstrapped values of the estimator,
#'   - UCL_boot: Upper Confidence Limit based on the quantile of bootstrapped values of the estimator.
#'   - LCL_normboot: Lower Confidence Limit based on the quantile of the standard normal distribution and the bootstrap standard error estimator.
#'   - UCL_normboot: Upper Confidence Limit based on the quantile of the standard normal distribution and the bootstrap standard error estimator.
#'   - estVar_HT: bootstrap variance estimator of the HT estimator of the population total.
#'   - estSE_HT: bootstrap standard error estimator of the HT estimator of the population total.
#'   - HT: value of HT estimator of the population total based on original sample data.
#'
#' @examples
#' x_sample <- 1:10
#' y_sample <- x_sample + rnorm(10)
#' pi_sample <- (1:10) / 100
#' set.seed(1)
#' bootstrap_samples <- boot_H_pips(x_sample, 1000, 4, UPbrewer_rcpp)
#' HT_boot(y_sample, pi_sample, bootstrap_samples, 0.95)
#'
#' @export
HT_boot <- function(y_sample, pi_sample, bootstrap_samples, confidence_level) {
  if (!is.numeric(y_sample)) {
    stop(paste("The first argument of HT_boot function should be a numeric vector, but we got:\n",
               paste(y_sample, collapse = " ")))
  }
  if (!is.numeric(pi_sample) || !all(pi_sample >= 0 & pi_sample <= 1)) {
    stop(paste("The second argument of HT_boot function should be a vector of probabilitites, but we got:\n",
               paste(pi_sample, collapse = " ")))
  }
  if (!is.numeric(confidence_level) || !(length(confidence_level) == 1) || !(confidence_level >= 0 & confidence_level <= 1)) {
    stop(paste("The fourth argument of HT_boot function should be a probability but we got:\n",
               paste(confidence_level, collapse = " ")))
  }
  if (!(length(y_sample) == ncol(bootstrap_samples))) {
    stop(paste("The number of columns of 'bootstrap_samples' must be equal to the length of 'y_sample' vector"))
  }
  if (!(length(y_sample) == length(pi_sample))) {
    stop(paste("Vectors 'y_sample' and 'pi_sample' should be of the same length"))
  }

  y_over_pi <- y_sample / pi_sample

  # value of HT estimator of the population total based on original sample data
  HT <- sum(y_over_pi)

  # bootstrapped values of HT estimator of the population total
  HT_b <- drop(bootstrap_samples %*% y_over_pi)

  calculate_stats(HT, HT_b, confidence_level)
}
