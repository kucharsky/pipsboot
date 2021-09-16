#' Precision of HT estimator of the domain total
#'
#' The function computes bootstrap estimators of the variance and the standard
#' error of Horvitz-Thompson estimator of the domain total as well as lower and
#' upper confidence limits.
#'
#' @param y_sample Numeric vector. Values of the variable of interest for sample elements.
#' @param pi_sample Numeric vector. Values of inclusion probabilities for sample elements.
#' @param domain_sample Identifiers of domains in the sample.
#' @param d Identifier of the domain of interest.
#' @param bootstrap_samples Integer matrix with `B` rows and `n` columns, where
#'   `B` is the number of bootstrap replications and `n` is the sample size.
#'   Each row is a bootstrap sample, and the number in i-th column denotes the
#'   number of times i-th unit from the original sample occurs in the bootstrap sample.
#' @param confidence_level The assumed confidence level.
#'
#' @return A list with the following elements:
#'   - LCL_boot Lower Confidence Limit based on the quantile of bootstrapped values of the estimator.
#'   - UCL_boot Upper Confidence Limit based on the quantile of bootstrapped values of the estimator.
#'   - LCL_normboot Lower Confidence Limit based on the quantile of the standard normal distribution and the bootstrap standard error estimator.
#'   - UCL_normboot Upper Confidence Limit based on the quantile of the standard normal distribution and the bootstrap standard error estimator.
#'   - estVar_HT bootstrap variance estimator of the HT estimator of the domain total.
#'   - estSE_HT bootstrap standard error estimator of the HT estimator of the domain total.
#'   - HT: value of HT estimator of the domain total based on original sample data.
#'
#' @examples
#' x_sample <- 1:10
#' y_sample <- x_sample + rnorm(10)
#' pi_sample <- (1:10) / 100
#' # First seven observations in the sample are from domain 1,
#' # last three observations are from domain 2:
#' domain_sample <- c(rep(1, 7), rep(2, 3))
#' set.seed(1)
#' bootstrap_samples <- boot_H_pips(x_sample, 1000, 4, UPbrewer_rcpp)
#' HTd_boot(y_sample,pi_sample,domain_sample, 2, bootstrap_samples, 0.95)
#' set.seed(111)
#' bootstrap_samples <- boot_H_pips(x_sample, 1000, 4, UPbrewer_rcpp)
#' bootstrap_samples
#' # We have obtained 0's in the first row and in the last three columns of
#' # 'bootstrap_samples'- the second domain (d=2) is not observed
#' # in the first bootstrap sample. Hence, we will obtain NA's:
#' HTd_boot(y_sample, pi_sample, domain_sample, 2, bootstrap_samples, 0.95)
#'
#' @export
HTd_boot <- function(y_sample, pi_sample, domain_sample, d, bootstrap_samples, confidence_level){
  if (!is.numeric(y_sample)) {
    stop(paste("The first argument of HTd_boot function should be a numeric vector, but we got:\n",
               paste(y_sample, collapse = " ")))
  }
  if (!is.numeric(pi_sample) || !all(pi_sample >= 0 & pi_sample <= 1)) {
    stop(paste("The second argument of HTd_boot function should be a vector of probabilitites, but we got:\n",
               paste(pi_sample, collapse = " ")))
  }
  if (!(length(d) == 1)) {
    stop(paste("The length of the fourth argument of HTd_boot function should be 1 but we got:\n",
               paste(d, collapse = " ")))
  }
  if (!is.numeric(confidence_level) || !(length(confidence_level) == 1) || !(confidence_level >= 0 & confidence_level <= 1)) {
    stop(paste("The sixth argument of HTd_boot function should be a probability but we got:\n",
               paste(confidence_level, collapse = " ")))
  }
  if (!(length(y_sample) == ncol(bootstrap_samples))) {
    stop(paste("The number of columns of 'bootstrap_samples' must be equal to the length of 'y_sample' vector"))
  }
  if (!(length(y_sample) == length(pi_sample))  ||
      !(length(pi_sample) == length(domain_sample)) ||
      !(length(y_sample) == length(domain_sample))) {
    stop(paste("Vectors 'y_sample' and 'pi_sample' and 'domain_sample' should be of the same length"))
  }

  d_domain <- as.integer(domain_sample == d)
  HT_b_d_ind <- bootstrap_samples %*% d_domain

  if (any(HT_b_d_ind == 0)) {
    print("In at least one of the bootstrap samples no elements of the domain of interest are observed")
    list(LCL_boot = NA,
         UCL_boot = NA,
         LCL_normboot = NA,
         UCL_normboot = NA,
         estVar_HT = NA,
         estSE_HT = NA)
  } else {
    d_y_over_pi <- d_domain * y_sample / pi_sample

    # value of HT estimator of the domain total based on original sample data
    HT <- sum(d_y_over_pi)

    # bootstrapped values of HT estimator of the domain total
    HT_b_d <- drop(bootstrap_samples %*% d_y_over_pi)

    calculate_stats(HT, HT_b_d, confidence_level)
  }
}
