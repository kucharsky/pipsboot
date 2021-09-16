#' Calculate statistics based on bootstrap estimates
#'
#' @param estimator_value Numeric. Value of estimator based on a sample
#' @param bootstrap_estimates Vector. Values of estimator based on bootstrap resamples
#' @param confidence_level Numeric. Confidence level used for calculations of confidence interval limits.
#'
#' @return List with the following elements:
#'  - LCL_boot, UCL_boot: lower and upper confidence limits based on quantiles in bootstrap distribution,
#'  - LCL_normboot, UCL_normboot: lower and upper confidence limits based on normal distribution and bootstrap SE estimator
#'  - estVar, estSE: bootstrap variance and standard error of the estimator.
#'  - est: estimator value (the first argument).
#'
#' @examples
#' calculate_stats(rnorm(1), rnorm(1000), 0.95)
#'
#' @export
calculate_stats <- function(estimator_value, bootstrap_estimates, confidence_level) {
  # bootstrap variance
  estVar <- stats::var(bootstrap_estimates)

  # bootstrap standard error estimator
  estSE <- sqrt(estVar)

  # Lower Confidence Limit & Upper Confidence Limit based on quantiles in bootstrap distribution:
  p1 <- (1 - confidence_level) / 2
  p2 <- p1 + confidence_level
  qz_times_estSE <- stats::qnorm(p2) * estSE
  CL_boot <- stats::quantile(bootstrap_estimates, probs = c(p1, p2))

  list(LCL_boot = CL_boot[1],
       UCL_boot = CL_boot[2],
       # Lower & Upper Confidence Limits based on normal distribution and bootstrap SE estimator
       LCL_normboot = estimator_value - qz_times_estSE,
       UCL_normboot = estimator_value + qz_times_estSE,
       estVar = estVar,
       estSE = estSE,
       est = estimator_value)
}
