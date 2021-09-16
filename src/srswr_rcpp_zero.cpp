#include <Rcpp.h>
using namespace Rcpp;

//' Simple Random Sampling With Replacement
//'
//' Samples `n` integers from the range `0` to `N-1`.
//' A call `srswr_rcpp_zero(N, n)` is statistically equivalent to R's
//' `sample(N, n, replace = TRUE) - 1` however the results of every individual
//' call are different due to differences in the implementation.
//' `srswr_rcpp_zero` is up to 5 times faster.
//'
//' @param N Integer. Number of elements to draw from.
//' @param n Integer. Number of elements in the sample.
//'
//' @return Vector of length `n` of integers from the range `0` to `N-1`.
//'
//' @examples
//' srswr_rcpp_zero(10, 20)
//'
//' @export
// [[Rcpp::export]]
IntegerVector srswr_rcpp_zero(const int & N, const int & n) {
  return as<IntegerVector>(Rcpp::runif(n, 0.0, N));
}
