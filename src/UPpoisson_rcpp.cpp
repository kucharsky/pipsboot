#include <Rcpp.h>
using namespace Rcpp;

//' Poisson sampling
//'
//' Draws a Poisson sample (without replacement, random sample size) using a
//' given vector of first-order inclusion probabilities. Rewriten in Rcpp based
//' on the `\link[sampling]{UPpoisson}` function.
//'
//' @param pik  Numeric vector. First-order inclusion probabilities.
//'
//' @return Integer vector of size \code{length(pik)} (population size) with
//' elements 1 (if the unit is selected in the sample) or 0 (if it is not).
//'
//' @examples
//' UPpoisson_rcpp((1:10) / 55 * 5)
//' UPpoisson_rcpp(inclusion_probabilities_rcpp(runif(10), 5))
//'
//' @export
// [[Rcpp::export]]
IntegerVector UPpoisson_rcpp(const NumericVector & pik) {
  int n = pik.size();
  IntegerVector s (n);
  NumericVector u = Rcpp::runif(n, 0, 1);
  for(int i = 0; i < n; i++){
    s[i] = (u[i] < pik[i]) ? 1 : 0;
  }
  return s;
}
