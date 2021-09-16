#include <Rcpp.h>
using namespace Rcpp;

//' Quatember \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
//'
//' The function generates bootstrap samples using Quatember (2014, p. 95)
//' bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling.
//'
//' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
//' @param x_pop_total Numeric. Sum of the auxiliary variable in the population.
//' @param n_boots Integer. A number of bootstrap samples to be generated.
//'
//' @return Integer matrix with `n_boots` rows and `length(x_sample)`
//'   columns. Each row is a bootstrap sample, and the number in i-th column
//'   denotes how many times a given unit is resampled.
//'
//' @references Quatember, A. (2014). The Finite Population Bootstrap - From the
//' Maximum Likelihood to the Horvitz-Thompson Approach. Austrian Journal of
//' Statistics, 43(2), 93-102. https://doi.org/10.17713/ajs.v43i2.10
//'
//' @examples
//' boot_Q_pips_rcpp(1:3, 100, 10)
//' boot_Q_pips_rcpp(1:10, 1000, 4)
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_Q_pips_rcpp(const NumericVector & x_sample,
                               const double & x_pop_total,
                               int n_boots) {
  const int n = x_sample.size();
  IntegerMatrix M (n_boots, n);
  // we implement formula (6) on page 97, nominator and denominator divided by n
  const double tx_over_n = x_pop_total / n;
  for(int row = 0; row < n_boots; row++) {
    NumericVector u = Rcpp::runif(n, 0, 1);
    double denominator = x_pop_total;
    for (int j = 0; j < n; j++) {
      int k = -1;
      double sn = 0;
      double tresh = denominator * u[j];
      while(sn < tresh) {
        k++;
        double numerator = tx_over_n - M(row, k) * x_sample[k];
        if (numerator > 0) {
          sn += numerator;
        }
      }
      M(row, k) += 1;
      denominator -= x_sample[k];
    }
  }
  return M;
}
