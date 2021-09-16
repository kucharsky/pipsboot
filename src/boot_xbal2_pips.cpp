#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(unwindProtect)]]

#include "inclusion_probabilities_rcpp.h"

//' Barbiero and Mecatti x-balanced 2 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
//'
//' The function generates bootstrap samples using Barbiero and Mecatti
//'   (2010, p. 63) x-balanced 2 bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling
//'
//' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
//' @param x_pop_total Number. Sum of the auxiliary variable in the population.
//' @param n_boots Integer. A number of bootstrap samples to be generated.
//' @param scheme Function. Implementation of a sampling scheme (takes a vector
//'   of inclusion probabilities as an argument and returns a sample: a vector of
//'   0's and 1's of the same length as argument).
//'
//' @return Integer matrix with `n_boots` rows and `length(x_sample)`
//'   columns. Each row is a bootstrap sample, and the number in i-th column
//'   denotes how many times a given unit is resampled.
//'
//' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
//'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
//'   in: Complex Data Modeling and Computationally Intensive Statistical
//'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
//'
//' @examples
//' boot_xbal2_pips_rcpp( 1:3, 100, 10, UPbrewer_rcpp)
//' boot_xbal2_pips_rcpp(1:10, 1000, 4, UPbrewer_rcpp)
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_xbal2_pips_rcpp(const NumericVector & x_sample,
                                   double x_pop_total,
                                   int n_boots,
                                   Function & scheme) {
  int n = x_sample.size();
  NumericVector inv_prs = x_pop_total / (x_sample * n);
  IntegerVector ci = no_init(n);
  NumericVector q = no_init(n);
  double sum_ci_times_x_sample = - x_pop_total;
  for(int i = 0; i < n; i++) {
    ci[i] = std::floor(inv_prs[i]);
    q[i] = inv_prs[i] / (inv_prs[i] + 1);       // calculate q
    sum_ci_times_x_sample += ci[i] * x_sample[i];
  }
  // determine order(r, decreasing = TRUE)
  IntegerVector idx = no_init(n); // index for order, declare without initialization
  std::iota(idx.begin(), idx.end(), static_cast<size_t>(0)); // fill with consecutive numbers, from 0
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return q[i] > q[j];});
  // order x_sample accordingly
  NumericVector x_sample_ord = no_init(n);
  // calculate cumsum(sort(x))
  NumericVector cum_x_sort = no_init(n);
  cum_x_sort[0] = x_sample[idx[0]];
  for(int i = 1; i < n; i++) {
    cum_x_sort[i] = cum_x_sort[i-1] + x_sample[idx[i]];
  }
  double min_pop_diff = 1e99;
  int which_min_pop_diff = n;
  NumericVector x_diff = no_init(n);
  for (int i = 0; i < n; i++) {
    x_diff[i] = std::abs(sum_ci_times_x_sample + cum_x_sort[i]);
    if (x_diff[i] < min_pop_diff) {
      min_pop_diff = x_diff[i];
      which_min_pop_diff = i;
    }
  }
  for(int i = 0; i <= which_min_pop_diff; i++) {
    ci[idx[i]] = ci[idx[i]] + 1;
  }
  int U_size = sum(ci);                  // size of pseudopopulation
  NumericVector x_b = no_init(U_size);   // x values in pseudopopulation
  IntegerVector group = no_init(U_size); // indicators of units from the original sample
  int k = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < ci[i]; j++) {
      x_b[k] = x_sample[i];
      group[k++] = i;
    }
  }
  // inclusion probabilities for pseudo-population
  NumericVector pr_b = inclusion_probabilities_rcpp(x_b, n);
  // matrix for results
  IntegerMatrix M (n_boots, n);
  // The bootstrap replications are here:
  for(int row = 0; row < n_boots; row++) {
    // sample from the pseudo-population
    IntegerVector sample = scheme(pr_b);
    // sum occurrences of units and put this into matrix row-wise
    for(int i = 0; i < U_size; i++) {
      M(row, group[i]) += sample[i];
    }
  }
  return M;
}
