#include <Rcpp.h>
using namespace Rcpp;

#include "inclusion_probabilities_rcpp.h"

//' Sverchkov and Pfeffermann \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
//'
//' The function generates bootstrap samples using Sverchkov and Pfeffermann
//'   (2004, p. 85) bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling
//'
//' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
//' @param pi_sample Numeric vector. Values of inclusion probabilities for sample elements.
//' @param N Integer. Population size.
//' @param n_boots Integer. A number of bootstrap samples to be generated.
//' @param scheme Function. Implementation of a sampling scheme (takes a vector of inclusion probabilities as an argument and returns a sample: a vector of 0's and 1's of the same length as argument).
//'
//' @return Integer matrix with `n_boots` rows and `length(x_sample)` columns. Each row is a bootstrap sample, and the number in i-th column denotes the number of times i-th unit from the original sample occurs in the bootstrap sample.
//'
//' @references Sverchkov M. and Pfeffermann D. (2004). Prediction of Finite Population Totals Based on the Sample Distribution, \emph{Survey Methodology}, 30(1):79-92.
//'
//' @examples
//' boot_SP_pips_rcpp(1:3, c(0.005, 0.01, 0.015), 100, 10, UPbrewer_rcpp)
//' boot_SP_pips_rcpp(1:10, (1:10) * 2e-5, 1000, 4, UPbrewer_rcpp)
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_SP_pips_rcpp(const NumericVector & x_sample,
                                const NumericVector & pi_sample,
                                const int & N,
                                const int & n_boots,
                                const Rcpp::Function & scheme) {
  int n = x_sample.size(); // sample size
  // ni = a number of repetitions in pseudo-population
  // don't try to use "sample", because we need the counts
  // ni <- sampling::UPmultinomial(N / sum(inv_pi_sample) * inv_pi_sample)
  // Obtaining namespace of Matrix package
  // Environment pkg = Environment::namespace_env("Matrix");
  // Picking up Matrix() function from Matrix package
  // Function rmultinom = pkg["Matrix"];
  Function rmultinom("rmultinom");
  IntegerMatrix ni = rmultinom(1, N, 1/pi_sample); // that's what sampling::UPmultinomial really is
  int nxb = 0;
  for(int i = 0; i < n; i++) {
    nxb += ni(i, 0);
  }
  NumericVector x_b (nxb);   // x values in pseudo-population
  int k = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < ni(i, 0); j++) {
      x_b[k++] = x_sample[i]; // replicate i-th element of x_sample ni[i] times
    }
  }
  // inclusion probabilities for pseudo-population
  NumericVector pr_b = inclusion_probabilities_rcpp(x_b, n);
  IntegerMatrix M (n_boots, n);   // matrix for results
  for(int row = 0; row < n_boots; row++) {
    IntegerVector s = scheme(pr_b);
    k = 0;
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < ni[i]; j++) {
        M(row, i) += s[k++]; // sum-up occurences in the sample
      }
    }
  }
  return M;
}
