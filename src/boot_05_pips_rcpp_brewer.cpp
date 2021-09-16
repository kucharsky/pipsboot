#include <Rcpp.h>
using namespace Rcpp;

#include <math.h>
#include "inclusion_probabilities_rcpp.h"
#include "UPbrewer_rcpp.h"

//' Barbiero and Mecatti 0.5 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
//'
//' The function generates bootstrap samples using Barbiero and Mecatti
//'   (2010, p. 62) "0.5 \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS algorithm".
//'   Uses Brewers' design for sampling.
//'
//' @param x_sample Numeric vector. Positive values of an auxiliary variable for
//'   sample elements.
//' @param x_pop_total Number. Sum of the auxiliary variable in the
//'   population.
//' @param n_boots Integer. A number of bootstrap samples to be generated.
//' @return Integer matrix with `n_boots` rows and `length(x_sample)`
//'   columns. Each row is a bootstrap sample, and the number in i-th column
//'   denotes the number of times i-th unit from the original sample occurs in
//'   the bootstrap sample.
//' @export
//' @references Barbiero A. and Mecatti F. (2010). Bootstrap algorithms for
//'   variance estimation in \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling,
//'   in: Complex Data Modeling and Computationally Intensive Statistical
//'   Methods, ed. P. Mantovan and P. Secchi, Springer, Milano.
//' @examples
//' boot_05_pips_rcpp(1:3, 100, 10, UPbrewer_rcpp)
//' boot_05_pips_rcpp(1:10, 1000, 5, UPbrewer_rcpp)
// [[Rcpp::export]]
IntegerMatrix boot_05_pips_rcpp_brewer(const NumericVector & x_sample,
                                       const double & x_pop_total,
                                       const int & n_boots) {
  int n = x_sample.size();   // sample size
  double xptn = x_pop_total / n;
  // number of repetitions in pseudo-population
  IntegerVector ni (n);
  int nxb = 0;
  for(int i = 0; i < n; i++) {
    ni[i] = (int) round(xptn / x_sample[i]);
    nxb += ni[i];
  }
  NumericVector x_b (nxb);   // x values in pseudo-population
  int k = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < ni[i]; j++) {
      x_b[k++] = x_sample[i]; // replicate i-th element of x_sample ni[i] times
    }
  }
  // inclusion probabilities for pseudo-population
  NumericVector pr_b = inclusion_probabilities_rcpp(x_b, n);
  IntegerMatrix M (n_boots, n);   // matrix for results
  for(int row = 0; row < n_boots; row++) {
    IntegerVector s = UPbrewer_rcpp(pr_b);
    k = 0;
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < ni[i]; j++) {
        M(row, i) += s[k++]; // sum-up occurences in the sample
      }
    }
  }
  return M;
}
