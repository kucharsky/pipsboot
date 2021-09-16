#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(unwindProtect)]]

#include "inclusion_probabilities_rcpp.h"
#include "sampling_with_overreplacement_rcpp_scalar.h"
#include "sampling_oneone_rcpp.h"

// #include <functional> // for passing functions as parameters as follows:
// std::function <return_type (input_type)> func_as_arg_name
// eg: std::function <IntegerVector (int, int)> scheme

// #include "fill_samples.h"
// helper function for boot_AT2011_alg45_rcpp.cpp
inline void fill_samples (IntegerMatrix & S,
                          IntegerVector & w,
                          const int & Scol,
                          const int & ns,
                          const int & row) {
  for(int i = 0; i < Scol; i++) {
    S(row, i) = w[i];
  }
  // replace 0's in S with values of v
  IntegerVector v = sampling_oneone_rcpp(ns);
  int k = 0;
  int i = 0;
  while (k < ns) {
    if (w(i) == 0) {
      S(row, i) = v[k++];
    }
    i++;
  }
}

//' Antal and Tille (2011) bootstrap for complex sampling designs
//'
//' Function implements a method of resampling from a sample selected with unequal
//' probabilities without replacement as described in algorithms 4 and 5 in
//' Antal & Tille (2011), pp. 539-540.
//'
//' @param pi_sample Numeric vector. Inclusion probabilities for elements in the original sample.
//' @param n_boots Integer. A number of bootstrap samples to be generated.
//' @param scheme Function. Implementation of a sampling scheme (takes a vector
//'   of inclusion probabilities as an argument and returns a sample: a vector of
//'   0's and 1's of the same length as the argument) used to select the original sample.
//'
//' @return Integer Matrix with \code{n_boots} rows and \code{length(pi_sample)}
//' columns. Each row is a bootstrap sample with entries denoting how many
//' times given unit is resampled.
//'
//' @references Erika Antal and Yves Tille (2011) A Direct Bootstrap Method for
//' Complex Sampling Designs From a Finite Population, Journal of the American
//' Statistical Association, 106:494, 534-543, DOI: 10.1198/jasa.2011.tm09767
//'
//' @examples
//' pi_sample <- inclusion_probabilities_rcpp(1:10, 6)
//' boot_AT2011_alg45_rcpp(pi_sample, 5, UPbrewer_rcpp)
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_AT2011_alg45_rcpp(const NumericVector & pi_sample,
                                     const int & n_boots,
                                     const Rcpp::Function & scheme) {
  int n = pi_sample.size();
  IntegerMatrix S = no_init(n_boots, n);
  double n_star = 0;
  for (int i = 0; i < n; i++) {
    n_star += pi_sample[i];
  }
  double n_rem = n - n_star;
  NumericVector psik (n);
  if (n_rem >= 2) { // Algorithm 4
    if (abs(n_star - std::round(n_star)) < 1e-6) {  // n_star is an integer
      // psik = pi_sample; ns = n_rem;
      int ns = (int) std::round(n_rem);
      for (int row = 0; row < n_boots; row++) {
        IntegerVector w = scheme(pi_sample);
        fill_samples(S, w, n, ns, row);
      } // end of `row` loop
    } else {                                        // n_star is not an integer
      int m2 = std::floor(n_star) + 1;
      double q = m2 - n_star;
      // increase phik2 iteratively to have: sum(phik2) == m2
      NumericVector phik2 = inclusion_probabilities_rcpp(pi_sample, m2);
      // calculate phik1 from the equation: q * phik1 + (1-q) * phik2 = pi_sample
      NumericVector phik1 = (pi_sample - (1 - q) * phik2) / q;
      for(int i = 0; i < n; i++) { // for some strange numerical reason they may be negative
        if (phik1[i] < 0) {
          phik1[i] = 0;
        }
      }
      int ns2 = n - m2;
      int ns1 = ns2 + 1; // = n - (m2 - 1)
      NumericVector u = Rcpp::runif(n_boots, 0, 1);
      for (int row = 0; row < n_boots; row++) {
        if (u[row] < q) {
          // psik = phik1; ns = ns1;
          IntegerVector w = scheme(phik1);
          fill_samples(S, w, n, ns1, row);
        } else {
          // psik = phik2; ns = ns2;
          IntegerVector w = scheme(phik2);
          fill_samples(S, w, n, ns2, row);
        }
      }
    } // end of "else" `n_star is not an integer
  } else { // Algorithm 5
    NumericVector u = Rcpp::runif(n_boots, 0, 1);
    double half_n_rem = n_rem / 2.0;
    NumericVector psik = 1 - inclusion_probabilities_rcpp(1 - pi_sample, 2);
    for(int row = 0; row < n_boots; row++) {
      if (u[row] > half_n_rem) { // "q" for this case is hidden here
        for(int i = 0; i < n; i++) {
          S(row, i) = 1;
        }
      } else {
        IntegerVector w = scheme(psik);
        fill_samples(S, w, n, 2, row);
      }
    }
  }
  return S;
}
