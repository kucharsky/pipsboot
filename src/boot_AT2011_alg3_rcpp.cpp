#include <Rcpp.h>
using namespace Rcpp;

#include "srswor_rcpp_pos.h"
#include "sampling_oneone_rcpp.h"

//' Antal and Tille (2011) bootstrap for simple random sampling without replacement
//'
//' This function implements method of resampling described in the algorithm 3 of
//' Antal & Tille (2011), p. 539.
//'
//' @param N Integer. Original population size.
//' @param n Integer. Sample size.
//' @param n_boots Integer. Number of bootstrap samples to be generated.
//'
//' @return Integer matrix with \code{n_boots} rows and \code{n} columns.
//' Each row is a bootstrap sample, with entries in columns denoting how many
//' times a given unit is resampled.
//'
//' @references Erika Antal and Yves Tille (2011) A Direct Bootstrap Method for
//' Complex Sampling Designs From a Finite Population, Journal of the American
//' Statistical Association, 106:494, 534-543, DOI: 10.1198/jasa.2011.tm09767
//'
//' @examples
//' boot_AT2011_alg3_rcpp(10, 5, 20)
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_AT2011_alg3_rcpp(int N, int n, int n_boots) {
  IntegerMatrix S (n_boots, n);
  if (n == 1) {
    std::fill(S.begin(), S.end(), 1);
  }
  if (n - n * n / N < 2) { // N < n + 2 + 4 / (n - 2)
    // n = 2 -> for every N
    // n = 3, 5 -> for N <= 8
    // n = 4 -> for N <= 7
    // n = 6, 7, ... -> for N <= n + 2
    // for "large" N, it holds only for n = 1, 2, N - 2, N - 1
    std::fill(S.begin(), S.end(), 1);
    double q = n * (N - n) / (2.0 * N);
    NumericVector ru = Rcpp::runif(n_boots, 0, 1);
    for(int k = 0; k < n_boots; k++) {
      if (ru[k] < q) {
        IntegerVector ij = srswor_rcpp_pos(n, 2);
        S(k, ij[0]) = 2;
        S(k, ij[1]) = 0;
      }
    }
  } else {
    double a = n * n / (1.0 * N);
    int mb = std::floor(a);
    int mbp = mb + 1;
    double q = mbp - a;
    int m;
    NumericVector ru = Rcpp::runif(n_boots, 0, 1);
    for(int k = 0; k < n_boots; k++) {
      if (ru[k] < q) {
        m = mb;
      } else {
        m = mbp;
      }
      IntegerVector v = srswor_rcpp_pos(n, m);
      for (int i = 0; i < m; i++) {
        S(k, v[i]) = 1;
      }
      IntegerVector w = sampling_oneone_rcpp(n - m);
      int j = 0;
      for (int i = 0; i < n; i++) {
        if (S(k, i) == 0) {
          S(k, i) = w[j++];
        }
      }
    }
  }
  return S;
}
