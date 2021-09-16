#include <Rcpp.h>
using namespace Rcpp;

//' Simple Random Sampling Without Replacement
//'
//' Function draws a simple random sample without replacement. This is a fast
//' replacement for R's \code{sample(N, n, replace = FALSE)}.
//'
//' @param N Integer. Population size - number of elements to draw from.
//' @param n Integer. Sample size - number of elements to draw.
//'
//' @return Integer vector of length \code{n} of indices of selected units.
//' Indices have values from the range 0 to N-1 (C++ convention)
//'
//' @examples
//' srswor_rcpp_pos(10, 5)
//'
//' @export
// [[Rcpp::export]]
IntegerVector srswor_rcpp_pos(const int & N, const int & n) {
  // first part is exactly the same as in srswor_rcpp_01
  int nn = n;
  if (2 * n > N) {
    // If we have to draw more than a half of the population,
    // then it is faster to draw elements that will not go into the sample.
    nn = N - n;
  }
  IntegerVector S (N); // 0-1 vector
  int k = 0;
  while (k < nn) {
    int kk = nn - k; // how many more units we have to draw
    IntegerVector r = as<IntegerVector>(Rcpp::runif(kk, 0, N)); // drawing
    for (int i = 0; i < kk; i++) {
      k += 1 - S[r[i]]; // if unit was not selected before, increment k
      S[r[i]] = 1; // mark as selected anyhow
    }
  }
  IntegerVector K (n); // vector of indices
  int collect = 1; // "default": collect indices if S == 1
  if (2 * n > N) {
    collect = 0; // collect indices if S == 0
  }
  k = 0;
  int i = -1;
  while (k < n) {
    if (S[++i] == collect) {
      K[k++] = i;
    }
  }
  return K;
}
