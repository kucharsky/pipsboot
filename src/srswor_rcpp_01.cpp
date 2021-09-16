#include <Rcpp.h>
using namespace Rcpp;

//' Simple Random Sampling Without Replacement
//'
//' Function draws a simple random sample without replacement. This is an Rcpp
//' replacement for \code{\link[sampling]{srswor}}, however with different
//' implementation.
//'
//' @param N Integer. Population size - number of elements to draw from.
//' @param n Integer. Sample size - number of elements to draw.
//'
//' @return Integer vector of size \code{N} with exactly \code{n} elements
//' equal to 1 (units selected into the sample) and others equal to 0 (units not
//' selected).
//'
//' @examples
//' srswor_rcpp_01(10, 5)
//'
//' @export
// [[Rcpp::export]]
IntegerVector srswor_rcpp_01(const int & N, const int & n) {
  int nn = n;
  if (2 * n > N) {
    // If we have to draw more than a half of the population,
    // then it is faster to draw elements that will not go into the sample.
    nn = N - n;
  }
  IntegerVector S (N, 0);
  int k = 0; // number of already sampled elements
  while (k < nn) {
    int kk = nn - k; // how many more units we have to draw
    IntegerVector r = as<IntegerVector>(Rcpp::runif(kk, 0, N)); // drawing
    for (int i = 0; i < kk; i++) {
      k += 1 - S[r[i]]; // if unit was not selected before, increment k
      S[r[i]] = 1; // mark as selected anyhow
    }
  }
  if (2 * n > N) { // in this case simply revert selection
    // this step can be avoided if we would accordingly mark elements in lines 34-35
    for (int i = 0; i < N; i++) {
      S[i] = 1 - S[i];
    }
  }
  return S;
}

