#include <Rcpp.h>
using namespace Rcpp;

#include "inverse_hypergeometric_rcpp_scalar.h"

//' Sampling with over-replacement
//'
//' Implements sampling design described in Antal and Tille (2011),
//' and called there "Simple random sampling with over-replacement".
//'
//' @param N Integer. Size of the population.
//' @param n Integer. Number of elements in the sample.
//'
//' @return Integer vector of length `N`.
//' The value of `i`-th coordinate tells how many times `i`-th unit appears
//' in the sample.
//'
//' @references E. Antal and Y. Tille (2011)
//' Simple random sampling with over-replacement,
//' Journal of Statistical Planning and Inference 141, 597-601.
//'
//' @details
//' A sample of fixed size `n` is drawn with replacement from population of
//' size `N`. All samples have the same probability of being selected:
//' \eqn{\left(\dbinom{N + n - 1}{n}\right)^{-1}}.
//'
//' After Antal and Tille (2011) we describe the derivation of this design.
//' Consider a sequence of `N` independent geometric random variables
//' \eqn{X_1, \dots, X_N}. Conditioning on \eqn{X_1 + \dots + X_N = n}, we
//' obtain a distribution such that
//' \deqn{P(X_1 = x_1, \dots, X_N = x_n) = 1 / \dbinom{N+n-1}{n}} for every
//' vector of nonnegative numbers
//' \eqn{(x_1,\dots, x_N)} with \eqn{x_1 + \dots + x_N = n}. Marginals of this
//' distribution are inverse hypergeometric. We have
//' \eqn{\mathbb{E}(X_k) = n/N},
//' \eqn{\mathrm{var}(X_k) = \dfrac{n(N-1)(N+n)}{N^2(N+1)}} and
//' \eqn{\mathrm{cov}(X_k, X_l) = -\mathrm{var}(X_k)/(N-1)} if \eqn{k\not=l}.
//'
//' @examples
//' sampling_with_overreplacement_rcpp_scalar(10, 8)
//'
//' @export
// [[Rcpp::export]]
IntegerVector sampling_with_overreplacement_rcpp_scalar(const int & N, const int & n) {
  IntegerVector S (N);
  int nk = n;
  for (int k = 0; k < N; k++) {
    S[k] = inverse_hypergeometric_rcpp_scalar(N - k, nk);
    nk = nk - S[k];
  }
  return S;
}
