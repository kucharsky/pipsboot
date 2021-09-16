#include <Rcpp.h>
using namespace Rcpp;

#include <math.h>
#include "srswr_rcpp_zero.h"
#include "sampling_with_overreplacement_rcpp_scalar.h"

//' One-one resampling design
//'
//' Implements sampling design described in
//' Antal and Tille (2011, page 537, Algorithm 2).
//'
//' @param n Integer. Sample size.
//'
//' @return Integer vector of size `n. The value of `i`-th element tells how
//' many times `i`-th unit appears in the sample.
//'
//' @references Erika Antal and Yves Tille (2011),
//' A Direct Bootstrap Method for Complex Sampling Designs From a Finite
//' Population, Journal of the American Statistical Association 106:494,
//' 534-543, DOI: 10.1198/jasa.2011.tm09767
//'
//' @details
//' Fixed sample size resampling design with replacement. Citing
//' Antal and Tille (2011) p. 537: "an ad hoc construction used to randomly
//' select `n` units from a sample of size `n` in such a way that the
//' expectation and the variance of \eqn{S^*_k} are equal to 1 (...). This
//' sampling design is a mixture between a simple random sampling with
//' replacement and a simple random sampling with  over-replacement".
//'
//' @examples
//' sampling_oneone_rcpp(10)
//'
//' @export
// [[Rcpp::export]]
IntegerVector sampling_oneone_rcpp(const int & n) { // version scalar_srswr
  IntegerVector S (n);
  if (n == 2) {
    if (R::runif(0, 1) < 0.5) {
      S = {2, 0};
    } else {
      S = {0, 2};
    }
  } else {
    int m = std::floor((1 + std::sqrt(4 * n + 9 + 8 / (n - 1.0))) / 2.0);
    double alpha = (m + 1 - n * (n + 1) / (m * (n - 1.0))) / 2.0;
    if (R::runif(0, 1) >= alpha) {
      m++;
    }
    S = sampling_with_overreplacement_rcpp_scalar(n, m);
    int nnt = n - m;
    IntegerVector v = srswr_rcpp_zero(n, nnt);
    for (int i = 0; i < nnt; i++) {
      S[v[i]] += 1;
    }
  }
  return S;
}
