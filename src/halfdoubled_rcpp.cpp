#include <Rcpp.h>
using namespace Rcpp;

#include "srswor_rcpp_01.h"

//' Doubled half sampling
//'
//' Draws a sample using "doubled half sampling" described in
//' Antal and Tille (2014), pp. 1351-1352. Samples are drawn with replacement.
//'
//' @param n Integer. Number of elements to sample from.
//'
//' @return Integer vector of size \code{n}. Values indicates how many times
//'         units are selected into sample, 0 means that the unit is not
//'         in the sample.
//'
//' @references Antal, E. and Tille, Y. (2014) A new resampling method for
//' sampling designs without replacement: the doubled half bootstrap.
//' Comput Stat 29, 1345-1363, https://doi.org/10.1007/s00180-014-0495-0
//'
//' @examples
//' halfdoubled_rcpp(5)
//' halfdoubled_rcpp(20)
//'
//' @export
// [[Rcpp::export]]
IntegerVector halfdoubled_rcpp(const int & n) {
  IntegerVector s (n);
  if (n % 2 == 0) {
    s = 2 * srswor_rcpp_01(n, n / 2); // even case is simple
  } else {
    s = 2 * srswor_rcpp_01(n, (n - 1) / 2); // odd case
    if (R::runif(0, 1) < 0.25) { // sometimes we do something
      int r = (int) R::runif(0, (n - 1) / 2); // this is srswor(1, (n - 1) / 2)
      int i = 0;
      int k = 0;
      while (k < r) {
        if (s[i++] == 2) {
          k++;
        }
      }
      while (k == r) {
        if (s[i] == 2) {
          s[i]++;
          break;
        }
        i++;
      }
    } else { // and sometimes we do something different
      int r = (int) R::runif(0, (n + 1) / 2); // this is srswor(1, (n + 1) / 2)
      int i = 0;
      int k = 0;
      while (k < r) {
        if (s[i++] == 0) {
          k++;
        }
      }
      while (k == r) {
        if (s[i] == 0) {
          s[i]++;
          break;
        }
        i++;
      }
    }
  }
  return s;
}
