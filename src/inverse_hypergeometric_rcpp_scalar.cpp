#include <Rcpp.h>
using namespace Rcpp;

//' Inverse Hypergeometric random numbers
//'
//' Generates random number from a special case of the inverse hypergeometric distribution.
//'
//' Generates random number from a discrete distribution with the probability mass function:
//' \deqn{\mathbb{P}(X = j) = \dfrac{\dbinom{N + n - j - 2}{n - j}}{\dbinom{N + n - 1}{n}} \qquad j = 0,\dots,n,}
//' being a special case of the inverse hypergeometric distribution.
//'
//' This distribution may be described as follows: we have `n` red balls and `N-1` black balls.
//' We draw balls without replacement, and `X` is the number of red balls drawn before first black ball.
//'
//' The expectation of this distribution is \eqn{\mathbb{E}(X) = n/N}
//' and variance \eqn{\mathrm{Var}(X) = \dfrac{n(N-1)(N+n)}{N^2(N+1)}}.
//' Corner cases are covered in the following way:
//' if \eqn{n = 0} then \eqn{X \equiv 0} and if \eqn{N = 1} then \eqn{X \equiv n}.
//'
//' @param N Positive integer
//' @param n Vector of nonnegative integers.
//'
//' @return Integer.
//'
//' @examples
//' inverse_hypergeometric_rcpp_scalar(10, 20)
//'
//' @export
// [[Rcpp::export]]
int inverse_hypergeometric_rcpp_scalar(const int & N, const int & n) {
  if (N == 1) {
    return n;
  } else {
    double u = R::runif(0, 1);
    double d = (N - 1) / (N + n - 1.0);
    int j = 0;
    double cum = d;
    while (u >= cum) {
      d = d * (n - j) / (N + n - j - 2.0);
      cum += d;
      j++;
    }
    return j;
  }
}
