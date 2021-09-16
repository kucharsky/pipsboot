#include <Rcpp.h>
using namespace Rcpp;

//' Bootstrap replications of sample data
//'
//' The function computes a matrix with bootstrap replications of sample data.
//'
//' @param bootstrap_samples Integer matrix with `B` rows and `n` columns,
//'   where `B` is the number of bootstrap replications and `n` is the sample
//'   size. Each row is a bootstrap sample, with entries in columns denoting how
//'   many times given unit is resampled.
//' @param x_sample Numeric vector. Values of the variable for sampled elements.
//'
//' @return Matrix with `B` rows and `n` columns. Each row is a bootstrap sample of values of `x_sample`.
//'
//' @examples
//' x_sample <- 1:10
//' set.seed(111)
//' bootstrap_samples <- boot_H_pips(x_sample, 1000, 4, UPbrewer_rcpp)
//' sampledata_bootrep_rcpp(bootstrap_samples, x_sample)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix sampledata_bootrep_rcpp(const IntegerMatrix & bootstrap_samples,
                                       const NumericVector & x_sample) {
  std::size_t B = bootstrap_samples.nrow();
  std::size_t n = bootstrap_samples.ncol();
  NumericMatrix X = no_init(B, n);
  for (std::size_t row = 0; row < B; row++) {
    int k = 0;
    for (std::size_t i = 0; i < n; i++) { // n = sum(bootstrap_samples(row, _))
      for (int j = 0; j < bootstrap_samples(row, i); j++)
        X(row, k++) = x_sample[i];
    }
  }
  return X;
}
