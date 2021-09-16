#include <Rcpp.h>
using namespace Rcpp;

//' Generalized Quatember \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
//'
//' The function generates bootstrap samples using generalized Quatember (2014, p. 95)
//' bootstrap \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS sampling where additional weights
//' are taken into account to compute the number of replicates in the pseudo-population.
//'
//' @param x_sample Numeric vector. Positive values of an auxiliary variable for sample elements.
//' @param x_pop_total Number. Sum of the auxiliary variable in the population.
//' @param g_weights Numeric vector of the length `length(x_sample)`. Positive
//' weights for sample elements. Should sum up to the sample size to ensure that
//' products of 'x_sample' and calibration weights (defined as products of
//' 'g_weights` and inverses of first order inclusion probabilities) are equal
//' to 'x_pop_total'. If this condition is not met, they are rescaled.
//' @param n_boots Integer. A number of bootstrap samples to be generated.
//'
//' @return Integer matrix with `n_boots` rows and `length(x_sample)`
//'   columns. Each row is a bootstrap sample, with entries in columns denoting
//'   how many times given unit is resampled.
//'
//' @references Quatember, A. (2014). The Finite Population Bootstrap - From the
//' Maximum Likelihood to the Horvitz-Thompson Approach. Austrian Journal of
//' Statistics, 43(2), 93-102. https://doi.org/10.17713/ajs.v43i2.10
//'
//' @references Zadlo T. (2021) On the generalisation of Quatember's
//' bootstrap, Statistics in Transition new series, 22(1), 163-178,
//' https://doi.org/10.21307/stattrans-2021-009.
//'
//' @examples
//' N <- 200 # population size
//' x_population <- rgamma(N, 1)
//' n <- 10 # sample size
//' pi <- inclusion_probabilities_rcpp(x_population, n)
//' set.seed(123)
//' sample <- UPbrewer_rcpp(pi)
//' x_s <- x_population[sample == 1]
//' pi_s <- pi[sample == 1]
//' # g-weigths
//' g <- sampling::calib(Xs = cbind(x_s, rep(1, n)),
//'                      d = 1 / pi_s,
//'                      total = c(sum(x_population), N),
//'                      method = "truncated",
//'                      bounds = c(low = 0, upp = 10))
//' boot_Qg_pips_rcpp(x_s, sum(x_population), g, 5)
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_Qg_pips_rcpp(const NumericVector & x_sample, // fast version
                                const double & x_pop_total,
                                const NumericVector & g_weights,
                                int n_boots) {
  const int n = x_sample.size();

  // NumericVector gw = clone(g_weights); // clone since we probably modify it
  // check if all g_weights are non-negative
  double sg = 0;
  for(int i = 0; i < n; i++) {
    if (g_weights[i] < 0) {
      // sugar for: throw std::range_error();
      stop("g_weights must be a vector of nonnegative numbers.\n");
    }
    sg += g_weights[i];
  }
  // scaling is done below, when computing w
  // if (std::abs(sg - n) > 1e-6) {
  //   // stop("g_weights should sum up to the sample size.\n");
  //   double corection = n / sg;
  //   for(int i = 0; i < n; i++) {
  //     gw[i] = gw[i] * corection;
  //   }
  // }
  // check if scaling works:
  // sg = 0;
  // for(int i = 0; i < n; i++) sg += gw[i];
  // Rcpp::Rcout << (n - sg) << "\n";

  // we implement the formula (9), page 168, from [Zadlo 2021]
  IntegerMatrix M (n_boots, n);
  NumericVector w (n);
  // it was: tx_over_n = x_pop_total / n, but we scale it by: n / sg,
  const double tx_over_n = x_pop_total / sg;
  for (int k = 0; k < n; k++) {
    w[k] = tx_over_n * g_weights[k] / x_sample[k];
  }
  for(int row = 0; row < n_boots; row++) {
    NumericVector u = Rcpp::runif(n, 0, 1);
    for (int i = 0; i < n; i++) {
      // compute ctp's and their sum
      double S = 0;
      NumericVector ctp (n);
      for (int j = 0; j < n; j++) {
        ctp[j] = w[j] - M(row, j);
        ctp[j] = (ctp[j] > 0) ? (ctp[j] * x_sample[j]) : 0;
        S += ctp[j];
      }
      // sample one unit based on u, S and ctp
      double tresh = S * u[i];
      int k = 0;
      double sn = ctp[k];
      while(sn < tresh) {
        sn += ctp[++k];
      }
      M(row, k) += 1;
    }
  }
  return M;
}
