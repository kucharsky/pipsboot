#include <Rcpp.h>
using namespace Rcpp;

//' Zadlo \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS bootstrap
//'
//' The function generates bootstrap samples mimicking Brewer's \ifelse{html}{\out{&pi;}}{\eqn{\pi}}PS
//'   sampling scheme (Tille (2006), 112-114) for pseudo-population where numbers
//'   of replications of sample elements do not have to be integers and the
//'   pseudo-population is not physically constructed.
//'
//' @param n_replications Numeric vector. Positive values (not necessarily
//' integers) of number of replications of sample elements in the pseudo-populations
//' e.g. inverses of first order inclusion probabilities for sample elements or
//' calibration weights for sample elements.
//' @param pi_sample Numeric vector. Values of first order inclusion probabilities for sample elements.
//' @param n_boots Integer. A number of bootstrap samples to be generated.
//'
//' @return Integer matrix with `n_boots` rows and `length(x_sample)` columns.
//' Each row is a bootstrap sample, and the number in i-th column denotes how many
//' times a given unit is resampled.
//'
//' @references Tille, Y. (2006) Sampling algorithms, Springer, New York.
//'
//' @examples
//' N <- 200 # population size
//' x_population <- rgamma(N, 1)
//' n <- 10  # sample size
//' pi <- inclusion_probabilities_rcpp(x_population, n)
//' set.seed(123)
//' sample <- UPbrewer_rcpp(pi)
//' x_s <- x_population[sample == 1]
//' pi_s <- pi[sample == 1]
//' # g-weigths
//' g <- sampling::calib(matrix(c(x_s, rep(1, n)), ncol = 2),
//'                      d = 1 / pi_s, c(sum(x_population), N), method = "linear")
//' # numbers of replications of sample elements in the pseudo-population
//' # equal inverses of first order inclusion probabilities:
//' boot_Z_pips_rcpp(1 / pi_s, pi_s, 5)
//' # numbers of replications of sample elements in the pseudo-population equal to calibration weights:
//' boot_Z_pips_rcpp(g / pi_s, pi_s, 5)
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_Z_pips_rcpp(const NumericVector & n_replications,
                               const NumericVector & pi_sample,
                               int n_boots) {
  int n = pi_sample.size();  // population size
  IntegerMatrix sb (n_boots, n);      // sample: 1 - in the sample, 0 - out of the sample
  for(int i = 0; i < n; i++) {
    // Check for NaN
    if (R_IsNaN(pi_sample[i])) {
      Rcpp::Rcout << "NaN value in inclusion probabilities vector.\n";
      return IntegerMatrix (1, 1);
    }
    // Check for NA
    if (R_IsNA(pi_sample[i])) {
      Rcpp::Rcout << "NA value in inclusion probabilities vector.\n";
      return IntegerMatrix (1, 1);
    }
    // Check for pik > 1
    if(pi_sample[i] > 1) {
      Rcpp::Rcout << "Inclusion probability greater than 1: pi_sample[" << i + 1 << "] = " << pi_sample[i] << "\n";
      return IntegerMatrix (1, 1);
    }
    // Check for pik < 0
    if(pi_sample[i] < 0) {
      Rcpp::Rcout << "Negative inclusion probability: pi_sample[" << i + 1 << "] = " << pi_sample[i] << "\n";
      return IntegerMatrix (1, 1);
    }
  }
  // no negative n_replications allowed
  NumericVector nrep = clone(n_replications);
  double sum_n_rep = 0;
  for(int j = 0; j < n; j++){
    if(nrep[j] < 0) {
      nrep[j] = 0;
    }
    sum_n_rep += nrep[j];
  }
  if(sum_n_rep == 0) {
    Rcpp::Rcout << "No positive n_replications.\n";
    return IntegerMatrix (1, 1);
  }
  NumericVector p = no_init(n);
  for(int row = 0; row < n_boots; row++) {
    double nma = n; // n - a, in the beginning a == 0
    NumericVector n_rep_min_sb = clone(nrep); // n_replications - sum(sb(row, _))
    NumericVector u = Rcpp::runif(n, 0, 1);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        double part_b = nma - pi_sample[j];
        double denom = (nma - pi_sample[j] * (n - i));
        p[j] = n_rep_min_sb[j] * pi_sample[j] * part_b / denom;
      }
      for(int j = 1; j < n; j++) {
        p[j] = p[j - 1] + p[j];     // p = cumsum(p);
      }
      u[i] = u[i] * p[n - 1];
      int j = 0;
      while(u[i] > p[j]) {
        j++;
      }
//      Rcout << i << " " << j << " " << u[i] << " " << p[j] << "\n";
      sb(row, j)++;
      nma -= pi_sample[j];
      n_rep_min_sb[j]--;
      if(n_rep_min_sb[j] < 0) {
        n_rep_min_sb[j] = 0;
      }
    }
  }
  return sb;
}
