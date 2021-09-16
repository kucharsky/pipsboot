#include <Rcpp.h>
using namespace Rcpp;

#include <math.h>

//' Brewer sampling
//'
//' Uses the Brewer's method to select a sample of units (unequal probabilities,
//' without replacement, fixed sample size).
//' We follow the Algorithm 6.10 from Tille (2006), pp. 112-114, however details
//' of implementation differs a little from that of the \code{sampling} package.
//'
//' @param pik Numeric vector of  the first order inclusion probabilities.
//'
//' @return Numeric vector of the same length as \code{pik}. Value of 1
//'         indicates that unit is selected into sample, 0 means that unit is
//'         not in the sample.
//'
//' @references Tille, Y. (2006) Sampling algorithms, Springer, New York.
//'
//' @examples
//' UPbrewer_rcpp((0:5) / 5)
//' UPbrewer_rcpp(c(10:0, 5) / 10)
//'
//' @export
// [[Rcpp::export]]
IntegerVector UPbrewer_rcpp(const NumericVector & pik) {
  NumericVector pikc = clone(pik); // we need a copy of pik to be modified
  int N = pik.size();      // population size
  IntegerVector sample(N); // sample: 1 - in the sample, 0 - out of the sample
  double nsum = 0;         // sum of inclusion probabilities
  for(int i = 0; i < N; i++) {
    if (R_IsNaN(pik[i])) { // Check for NaN
      stop("NaN value in inclusion probabilities vector.\n");
    }
    if (R_IsNA(pik[i])) {  // Check for NA
      stop("NA value in inclusion probabilities vector.\n");
    }
    if(pik[i] > 1) {       // Check for pik > 1
      stop("Inclusion probability greater than 1: pik[ %u ] = %f\n", i+1, pik[i]);
    }
    if(pik[i] < 0) {       // Check for pik < 0
      stop("Negative inclusion probability: pik[ %u ] = %f\n", i+1, pik[i]);
    }
    if(pik[i] == 1) {
      // element is for sure in the sample,
      // but has to be removed from the main part of the algorithm
      sample[i] = 1;
      // we need to have pikc[j] == 0 for already sampled units
      // for reasons explained inside the loop
      // since we modify pikc, we need to do this on the copy!
      pikc[i] = 0;
    }
    nsum += pikc[i]; // sum up to determine size of the probabilistic part of the sample
  }
  // Check whether nsum is (close to) an integer
  if (abs(nsum - round(nsum)) > 1e-6) {
    stop("Inclusion probabilities do not sum to an integer (diff = %f).\n",
         abs(nsum - round(nsum)));
  }
  int n = (int) round(nsum); // sample size
  if (n == 0) { // nothing to do
    return sample;
  }
  // here we declare something equal to "n - a" in sampling package,
  // that is "n minus sum of sampled piks"
  // it will be modified later; nsum has to stay fixed
  double nma = nsum;
  NumericVector p = no_init(N); // for calculations
  NumericVector u = Rcpp::runif(n, 0, 1); // for sampling
  for (int i = 0; i < n; i++) {
    double sp = 0;
    double nsumi = nsum - i; // not "i - 1", since index starts from 0
    for (int j = 0; j < N; j++) {
      // we need to have pikc[j] == 0 for already sampled units
      // to avoid 0/0 = nan values below
      p[j] = (1 - sample[j]) * pikc[j] * (nma - pikc[j]) / (nma - pikc[j] * nsumi);
      sp += p[j]; // sum of p[j] - for normalization
      // for debugging
      // Rcout << "i = " << i << "\tj = " << j << "\tp[j] = " << p[j] << "\tsp = " << sp << "\n";
    }
    // we look for the first j such that sum_{k=1}^j p_k > u[i] * sum_{k=1}^N p_k
    // its equivalent to sampling with probabilities p / sum(p)
    u[i] = u[i] * sp; // p sums to sp
    double acc = p[0];
    int k = 0;
    while (acc <= u[i]) {
      acc += p[++k];
    }
    sample[k] = 1;
    nma -= pikc[k];
    // Rcpp::checkUserInterrupt();
  }
  return sample;
}
