#include <Rcpp.h>
using namespace Rcpp;

#include "halfdoubled_rcpp.h"
#include "UPpoisson_rcpp.h"
#include "UPbrewer_rcpp.h"
#include "inclusion_probabilities_rcpp.h"

// version 4 - vectorize inside function
// it's faster when non_selected_num == 1)

//' Antal and Tille (2014) bootstrap for complex sampling designs
//'
//' This function implements resampling algorithm described in Antal and Tille
//' (2014), pp. 1355-1356 for sampling with unequal probabilities without replacement.
//'
//' @param pi_sample Numeric vector of inclusion probabilities of the first kind
//' (of elements in the original sample).
//' @param n_boots Integer. Number of bootstrap samples to create.
//'
//' @return Integer matrix with `n_boots` rows and `length(phik)` columns.
//' Each row is a bootstrap sample, with entries in columns denoting how many
//' times given unit is resampled.
//'
//' @references Erika Antal and Yves Tille (2014) A new resampling method for
//' sampling designs without replacement: the doubled half bootstrap.
//' Comput Stat 29, 1345-1363, https://doi.org/10.1007/s00180-014-0495-0.
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix boot_AT2014_rcpp(const NumericVector & pi_sample, int n_boots) {
  int n = pi_sample.size();
  IntegerMatrix S (n_boots, n);
  NumericVector pik_odds = (1 - pi_sample) / pi_sample;
  pik_odds = 1 - inclusion_probabilities_rcpp(pik_odds / sum(pik_odds), 2);
  for (int row = 0; row < n_boots; row++) {
    IntegerVector s = UPpoisson_rcpp(pi_sample);
    int non_selected_num = n;
    for (int j = 0; j < n; j++) {
      S(row, j) = s[j];
      non_selected_num -= s[j];
    }
    if (non_selected_num >= 2) {
      IntegerVector v = halfdoubled_rcpp(non_selected_num);
      int i = 0, k = 0;
      while (k < non_selected_num) {
        if (S(row, i) == 0) {
          S(row, i) = v[k++];
        }
        i++;
      }
    } else {
      if (non_selected_num == 1) {
        if (R::runif(0, 1) < 0.5) {
          for (int j = 0; j < n; j++) {
            S(row, j) = 1;
          }
        } else {
          IntegerVector b = UPbrewer_rcpp(pik_odds);
          for (int j = 0; j < n; j++) {
            S(row, j) = b[j];
          }
          IntegerVector v = halfdoubled_rcpp(2);
          int i = 0, k = 0;
          while (k < 2) {
            if (S(row, i) == 0) {
              S(row, i) = v[k++];
            }
            i++;
          }
        }
      }
    }
  }
  return S;
}
