#include <Rcpp.h>
using namespace Rcpp;

//' Inclusion probabilities
//'
//' Calculates first order inclusion probabilities for a probability
//' proportional-to-size sampling design based on positive values of an
//' auxiliary variable. We use the algorithm described
//' in Tille (2006), pp. 18-19.
//'
//' @param values Numeric vector. Positive values of an auxiliary variable.
//' @param size   Integer. Sample size.
//'
//' @return Vector of size \code{length(values)} of inclusion probabilities.
//'
//' @references Tille, Y. (2006) Sampling algorithms, Springer, New York.
//'
//' @examples
//' inclusion_probabilities_rcpp(1:10, 5)
//' inclusion_probabilities_rcpp(1:10, 6)
//' inclusion_probabilities_rcpp(c(1:8, 50, 100), 5)
//'
//' @export
// [[Rcpp::export]]
NumericVector inclusion_probabilities_rcpp(const NumericVector & values, const int & size) {
  int N = values.size();
  double sum_of_values = 0;
  for(int i = 0; i < N; i++){
    sum_of_values += values[i];
  }
  int ready = 0; // how many probs are fixed (== 1)
  NumericVector probs (N);
  bool there_is_prob_greater_than_one = true;
  while (there_is_prob_greater_than_one) {
    there_is_prob_greater_than_one = false;
    double scale_factor = (size - ready) / sum_of_values;
    sum_of_values = 0;
    for (int i = 0; i < N; i++) {
      if (probs[i] < 1) {
        probs[i] = scale_factor * values[i];
        if (probs[i] < 1) {
          sum_of_values += values[i];
          continue;
        }
        if (probs[i] > 1) {
          there_is_prob_greater_than_one = true;
          ready++;
          probs[i] = 1;
        }
      }
    }
  }
  return probs;
}
