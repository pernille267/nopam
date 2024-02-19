#include <Rcpp.h>
using namespace Rcpp;

//' Moving Median Smoothing Method
//'
//' @title Moving Median Smoothing Method
//' @name moving_median
//'
//' @param date A \code{numeric} vector containing dates for the smoothing.
//' @param median A \code{numeric} vector contaning measurements for the relevant dates. Can contain NA-values, which will be removed if present.
//' @param bandwidth A \code{double}. The bandwidth used in the kernel smoothing
//'
//' @description Moving median smoothing of a set of measure values
//'
//' @details No details.
//'
//' @return A \code{list} containing moving median smoothed values.
//'
//' @examples \dontrun{
//'
//'   print(1)
//'
//' }

// [[Rcpp::export]]
List moving_median(NumericVector date, NumericVector median, int bandwidth = 30) {
  int n = date.size();
  NumericVector moving_median(n);
  IntegerVector number_of_results(n);
  for (int i = 0; i < n; i++) {
    int current_date = date[i];
    int start_date = current_date - bandwidth + 1;

    // Collecting observations within the date window
    std::vector<double> window;
    for (int j = 0; j <= i; j++) {
      if (date[j] >= start_date && date[j] <= current_date) {
        bool median_j_is_na = ISNAN(median[j]);
        if(!median_j_is_na){
          window.push_back(median[j]);
        }
      }
    }

    // Calculate the median of the window
    std::sort(window.begin(), window.end());
    int window_length = window.size();
    number_of_results[i] = window_length;

    if(window_length == 0){
      moving_median[i] = NA_REAL; // Assign NA if there are no observations in the window
    }
    else if(window_length % 2 == 0){
      moving_median[i] = (window[window_length / 2 - 1] + window[window_length / 2]) / 2.0;
    }
    else{
      moving_median[i] = window[window_length / 2];
    }
  }

  return List::create(Named("date") = date, Named("moving_median") = moving_median, Named("number_of_results") = number_of_results);
}


