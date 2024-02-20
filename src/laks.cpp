#include <Rcpp.h>
using namespace Rcpp;

// Gaussian kernel function
double gaussian_kernel(double x1, double x2, double b) {
  return (1.0 / sqrt(2.0 * M_PI * b)) * exp(-0.5 * pow(x1 - x2, 2.0) / pow(b, 2.0));
}

// Calculate gradients using second differences
NumericVector laks_gradients(NumericVector date, NumericVector smoothed_median){
  int n = smoothed_median.size();
  NumericVector gradients(n);
  for(int i = 0; i < n; ++i){
    if(i + 2 <= n - 1){
      gradients[i + 1] = (smoothed_median[i + 2] - smoothed_median[i]) / 2.0;
      //gradients[i + 1] = (smoothed_median[i + 2] - smoothed_median[i]) / (date[i + 2] - date[i]);
      gradients[i + 1] = std::atan(gradients[i + 1]) * (180.0 / M_PI);
    }
  }
  return gradients;
}


//' Local Average Kernel Smoothing Method Using A Gaussian Kernel Function
//'
//' @title Local Average Kernel Smoothing Method Using A Gaussian Kernel Function
//' @name laks
//'
//' @param date A \code{numeric} vector containing dates for the smoothing.
//' @param median A \code{numeric} vector contaning measurements for the relevant dates. Can contain NA-values, which will be removed if present.
//' @param bandwidth A \code{double}. The bandwidth used in the kernel smoothing
//' @param average A \code{double} representing the sample mean of relevant measurements. If set to below \code{0}, it will be calculated in this function instead.
//' @param standard_deviation A \code{double} representing the sample standard deviation of relevant measurements. If set to below \code{0}, it will be calculated in this function instead.
//' @param approximate A \code{logical} value. If set to \code{TRUE}, gradients are estimated via second differences. If set to \code{FALSE}, a plug-in estimator of the gradients is used.
//' @param tol A \code{double} value. Minimum relative weight for observation to be included in calculation
//'
//' @description Local average kernel smoothing of a set of measure values
//'
//' @details Measurement values with weight < 0.01 is not included for efficiency
//'
//' @return A \code{list} containing kernel smoothed values.
//'
//' @examples \dontrun{
//'
//'   print(1)
//'
//' }

// [[Rcpp::export]]
List laks(NumericVector date, NumericVector median, double bandwidth, double average = -0.00001, double standard_deviation = -0.00001, bool approximate = false, double tol = 0.01){

  // Obtain the number of valid values (i.e., not NA-values)
  int number_of_valid = 0;
  for(int i = 0; i < median.size(); ++i){
    bool median_i_is_na = ISNAN(median[i]);
    if(!median_i_is_na){
      ++number_of_valid;
    }
  }

  // Remove NA-values from the 'median' vector.
  int n = number_of_valid;
  int counter = 0;
  NumericVector cleansed_median(n);
  NumericVector cleansed_date(n);
  for(int i = 0; i < median.size(); ++i){
    bool median_i_is_na = ISNAN(median[i]);
    if(!median_i_is_na){
      cleansed_median[counter] = median[i];
      cleansed_date[counter] = date[i];
      ++counter;
    }
  }

  // Calculate 'average' & 'standard deviation' if they are triggered
  if(average <= -0.00001){
    average = mean(cleansed_median);
  }
  if(standard_deviation <= -0.00001){
    standard_deviation = sd(cleansed_median);
  }

  // Vectors to be filled in for loop
  NumericVector smoothed_median(n);
  NumericVector normalized_smoothed_median(n);
  NumericVector gradients(n);
  NumericVector normalized_median = (cleansed_median - average) / standard_deviation;


  if(!approximate){
    double delta = 1e-10;  // small perturbation for numerical gradient
    for (int i = 0; i < n; ++i) {
      double weighted_sum_normalized = 0.0;
      double weighted_sum = 0.0;
      double weight_total = 0.0;
      double weighted_sum_delta = 0.0;
      double weight_total_delta = 0.0;
      double self_weight = gaussian_kernel(cleansed_date[i], cleansed_date[i], bandwidth);
      double self_weight_delta = gaussian_kernel(cleansed_date[i], cleansed_date[i] + delta, bandwidth);

      weighted_sum_normalized += self_weight * normalized_median[i];
      weighted_sum_delta += self_weight_delta * normalized_median[i];
      weighted_sum += self_weight * cleansed_median[i];
      weight_total += self_weight;
      weight_total_delta += self_weight_delta;

      // Lower points
      for(int j = i - 1; j > -1; --j){
        double weight = gaussian_kernel(cleansed_date[j], cleansed_date[i], bandwidth);
        double weight_delta = gaussian_kernel(cleansed_date[j], cleansed_date[i] + delta, bandwidth);

        if(weight / self_weight < tol){
          break;
        }

        weighted_sum_normalized += weight * normalized_median[j];
        weighted_sum_delta += weight_delta * normalized_median[j];
        weighted_sum += weight * cleansed_median[j];
        weight_total += weight;
        weight_total_delta += weight_delta;
      }

      // Upper points
      for(int j = i + 1; j < n; ++j){
        double weight = gaussian_kernel(cleansed_date[j], cleansed_date[i], bandwidth);
        double weight_delta = gaussian_kernel(cleansed_date[j], cleansed_date[i] + delta, bandwidth);

        if(weight / self_weight < tol){
          break;
        }

        weighted_sum_normalized += weight * normalized_median[j];
        weighted_sum_delta += weight_delta * normalized_median[j];
        weighted_sum += weight * cleansed_median[j];
        weight_total += weight;
        weight_total_delta += weight_delta;
      }

      normalized_smoothed_median[i] = weighted_sum_normalized / weight_total;
      smoothed_median[i] = weighted_sum / weight_total;
      double smoothed_median_delta = weighted_sum_delta / weight_total_delta;

      // Numerical gradient estimation
      gradients[i] = (smoothed_median_delta - normalized_smoothed_median[i]) / delta;
      gradients[i] = std::atan(gradients[i]) * (180.0 / M_PI);
    }
  }
  else{
    for (int i = 0; i < n; ++i) {
      double weighted_sum_normalized = 0.0;
      double weighted_sum = 0.0;
      double weight_total = 0.0;
      double self_weight = gaussian_kernel(cleansed_date[i], cleansed_date[i], bandwidth);

      weighted_sum_normalized += self_weight * normalized_median[i];
      weighted_sum += self_weight * cleansed_median[i];
      weight_total += self_weight;

      // Lower points
      for(int j = i - 1; j > -1; --j){
        double weight = gaussian_kernel(cleansed_date[j], cleansed_date[i], bandwidth);

        if(weight / self_weight < tol){
          break;
        }

        weighted_sum_normalized += weight * normalized_median[j];
        weighted_sum += weight * cleansed_median[j];
        weight_total += weight;
      }

      // Upper points
      for(int j = i + 1; j < n; ++j){
        double weight = gaussian_kernel(cleansed_date[j], cleansed_date[i], bandwidth);

        if(weight / self_weight < tol){
          break;
        }

        weighted_sum_normalized += weight * normalized_median[j];
        weighted_sum += weight * cleansed_median[j];
        weight_total += weight;
      }
      normalized_smoothed_median[i] = weighted_sum_normalized / weight_total;
      smoothed_median[i] = weighted_sum / weight_total;
    }
    gradients = laks_gradients(cleansed_date, normalized_smoothed_median);
  }
  //Rcout << "normalized smoothed: " << normalized_smoothed_median << "\n";
  return List::create(Named("date") = cleansed_date, Named("smoothed_median") = smoothed_median, Named("gradients") = gradients);
}

//' Local Average Kernel Smoothing Leave One Out Cross Validation
//'
//' @title Local Average Kernel Smoothing Leave One Out Cross Validation
//' @name laks_loo
//'
//' @param date A \code{numeric} vector containing dates for the smoothing.
//' @param median A \code{numeric} vector contaning measurements for the relevant dates. Cannot contain NA-values.
//' @param bandwidth A \code{double}. The bandwidth used in the kernel smoothing
//'
//' @description Local average kernel smoothing loo cv of a set of measure values
//'
//' @details Measurement values with weight < 0.01 is not included for efficiency
//'
//' @return A \code{double}, which is the cross-validation value for the fit given the \code{bandwidth}.
//'
//' @examples \dontrun{
//'
//'   print(1)
//'
//' }

// [[Rcpp::export]]
double laks_loo(NumericVector date, NumericVector median, double bandwidth) {
  int n = median.size();
  NumericVector smoothed_median(n);
  NumericVector gradients(n);
  NumericVector cv_vector(n);
  for (int i = 0; i < n; ++i) {
    double weighted_sum = 0.0;
    double weight_total = 0.0;
    double i_weight = 0;
    double i_median = 0;
    for (int j = 0; j < n; ++j) {
      double self_weight = gaussian_kernel(date[i], date[i], bandwidth);
      double weight = gaussian_kernel(date[j], date[i], bandwidth);
      // Only consider weights >= 1%
      if (weight/self_weight >= 0.01){
        if(i==j){
          i_weight += weight;
          i_median += median[j];
        }
        weighted_sum += weight * median[j];
        weight_total += weight;
      }
    }
    smoothed_median[i] = weighted_sum / weight_total;
    cv_vector[i] = pow((i_median - smoothed_median[i]) / (1.0 - i_weight / weight_total), 2);
  }
  double mean_cv = mean(cv_vector);
  return mean_cv;
}
