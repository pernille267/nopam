#include <RcppArmadillo.h>

using namespace Rcpp;

// Gaussian kernel function
double gaussian_kernel_1(double x) {
  return 1 / sqrt(2 * M_PI) * exp(-0.5 * x * x);
}

// Calculate S_0, S_1, S_2, ...
double weight_function_r (std::vector<double> x, std::vector<double> w, int r = 0){
  double out = 0.0;
  int m = x.size();
  for(int i = 0; i < m; ++i){
    out += pow(x[i], r) * w[i];
  }
  return out;
}

// Calculate gradients using second differences
NumericVector llks_gradients(NumericVector date, NumericVector smoothed_median){
  int n = smoothed_median.size();
  NumericVector gradients(n);
  for(int i = 0; i < n; ++i){
    if(i + 2 <= n - 1){
      gradients[i + 1] = (smoothed_median[i + 2] - smoothed_median[i]) / 2;
      gradients[i + 1] = std::atan(gradients[i + 1]) * (180.0 / M_PI);
    }
  }
  return gradients;
}

//' Local Linear Kernel Smoothing Method Using A Gaussian Kernel Function
//'
//' @title Local Linear Kernel Smoothing Method Using A Gaussian Kernel Function
//' @name llks
//'
//' @param date A \code{numeric} vector containing dates for the smoothing.
//' @param median A \code{numeric} vector contaning measurements for the relevant dates. Can contain NA-values, which will be removed if present.
//' @param bandwidth A \code{double}. The bandwidth used in the kernel smoothing
//' @param average A \code{double} representing the sample mean of relevant measurements. If set to below \code{0}, it will be calculated in this function instead.
//' @param standard_deviation A \code{double} representing the sample standard deviation of relevant measurements. If set to below \code{0}, it will be calculated in this function instead.
//' @param approximate A \code{logical} value. If set to \code{TRUE}, gradients are estimated via second differences. If set to \code{FALSE}, a plug-in estimator of the gradients is used.
//' @param matrix_approach A \code{loigcal} value. Should matrix calculus be used to smooth. Only relevant if \code{method} is \code{ll}. Matrix calculus is faster, but may not be stable in some cases.
//'
//' @description Local linear kernel smoothing of a set of measure values
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

// Local linear kernel smoothing function
// [[Rcpp::export]]
List llks(NumericVector date, NumericVector median, double bandwidth = 11, double average = -0.00001, double standard_deviation = -0.00001, bool approximate = false, bool matrix_approach = true){

  int number_of_valid = 0;
  for(int i = 0; i < median.size(); ++i){
    bool median_i_is_na = ISNAN(median[i]);
    if(!median_i_is_na){
      ++number_of_valid;
    }
  }

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
  if(n == 0){
    return List::create(Named("date") = NA_REAL, Named("smoothed_median") = NA_REAL, Named("gradients") = NA_REAL);
  }
  else if(n == 1){
    return List::create(Named("date") = cleansed_date, Named("smoothed_median") = cleansed_median, Named("gradients") = 0.0);
  }

  if(average <= -0.00001){
    average = mean(cleansed_median);
  }
  if(standard_deviation <= -0.00001){
    standard_deviation = sd(cleansed_median);
  }
  NumericVector normalized_median = (cleansed_median - average) / standard_deviation;
  NumericVector smoothed_normalized_median(n);
  NumericVector smoothed_median(n);
  NumericVector gradients(n);

  // Use matrix algebra to fit the local-linear model (may fail if X^TWX is computionally singular)
  if(matrix_approach){
    for(int i = 0; i < n; ++i){
      std::vector<double> weights, relevant_dates, relevant_medians, relevant_normalized_medians;
      for (int j = 0; j < n; ++j) {
        double self_weight = gaussian_kernel_1(0.0);
        double weight = gaussian_kernel_1((cleansed_date[j] - cleansed_date[i]) / bandwidth);

        // Only consider weights >= 1%
        if(weight/self_weight >= 0.01){
          weights.push_back(weight);
          relevant_dates.push_back(cleansed_date[j]);
          relevant_medians.push_back(cleansed_median[j]);
          relevant_normalized_medians.push_back(normalized_median[j]);
        }
      }

      // Initialize matrices and vectors for LS fitting
      int m = weights.size();
      arma::mat X(m, 2);
      arma::colvec Y(m);
      arma::colvec Y2(m);
      arma::colvec W(m);

      // Fill in values
      for (int j = 0; j < m; ++j) {
        X(j, 0) = 1.0;
        X(j, 1) = relevant_dates[j];
        Y(j) = relevant_medians[j];
        Y2(j) = relevant_normalized_medians[j];
        W(j) = weights[j];
      }

      // Weighted LS regression fit for the i-th median (based on neighbor points)
      arma::mat Xt = trans(X);
      arma::mat XtWX = Xt * arma::diagmat(W) * X;
      arma::colvec XtWy = Xt * arma::diagmat(W) * Y;
      arma::colvec XtWy2 = Xt * arma::diagmat(W) * Y2;
      arma::colvec beta = arma::solve(XtWX, XtWy);
      arma::colvec beta2 = arma::solve(XtWX, XtWy2);

      // The fitted median
      smoothed_median[i] = beta(0) + beta(1) * cleansed_date[i];

      // Calculate gradient as the slope of the regression line at the i-th median
      if(!approximate){
        gradients[i] = beta2(1);
        gradients[i] = std::atan(gradients[i]) * (180.0 / M_PI);
      }
      // Needed for using second differences approach (not advised here)
      else if(approximate){
        smoothed_normalized_median[i] = beta2(0) + beta2(1) * cleansed_date[i];
      }
    }
  }

  // Use ordinary algebra to fit the local-linear model (slower)
  else{
    for (int i = 0; i < n; ++i) {

      double weighted_sum_beta0_normalized = 0.0;
      double weighted_sum_beta1_normalized = 0.0;
      double weighted_sum_beta0 = 0.0;
      double weighted_sum_beta1 = 0.0;
      double weight_total = 0.0;

      std::vector<double> weights, relevant_dates, relevant_medians, relevant_normalized_medians;
      for (int j = 0; j < n; ++j) {
        double self_weight = gaussian_kernel_1(0.0);
        double weight = gaussian_kernel_1((cleansed_date[j] - cleansed_date[i]) / bandwidth);
        // Only consider weights >= 1%
        if(weight/self_weight >= 0.01){
          weights.push_back(weight);
          relevant_dates.push_back(cleansed_date[j]);
          relevant_medians.push_back(cleansed_median[j]);
          relevant_normalized_medians.push_back(normalized_median[j]);
        }
      }

      int m = weights.size();

      for (int j = 0; j < m; ++j) {
        weighted_sum_beta0_normalized += (weight_function_r(relevant_dates, weights, 2) - weight_function_r(relevant_dates, weights, 1) * relevant_dates[j]) * weights[j] * relevant_normalized_medians[j];
        weighted_sum_beta1_normalized += (weight_function_r(relevant_dates, weights, 0) * relevant_dates[j] - weight_function_r(relevant_dates, weights, 1)) * weights[j] * relevant_normalized_medians[j];
        weighted_sum_beta0 += (weight_function_r(relevant_dates, weights, 2) - weight_function_r(relevant_dates, weights, 1) * relevant_dates[j]) * weights[j] * relevant_medians[j];
        weighted_sum_beta1 += (weight_function_r(relevant_dates, weights, 0) * relevant_dates[j] - weight_function_r(relevant_dates, weights, 1)) * weights[j] * relevant_medians[j];
      }

      weight_total += weight_function_r(relevant_dates, weights, 0) * weight_function_r(relevant_dates, weights, 2) - pow(weight_function_r(relevant_dates, weights, 1), 2);

      double i_beta0_normalized = weighted_sum_beta0_normalized / weight_total;
      double i_beta1_normalized = weighted_sum_beta1_normalized / weight_total;
      double i_beta0 = weighted_sum_beta0 / weight_total;
      double i_beta1 = weighted_sum_beta1 / weight_total;

      smoothed_median[i] = i_beta0 + i_beta1 * cleansed_date[i];

      // Calculate gradient as the slope of the regression line at the i-th median
      if(!approximate){
        gradients[i] = i_beta1_normalized;
        gradients[i] = std::atan(gradients[i]) * (180.0 / M_PI);
      }
      // Needed for using second differences approach (not advised here)
      else if(approximate){
        smoothed_normalized_median[i] = i_beta0_normalized + i_beta1_normalized * cleansed_date[i];
      }
    }
  }

  // Estimate all gradients using second differences approach
  if(approximate){
    gradients = llks_gradients(cleansed_date, smoothed_normalized_median);
  }

  return List::create(Named("date") = cleansed_date, Named("smoothed_median") = smoothed_median, Named("gradients") = gradients);
}
