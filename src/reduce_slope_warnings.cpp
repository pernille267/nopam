#include <Rcpp.h>
using namespace Rcpp;

//' Reduce Slope Type Warnings
//'
//' @title Reduce Slope Type Warnings
//' @name reduce_slope_warnings
//'
//' @param data A \code{list} containing slope warning data.
//' @param stringent A \code{logical} value. Should slope warnings be reduced?
//'
//' @description Reduce slope type warnings using the slope functionality
//'
//' @details No details.
//'
//' @return A \code{list} with updated warnings
//'
//' @examples \dontrun{
//'
//'   print(1)
//'
//' }

// Evaluate slope warnings
// [[Rcpp::export]]
List reduce_slope_warnings(List data, bool stringent = true){
  DateVector date = data["Measured At"];
  NumericVector mean_slope = data["Mean Slope Magnitude"];
  NumericVector evaluated_average = data["Evaluated Average"];
  NumericVector evaluated_standard_deviation = data["Evaluated Standard Deviation"];
  LogicalVector condition_met = data["Condition Met"];
  CharacterVector evaluated_days = data["Evaluated Days"];
  CharacterVector evaluated_slopes = data["Evaluated Slopes"];
  if(!stringent){
    LogicalVector is_warn = condition_met;
    return List::create(Named("Evaluated Days") = evaluated_days,
                        Named("Evaluated Slopes") = evaluated_slopes,
                        Named("Evaluated Average") = evaluated_average,
                        Named("Evaluated Standard Deviation") = evaluated_standard_deviation,
                        Named("Measured At") = date,
                        Named("Mean Slope Magnitude") = mean_slope,
                        Named("Condition Met") = condition_met,
                        Named("Is Warning") = is_warn);
  }
  bool can_give_warning = true;
  int last_warning_dir = 0;
  int current_dir = 0;
  int n = condition_met.size();
  LogicalVector is_warn(n);
  is_warn[0] = false;
  for(int i = 1; i < n; ++i){
    if(can_give_warning & condition_met[i]){
      is_warn[i] = true;
      can_give_warning = false; // avoids that warnings are delivered each day
      if(mean_slope[i] < 0){
        last_warning_dir = -1;
      }
      else{
        last_warning_dir = 1;
      }
    }
    else if((!can_give_warning) & condition_met[i]){
      if(mean_slope[i] < 0){
        current_dir = -1;
      }
      else if(mean_slope[i] > 0){
        current_dir = 1;
      }
      else{
        current_dir = 0;
      }
      if(current_dir == last_warning_dir){
        is_warn[i] = false;
      }
      // If the slope has changed signs, we deliver a new warning even though can_give_warning is false.
      else if(current_dir != last_warning_dir){
        is_warn[i] = true;
      }
    }
    else if(can_give_warning & !condition_met[i]){
      is_warn[i] = false;
    }
    // If the slope magnitude is less than the upper limit, open avenues for new warnings
    else if(!can_give_warning & !condition_met[i]){
      is_warn[i] = false;
      can_give_warning = true;
    }
  }
  return List::create(Named("Evaluated Days") = evaluated_days,
                      Named("Evaluated Slopes") = evaluated_slopes,
                      Named("Evaluated Average") = evaluated_average,
                      Named("Evaluated Standard Deviation") = evaluated_standard_deviation,
                      Named("Measured At") = date,
                      Named("Mean Slope Magnitude") = mean_slope,
                      Named("Condition Met") = condition_met,
                      Named("Is Warning") = is_warn);
}
