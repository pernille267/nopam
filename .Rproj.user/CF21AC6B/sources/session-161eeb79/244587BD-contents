#include <Rcpp.h>
using namespace Rcpp;

//' Reduce Bias Type Warnings
//'
//' @title Reduce Bias Type Warnings
//' @name reduce_bias_warnings
//'
//' @param warnings A \code{logical} vector containing raw bias type warnings.
//' @param snooze A \code{integer}. How many days should be used in snoozing.
//'
//' @description Reduce bias type warnings using the snoozing functionality
//'
//' @details No details.
//'
//' @return A \code{logical} vector containing a reduced set of warnings after snoozing.
//'
//' @examples \dontrun{
//'
//'   print(1)
//'
//' }

// [[Rcpp::export]]
LogicalVector reduce_bias_warnings(LogicalVector warnings, int snooze = 5){
  if(snooze == 0){
    return warnings;
  }
  bool can_give_warning = true;
  int snooze_tally = 0;
  int n = warnings.size();
  LogicalVector reduced_warnings(n);
  for(int i = 0; i < n; ++i){
    if(can_give_warning){
      if(warnings[i]){
        reduced_warnings[i] = true;
        can_give_warning = false;
        snooze_tally = 1;
      }
    }
    else{
      if(snooze_tally >= snooze){
        can_give_warning = true;
      }
      else if(!warnings[i]){
        can_give_warning = true;
      }
      else{
        ++snooze_tally;
      }
    }
  }
  return reduced_warnings;
}

