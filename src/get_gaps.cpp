#include <Rcpp.h>
using namespace Rcpp;

//' Get Gaps for Smoothing Window
//'
//' @title Get Gaps for Smoothing Window
//' @name get_gaps
//'
//' @param x A \code{numeric} vector containing either medians, hyper percentages or hypo percentages used for the smoothing window.
//' @param get_grouping A non-missing \code{logical} value. If \code{TRUE}, the function returns an \code{integer} vector that describes the grouping of the x-values based on the gaps observed. Overrides \code{get_last_group}.
//' @param get_last_group A non-missing \code{logical} value. If \code{TRUE}, the function returns an \code{integer} vector containing the indicies of the last group based on the gaps observed in the x-values. \code{get_grouping} must be set to \code{FALSE} in order for this to work.
//'
//' @description Efficiently provide information about the gaps in the measurements used for kernel smoothing. Used to validate whether we have enough information to give a slope warning to the user.
//'
//' @details Cannot locate gaps in x if dates with missing values are stripped. If both \code{get_grouping} and \code{get_last_group} is \code{FALSE}, only the gap counts are returned.
//'
//' @return A \code{integer} vector containing relevant gap information based on \code{get_grouping} and \code{get_last_group} arguments.
//'
//' @examples \dontrun{
//'
//'   print(1)
//'
//' }

// [[Rcpp::export]]
IntegerVector get_gaps(NumericVector x, bool get_grouping = false, bool get_last_group = false){
  int n = x.size();
  int count = 0;
  int group = 1;
  IntegerVector gaps(n);
  IntegerVector grouping(n);
  for(int i = 0; i < n; ++i){
    bool x_i_is_na = ISNAN(x[i]);
    if(x_i_is_na){
      count++;
      gaps[i] = count;
    }
    else{
      count = 0;
      gaps[i] = 0;
    }
  }
  grouping[0] = group;
  for(int i = 1; i < n; ++i){
    if(gaps[i - 1] > 3 and gaps[i] == 0){
      group++;
      grouping[i] = group;
    }
    else if(gaps[i - 1] > 3 and gaps[i] > 3){
      grouping[i] = 0;
    }
    else if(gaps[i - 1] == 3 and gaps[i] > 3){
      grouping[i] = 0;
    }
    else{
      grouping[i] = group;
    }
  }

  int m = 0;
  for(int i = 0; i < n; ++i){
    if(grouping[i] == group){
      m++;
    }
  }

  IntegerVector last_group(m);
  m = 0;
  for(int i = 0; i < n; ++i){
    if(grouping[i] == group){
      last_group[m] = i + 1;
      m++;
    }
  }

  if(get_grouping){
    return grouping;
  }

  else if(get_last_group){
    return last_group;
  }

  return gaps;
}
