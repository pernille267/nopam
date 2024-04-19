#' Checks Whether a Smoothing Window is Valid - P
#'
#' @param window A \code{data.table} object containing the relevant smoothing data for a particular window.
#' @param measure A \code{character} value referring to the measure we want to smooth. Available choices include \code{median}, \code{hyper} and \code{hypo}. Nothihng else is accepted.
#' @param dur A \code{integer} value. Signifying the number of days the slope must exceed the threshold value to return a warning.
#'
#' @return A \code{logical} value. If the requirements are satisfied, returns \code{TRUE}.
#' @export
#'
#' @examples print(1)

valid_slope_window_p <- function(window, measure, dur = 3L){
  if(measure == "median"){
    minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$Median)) >= dur + 2)
    if(!minimum_requirement){
      return(minimum_requirement)
    }

    last_days <- order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]
    n_valid_last_days <- sum(!is.na(window$Median[last_days]))
    has_enough_last_values <- n_valid_last_days == dur + 2

    if(has_enough_last_values){
      return(has_enough_last_values)
    }
    else{
      minimum_requirement <- sum(!is.na(window$Median)) >= dur + 4
      if(!minimum_requirement){
        return(minimum_requirement)
      }
      last_days_extended <- order(window$`Measured At`, decreasing = TRUE)[1:(dur + 4)]
      last_days_extended_weekend_days <- window$`Is Weekend`[last_days_extended]
      last_days_extended_na_days <- is.na(window$Median[last_days_extended])
      last_days_extended_na_weekend_days <- last_days_extended_weekend_days & last_days_extended_na_days
      weekend_buffer <- sum(last_days_extended_na_weekend_days)
      has_enough_last_values <- sum(!last_days_extended_na_days) + weekend_buffer == dur + 4
    }

    return(minimum_requirement & has_enough_last_values)
  }
  else if(measure == "hyper"){
    minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hyper Percentage`)) >= dur + 2)
    if(!minimum_requirement){
      return(minimum_requirement)
    }

    last_days <- order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]
    n_valid_last_days <- sum(!is.na(window$`Hyper Percentage`[last_days]))
    has_enough_last_values <- n_valid_last_days == dur + 2

    if(has_enough_last_values){
      return(has_enough_last_values)
    }
    else{
      minimum_requirement <- sum(!is.na(window$`Hyper Percentage`)) >= dur + 4
      if(!minimum_requirement){
        return(minimum_requirement)
      }
      last_days_extended <- order(window$`Measured At`, decreasing = TRUE)[1:(dur + 4)]
      last_days_extended_weekend_days <- window$`Is Weekend`[last_days_extended]
      last_days_extended_na_days <- is.na(window$`Hyper Percentage`[last_days_extended])
      last_days_extended_na_weekend_days <- last_days_extended_weekend_days & last_days_extended_na_days
      weekend_buffer <- sum(last_days_extended_na_weekend_days)
      has_enough_last_values <- sum(!last_days_extended_na_days) + weekend_buffer == dur + 4
    }

    return(minimum_requirement & has_enough_last_values)
  }
  else if(measure == "hypo"){
    minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hypo Percentage`)) >= dur + 2)
    if(!minimum_requirement){
      return(minimum_requirement)
    }

    last_days <- order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]
    n_valid_last_days <- sum(!is.na(window$`Hypo Percentage`[last_days]))
    has_enough_last_values <- n_valid_last_days == dur + 2

    if(has_enough_last_values){
      return(has_enough_last_values)
    }
    else{
      minimum_requirement <- sum(!is.na(window$`Hypo Percentage`)) >= dur + 4
      if(!minimum_requirement){
        return(minimum_requirement)
      }
      last_days_extended <- order(window$`Measured At`, decreasing = TRUE)[1:(dur + 4)]
      last_days_extended_weekend_days <- window$`Is Weekend`[last_days_extended]
      last_days_extended_na_days <- is.na(window$`Hypo Percentage`[last_days_extended])
      last_days_extended_na_weekend_days <- last_days_extended_weekend_days & last_days_extended_na_days
      weekend_buffer <- sum(last_days_extended_na_weekend_days)
      has_enough_last_values <- sum(!last_days_extended_na_days) + weekend_buffer == dur + 4
    }

    return(minimum_requirement & has_enough_last_values)
  }
}
#' Checks Whether a Smoothing Window is Valid - X
#'
#' @param window A \code{data.table} object containing the relevant smoothing data for a particular window.
#' @param measure A \code{character} value referring to the measure we want to smooth. Available choices include \code{median}, \code{hyper} and \code{hypo}. Nothihng else is accepted.
#' @param dur A \code{integer} value. Signifying the number of days the slope must exceed the threshold value to return a warning.
#'
#' @return A \code{logical} value. If the requirements are satisfied, returns \code{TRUE}.
#' @export
#'
#' @examples print(1)

valid_slope_window_x <- function(window, measure, dur = 3L){
  if(measure == "median"){
    minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$Median)) >= dur + 2)
    if(!minimum_requirement){
      return(minimum_requirement)
    }

    last_group <- get_gaps(x = window$Median, get_grouping = FALSE, get_last_group = TRUE)
    last_group <- window[last_group, ]
    has_enough_last_values <- sum(!is.na(last_group$Median)) >= dur + 2

    return(minimum_requirement & has_enough_last_values)
  }
  else if(measure == "hyper"){
    minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hyper Percentage`)) >= dur + 2)
    if(!minimum_requirement){
      return(minimum_requirement)
    }

    last_group <- get_gaps(x = window$`Hyper Percentage`, get_grouping = FALSE, get_last_group = TRUE)
    last_group <- window[last_group, ]
    has_enough_last_values <- sum(!is.na(last_group$`Hyper Percentage`)) >= dur + 2

    return(minimum_requirement & has_enough_last_values)
  }
  else if(measure == "hypo"){
    minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hypo Percentage`)) >= dur + 2)
    if(!minimum_requirement){
      return(minimum_requirement)
    }

    last_group <- get_gaps(x = window$`Hypo Percentage`, get_grouping = FALSE, get_last_group = TRUE)
    last_group <- window[last_group, ]
    has_enough_last_values <- sum(!is.na(last_group$`Hypo Percentage`)) >= dur + 2

    return(minimum_requirement & has_enough_last_values)
  }
}
