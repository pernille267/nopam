#' Alternative Warning Process for One Instrument Within One Laboratory
#'
#' @param data A \code{data.table} containing measurements.
#' @param pg_data A \code{data.table} containing peer group measurements.
#' @param from A \code{IDate} value. The date before the start of the warning process.
#' @param to A \code{IDate} value. The date before the end of the warning process.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param method A \code{character} signifying which smoothing method that is used. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing). Only relevant if \code{warning} is set to \code{slope}.
#' @param warning A \code{character} signifying which warning process that is of interest to simulate.
#' @param tol A \code{double} value. What is the warning threshold for relevant \code{warning}?
#' @param dur A \code{integer} value. How many days must the gradient be above the threshold given in \code{tol} in order to trigger a warning. Only relevant if \code{warning} is set to \code{slope}.
#' @param bw1 A \code{double}. The bandwidth used for the kernel smoothing. Typically set to 11.
#' @param bw2 A \code{double}. The bandwidth used for monthly moving median. Typically set to 30.
#' @param bw3 A \code{double}. The bandwidth used for yearly moving median. Typically set to 365.
#' @param snooze A \code{integer} signifying the number of days recurring warnings are snoozed. Only relevant if \code{warning} is \code{bias} or \code{peer_group}. Set to \code{0} to disable snoozing.
#' @param stringent_slope_warning A \code{logical} value. If set to \code{FALSE}, all triggered warnings are displayed. If set to \code{TRUE}, we avoid excessive number of warnings. See \code{?reduce_slope_warnings()} for more information. Only relevant if \code{warning} is set to \code{slope}.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#' @param approximate A \code{logical} value. If set to \code{TRUE}, second differences are used to estimate gradients. If \code{FALSE} other methods are utilized. Only relevant if \code{warning} is set to \code{slope}.
#'
#' @return A \code{data.table} object containing the warning process information
#' @export
#'
#' @examples print(1)
xwarning_process0 <- function(data, pg_data = NULL, from = "2021-01-01", to = "2022-12-31", measure = c("median", "hyper", "hypo"), method = c("lc", "ll"), warning = c("slope", "bias", "peer_group"), tol = 4, dur = 3, bw1 = 11, bw2 = 30, bw3 = 365, snooze = 0L, stringent_slope_warning = FALSE, attach = TRUE, approximate = TRUE){

  `Measured At` <- `m` <- `Evaluated Slopes` <- `Manufacturer Name` <- `Instrument Model Name` <- NULL
  measure <- measure[1]
  warning <- warning[1]
  method <- method[1]

  if(warning == "slope"){
    slope_padding <- ceiling(2 * sqrt(log(10)) * bw1)
    slope_start_date <- as.IDate(from) - slope_padding - (dur + 2) + 1
    slope_warning_lower_end_window <- as.IDate(from)
    slope_warning_upper_end_window <- as.IDate(to)
    slope_warning_lower_start_window <- slope_start_date
    slope_warning_upper_start_window <- as.IDate(to) - slope_padding - (dur + 2) + 1
    slope_warning_start_window <- seq(from = slope_warning_lower_start_window, to = slope_warning_upper_start_window, by = "1 day")
    slope_warning_end_window <- seq(from = slope_warning_lower_end_window, to = slope_warning_upper_end_window, by = "1 day")

    if(length(slope_warning_start_window) != length(slope_warning_end_window)){
      stop("Window is not always of equal length...")
    }

    current_average <- sapply(1:length(slope_warning_start_window), FUN = function(i) mean(data[`Measured At` <= slope_warning_end_window[i]]$`Median`, na.rm = TRUE), simplify = TRUE)
    current_sd <- sapply(1:length(slope_warning_start_window), FUN = function(i) sd(data[`Measured At` <= slope_warning_end_window[i]]$`Median`, na.rm = TRUE), simplify = TRUE)
    sub_data <- sapply(1:length(slope_warning_start_window), FUN = function(i) data[`Measured At` >= slope_warning_start_window[i] & `Measured At` <= slope_warning_end_window[i]], simplify = FALSE)

    valid_slope_window <- function(window, measure, dur = 3L){
      if(measure == "median"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$Median)) >= dur + 2)
        if(!minimum_requirement){
          return(minimum_requirement)
        }
        weekend_days <- window$`Is Weekend`[order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]]
        safe_guard <- min(dur + 2 + 1, nrow(window))
        weekend_day_extra <- window$`Is Weekend`[order(window$`Measured At`, decreasing = TRUE)[safe_guard]]
        na_days <- is.na(window$Median[order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]])
        na_day_extra <- is.na(window$Median[order(window$`Measured At`, decreasing = TRUE)[safe_guard]]) & (nrow(window) >= dur + 2 + 1)
        na_weekend_days <- weekend_days & na_days
        number_of_na_weekend_days <- sum(na_weekend_days)
        if(na_days[dur + 2] & na_day_extra){
          number_of_na_weekend_days <- number_of_na_weekend_days + (na_day_extra & weekend_day_extra)
        }
        safe_guard <- min(dur + 2 + number_of_na_weekend_days, nrow(window))
        has_enough_last_values <- sum(!is.na(window$Median[order(window$`Measured At`, decreasing = TRUE)[1:safe_guard]])) == dur + 2

        return(minimum_requirement & has_enough_last_values)
      }
      else if(measure == "hyper"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hyper Percentage`)) >= dur + 2)
        if(!minimum_requirement){
          return(minimum_requirement)
        }
        weekend_days <- window$`Is Weekend`[order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]]
        safe_guard <- min(dur + 2 + 1, nrow(window))
        weekend_day_extra <- window$`Is Weekend`[order(window$`Measured At`, decreasing = TRUE)[safe_guard]]
        na_days <- is.na(window$`Hyper Percentage`[order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]])
        na_day_extra <- is.na(window$`Hyper Percentage`[order(window$`Measured At`, decreasing = TRUE)[safe_guard]]) & (nrow(window) >= dur + 2 + 1)
        na_weekend_days <- weekend_days & na_days
        number_of_na_weekend_days <- sum(na_weekend_days)
        if(na_days[dur + 2] & na_day_extra){
          number_of_na_weekend_days <- number_of_na_weekend_days + (na_day_extra & weekend_day_extra)
        }
        safe_guard <- min(dur + 2 + number_of_na_weekend_days, nrow(window))
        has_enough_last_values <- sum(!is.na(window$`Hyper Percentage`[order(window$`Measured At`, decreasing = TRUE)[1:safe_guard]])) == dur + 2
        return(minimum_requirement & has_enough_last_values)
      }
      else if(measure == "hypo"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hypo Percentage`)) >= dur + 2)
        if(!minimum_requirement){
          return(minimum_requirement)
        }
        weekend_days <- window$`Is Weekend`[order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]]
        safe_guard <- min(dur + 2 + 1, nrow(window))
        weekend_day_extra <- window$`Is Weekend`[order(window$`Measured At`, decreasing = TRUE)[safe_guard]]
        na_days <- is.na(window$`Hypo Percentage`[order(window$`Measured At`, decreasing = TRUE)[1:(dur + 2)]])
        na_day_extra <- is.na(window$`Hypo Percentage`[order(window$`Measured At`, decreasing = TRUE)[safe_guard]]) & (nrow(window) >= dur + 2 + 1)
        na_weekend_days <- weekend_days & na_days
        number_of_na_weekend_days <- sum(na_weekend_days)
        if(na_days[dur + 2] & na_day_extra){
          number_of_na_weekend_days <- number_of_na_weekend_days + (na_day_extra & weekend_day_extra)
        }
        safe_guard <- min(dur + 2 + number_of_na_weekend_days, nrow(window))
        has_enough_last_values <- sum(!is.na(window$`Hypo Percentage`[order(window$`Measured At`, decreasing = TRUE)[1:safe_guard]])) == dur + 2
        return(minimum_requirement & has_enough_last_values)
      }
    }

    valid_windows <- sapply(X = sub_data, FUN = valid_slope_window, measure, dur, simplify = TRUE)

    output <- lapply(1:length(valid_windows), FUN = function(i){
      if(measure == "median"){
        if(valid_windows[i]){
          if(is.na(current_average[i])){
            current_average[i] <- 0
          }
          if(is.na(current_sd[i])){
            current_sd[i] <- 1
          }
          current_kernel_smoothed <- kernel_smoothing0(data = sub_data[[i]], method = method, measure = measure, bw = bw1, average = current_average[i], standard_deviation = current_sd[i], attach = "x", approximate = approximate, matrix_approach = TRUE, diagnostics = TRUE, na_rm = FALSE)
          current_kernel_smoothed_cleansed <- current_kernel_smoothed[!is.na(m), ]
          current_date <- current_kernel_smoothed$`Measured At`[which.max(current_kernel_smoothed$`Measured At`)]
          evaluated_ids <- order(current_kernel_smoothed_cleansed$`Measured At`, decreasing = TRUE)[2:(dur + 1)]
          evaluated_slopes <- current_kernel_smoothed_cleansed$`Smoothed Gradient Degrees`[evaluated_ids]
          evaluated_days <- current_kernel_smoothed_cleansed$`Measured At`[evaluated_ids]
          all_above_limit <- all(abs(evaluated_slopes) > tol)
          all_same_sign <- all(sign(evaluated_slopes) == -1) | all(sign(evaluated_slopes) == 1)
          warning_data <- data.table("Measured At" = current_date,
                                     "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                     "Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                     "Evaluated Average" = current_average[i],
                                     "Evaluated Standard Deviation" = current_sd[i],
                                     "Mean Slope Magnitude" = round(mean(evaluated_slopes, na.rm = TRUE), 3L),
                                     "Condition Met" = all_above_limit & all_same_sign)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Evaluated Slopes" = NA_character_,
                                     "Evaluated Days" = NA_character_,
                                     "Evaluated Average" = NA_real_,
                                     "Evaluated Standard Deviation" = NA_real_,
                                     "Mean Slope Magnitude" = NA_real_,
                                     "Condition Met" = FALSE)
          return(warning_data)
        }
      }
      else if(measure == "hyper" | measure == "hypo"){
        if(valid_windows[i]){
          current_kernel_smoothed <- kernel_smoothing0(data = sub_data[[i]], method = method, measure = measure, bw = bw1, average = 0, standard_deviation = 1, attach = "x", approximate = approximate, matrix_approach = TRUE, diagnostics = TRUE, na_rm = FALSE)
          current_kernel_smoothed_cleansed <- current_kernel_smoothed[!is.na(m), ]
          current_date <- current_kernel_smoothed$`Measured At`[which.max(current_kernel_smoothed$`Measured At`)]
          evaluated_ids <- order(current_kernel_smoothed_cleansed$`Measured At`, decreasing = TRUE)[2:(dur + 1)]
          evaluated_slopes <- current_kernel_smoothed_cleansed$`Smoothed Gradient Degrees`[evaluated_ids]
          evaluated_days <- current_kernel_smoothed_cleansed$`Measured At`[evaluated_ids]
          all_above_limit <- all(abs(evaluated_slopes) > tol)
          all_same_sign <- all(sign(evaluated_slopes) == -1) | all(sign(evaluated_slopes) == 1)
          warning_data <- data.table("Measured At" = current_date,
                                     "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                     "Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                     "Evaluated Average" = 0,
                                     "Evaluated Standard Deviation" = 1,
                                     "Mean Slope Magnitude" = round(mean(evaluated_slopes, na.rm = TRUE), 3L),
                                     "Condition Met" = all_above_limit & all_same_sign)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Evaluated Slopes" = NA_character_,
                                     "Evaluated Days" = NA_character_,
                                     "Evaluated Average" = NA_real_,
                                     "Evaluated Standard Deviation" = NA_real_,
                                     "Mean Slope Magnitude" = NA_real_,
                                     "Condition Met" = FALSE)
          return(warning_data)
        }
      }
    })

    output <- rbindlist(output)
    output <- output[!is.na(`Evaluated Slopes`)]
    output <- reduce_slope_warnings(output, stringent = stringent_slope_warning) |> setDT()

    return(output)
  }

  else if(warning == "bias"){
    bias_month_start_date <- as.IDate(from) - bw2 + 1
    bias_year_start_date <- as.IDate(from) - bw3 + 1
    bias_month_warning_lower_end_window <- as.IDate(from)
    bias_year_warning_lower_end_window <- as.IDate(from)
    bias_month_warning_upper_end_window <- as.IDate(to)
    bias_year_warning_upper_end_window <- as.IDate(to)
    bias_month_warning_lower_start_window <- bias_month_start_date
    bias_year_warning_lower_start_window <- bias_year_start_date
    bias_month_warning_upper_start_window <- as.IDate(to) - bw2 + 1
    bias_year_warning_upper_start_window <- as.IDate(to) - bw3 + 1
    bias_month_warning_start_window <- seq(from = bias_month_warning_lower_start_window, to = bias_month_warning_upper_start_window, by = "1 day")
    bias_month_warning_end_window <- seq(from = bias_month_warning_lower_end_window, to = bias_month_warning_upper_end_window, by = "1 day")
    bias_year_warning_start_window <- seq(from = bias_year_warning_lower_start_window, to = bias_year_warning_upper_start_window, by = "1 day")
    bias_year_warning_end_window <- seq(from = bias_year_warning_lower_end_window, to = bias_year_warning_upper_end_window, by = "1 day")

    common_length <- length(bias_month_warning_end_window)

    if(all(bias_month_warning_start_window == common_length, bias_month_warning_end_window == common_length, bias_year_warning_start_window == common_length, bias_year_warning_end_window == common_length)){
      stop("Window is not always of equal length...")
    }

    sub_data_month <- sapply(1:common_length, FUN = function(i) data[`Measured At` >= bias_month_warning_start_window[i] & `Measured At` <= bias_month_warning_end_window[i]], simplify = FALSE)
    sub_data_year <- sapply(1:common_length, FUN = function(i) data[`Measured At` >= bias_year_warning_start_window[i] & `Measured At` <= bias_year_warning_end_window[i]], simplify = FALSE)

    valid_bias_window <- function(window, measure){
      if(measure == "median"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$Median)) >= 5)
        return(minimum_requirement)
      }
      else if(measure == "hyper"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hyper Percentage`)) >= 5)
        return(minimum_requirement)
      }
      else if(measure == "hypo"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hypo Percentage`)) >= 5)
        return(minimum_requirement)
      }
    }

    valid_windows_month <- sapply(X = sub_data_month, FUN = valid_bias_window, measure, simplify = TRUE)
    valid_windows_year <- sapply(X = sub_data_year, FUN = valid_bias_window, measure, simplify = TRUE)

    output <- lapply(1:length(valid_windows_month), FUN = function(i){
      if(measure == "median"){
        if(valid_windows_month[i] & valid_windows_year[i]){
          current_date <- sub_data_month[[i]]$`Measured At`[which.max(sub_data_month[[i]]$`Measured At`)]
          month_start_date <- sub_data_month[[i]]$`Measured At`[which.min(sub_data_month[[i]]$`Measured At`)]
          current_month_median <- median(sub_data_month[[i]]$Median, na.rm = TRUE)
          current_month_results <- sum(!is.na(sub_data_month[[i]]$Median))
          year_start_date <- sub_data_year[[i]]$`Measured At`[which.min(sub_data_year[[i]]$`Measured At`)]
          current_year_median <- median(sub_data_year[[i]]$Median, na.rm = TRUE)
          current_year_results <- sum(!is.na(sub_data_year[[i]]$Median))
          current_bias <- abs((current_month_median - current_year_median) / current_year_median) * 100
          current_warning <- current_bias >= tol
          warning_data <- data.table("Measured At" = current_date,
                                     "Monthly Start Date" = month_start_date,
                                     "Monthly Median" = current_month_median,
                                     "Monthly Number of Results" = current_month_results,
                                     "Yearly Median" = current_year_median,
                                     "Yearly Start Date" = year_start_date,
                                     "Yearly Number of Results" = current_year_results,
                                     "Bias Percentage" = round(current_bias, 3L),
                                     "Is Warning" = current_warning)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Monthly Start Date" = as.IDate("2000-01-01"),
                                     "Monthly Median" = NA_real_,
                                     "Monthly Number of Results" = NA_integer_,
                                     "Yearly Median" = NA_real_,
                                     "Yearly Start Date" = as.IDate("2000-01-01"),
                                     "Yearly Number of Results" = NA_integer_,
                                     "Bias Percentage" = NA_real_,
                                     "Is Warning" = FALSE)
          return(warning_data)
        }
      }
      else if(measure == "hyper"){
        if(valid_windows_month[i] & valid_windows_year[i]){

          current_date <- sub_data_month[[i]]$`Measured At`[which.max(sub_data_month[[i]]$`Measured At`)]
          month_start_date <- sub_data_month[[i]]$`Measured At`[which.min(sub_data_month[[i]]$`Measured At`)]
          current_month_median <- median(sub_data_month[[i]]$`Hyper Percentage`, na.rm = TRUE)
          current_month_results <- sum(!is.na(sub_data_month[[i]]$`Hyper Percentage`))
          year_start_date <- sub_data_year[[i]]$`Measured At`[which.min(sub_data_year[[i]]$`Measured At`)]
          current_year_median <- median(sub_data_year[[i]]$`Hyper Percentage`, na.rm = TRUE)
          current_year_results <- sum(!is.na(sub_data_year[[i]]$`Hyper Percentage`))
          current_bias <- abs((current_month_median - current_year_median) / current_year_median) * 100
          current_warning <- current_bias >= tol
          warning_data <- data.table("Measured At" = current_date,
                                     "Monthly Start Date" = month_start_date,
                                     "Monthly Median" = current_month_median,
                                     "Monthly Number of Results" = current_month_results,
                                     "Yearly Median" = current_year_median,
                                     "Yearly Start Date" = year_start_date,
                                     "Yearly Number of Results" = current_year_results,
                                     "Bias Percentage" = round(current_bias, 3L),
                                     "Is Warning" = current_warning)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Monthly Start Date" = as.IDate("2000-01-01"),
                                     "Monthly Median" = NA_real_,
                                     "Monthly Number of Results" = NA_integer_,
                                     "Yearly Median" = NA_real_,
                                     "Yearly Start Date" = as.IDate("2000-01-01"),
                                     "Yearly Number of Results" = NA_integer_,
                                     "Bias Percentage" = NA_real_,
                                     "Is Warning" = FALSE)
          return(warning_data)
        }
      }
      else if(measure == "hypo"){
        if(valid_windows_month[i] & valid_windows_year[i]){

          current_date <- sub_data_month[[i]]$`Measured At`[which.max(sub_data_month[[i]]$`Measured At`)]
          month_start_date <- sub_data_month[[i]]$`Measured At`[which.min(sub_data_month[[i]]$`Measured At`)]
          current_month_median <- median(sub_data_month[[i]]$`Hypo Percentage`, na.rm = TRUE)
          current_month_results <- sum(!is.na(sub_data_month[[i]]$`Hypo Percentage`))
          year_start_date <- sub_data_year[[i]]$`Measured At`[which.min(sub_data_year[[i]]$`Measured At`)]
          current_year_median <- median(sub_data_year[[i]]$`Hypo Percentage`, na.rm = TRUE)
          current_year_results <- sum(!is.na(sub_data_year[[i]]$`Hypo Percentage`))
          current_bias <- abs((current_month_median - current_year_median) / current_year_median) * 100
          current_warning <- current_bias >= tol
          warning_data <- data.table("Measured At" = current_date,
                                     "Monthly Start Date" = month_start_date,
                                     "Monthly Median" = current_month_median,
                                     "Monthly Number of Results" = current_month_results,
                                     "Yearly Median" = current_year_median,
                                     "Yearly Start Date" = year_start_date,
                                     "Yearly Number of Results" = current_year_results,
                                     "Bias Percentage" = round(current_bias, 3L),
                                     "Is Warning" = current_warning)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Monthly Start Date" = as.IDate("2000-01-01"),
                                     "Monthly Median" = NA_real_,
                                     "Monthly Number of Results" = NA_integer_,
                                     "Yearly Median" = NA_real_,
                                     "Yearly Start Date" = as.IDate("2000-01-01"),
                                     "Yearly Number of Results" = NA_integer_,
                                     "Bias Percentage" = NA_real_,
                                     "Is Warning" = FALSE)
          return(warning_data)
        }
      }
    })

    output <- rbindlist(output)
    output$`Is Warning` <- reduce_bias_warnings(warnings = output$`Is Warning`, snooze = snooze)
    return(output)
  }

  else if(warning == "peer_group"){
    bias_month_start_date <- as.IDate(from) - bw2 + 1
    bias_month_warning_lower_end_window <- as.IDate(from)
    bias_month_warning_upper_end_window <- as.IDate(to)
    bias_month_warning_lower_start_window <- bias_month_start_date
    bias_month_warning_upper_start_window <- as.IDate(to) - bw2 + 1
    bias_month_warning_start_window <- seq(from = bias_month_warning_lower_start_window, to = bias_month_warning_upper_start_window, by = "1 day")
    bias_month_warning_end_window <- seq(from = bias_month_warning_lower_end_window, to = bias_month_warning_upper_end_window, by = "1 day")

    common_length <- length(bias_month_warning_end_window)

    if(all(bias_month_warning_start_window == common_length, bias_month_warning_end_window == common_length)){
      stop("Window is not always of equal length...")
    }

    # Extract peer group data from the full dataset!
    pg_data <- pg_data[`Manufacturer Name` %in% data$`Manufacturer Name` & `Instrument Model Name` %in% data$`Instrument Model Name`]

    if(nrow(pg_data) < nrow(data)){
      stop("Something went wrong when extracting peer group data...")
    }

    sub_data_month <- sapply(1:common_length, FUN = function(i) data[`Measured At` >= bias_month_warning_start_window[i] & `Measured At` <= bias_month_warning_end_window[i]], simplify = FALSE)
    sub_data_month_peer_group <- sapply(1:common_length, FUN = function(i) pg_data[`Measured At` >= bias_month_warning_start_window[i] & `Measured At` <= bias_month_warning_end_window[i]], simplify = FALSE)

    valid_peer_group_window <- function(window, measure){
      if(measure == "median"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$Median)) >= 5)
        return(minimum_requirement)
      }
      else if(measure == "hyper"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hyper Percentage`)) >= 5)
        return(minimum_requirement)
      }
      else if(measure == "hypo"){
        minimum_requirement <- (nrow(window) >= 1) & (sum(!is.na(window$`Hypo Percentage`)) >= 5)
        return(minimum_requirement)
      }
    }

    valid_windows_month <- sapply(X = sub_data_month, FUN = valid_peer_group_window, measure, simplify = TRUE)
    valid_windows_month_peer_group <- sapply(X = sub_data_month_peer_group, FUN = valid_peer_group_window, measure, simplify = TRUE)



    output <- lapply(1:length(valid_windows_month), FUN = function(i){
      if(measure == "median"){
        if(valid_windows_month[i] & valid_windows_month_peer_group[i]){
          current_date <- sub_data_month[[i]]$`Measured At`[which.max(sub_data_month[[i]]$`Measured At`)]
          month_start_date <- sub_data_month[[i]]$`Measured At`[which.min(sub_data_month[[i]]$`Measured At`)]
          current_month_median <- median(sub_data_month[[i]]$Median, na.rm = TRUE)
          current_month_results <- sum(!is.na(sub_data_month[[i]]$Median))
          peer_group_start_date <- sub_data_month_peer_group[[i]]$`Measured At`[which.min(sub_data_month_peer_group[[i]]$`Measured At`)][1]
          current_peer_group_median <- median(sub_data_month_peer_group[[i]]$Median, na.rm = TRUE)
          current_peer_group_results <- sum(!is.na(sub_data_month_peer_group[[i]]$Median))
          current_bias <- abs((current_month_median - current_peer_group_median) / current_peer_group_median) * 100
          current_warning <- current_bias >= tol
          warning_data <- data.table("Measured At" = current_date,
                                     "Monthly Start Date" = month_start_date,
                                     "Monthly Median" = current_month_median,
                                     "Monthly Number of Results" = current_month_results,
                                     "Peer Group Monthly Median" = current_peer_group_median,
                                     "Peer Group Start Date" = peer_group_start_date,
                                     "Peer Group Monthly Number of Results" = current_peer_group_results,
                                     "Bias Percentage" = round(current_bias, 3L),
                                     "Is Warning" = current_warning)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Monthly Start Date" = as.IDate("2000-01-01"),
                                     "Monthly Median" = NA_real_,
                                     "Monthly Number of Results" = NA_integer_,
                                     "Peer Group Monthly Median" = NA_real_,
                                     "Peer Group Start Date" = as.IDate("2000-01-01"),
                                     "Peer Group Monthly Number of Results" = NA_integer_,
                                     "Bias Percentage" = NA_real_,
                                     "Is Warning" = FALSE)
          return(warning_data)
        }
      }
      else if(measure == "hyper"){
        if(valid_windows_month[i] & valid_windows_month_peer_group[i]){

          current_date <- sub_data_month[[i]]$`Measured At`[which.max(sub_data_month[[i]]$`Measured At`)]
          month_start_date <- sub_data_month[[i]]$`Measured At`[which.min(sub_data_month[[i]]$`Measured At`)]
          current_month_median <- median(sub_data_month[[i]]$`Hyper Percentage`, na.rm = TRUE)
          current_month_results <- sum(!is.na(sub_data_month[[i]]$`Hyper Percentage`))
          peer_group_start_date <- sub_data_month_peer_group[[i]]$`Measured At`[which.min(sub_data_month_peer_group[[i]]$`Measured At`)][1]
          current_peer_group_median <- median(sub_data_month_peer_group[[i]]$`Hyper Percentage`, na.rm = TRUE)
          current_peer_group_results <- sum(!is.na(sub_data_month_peer_group[[i]]$`Hyper Percentage`))
          current_bias <- abs((current_month_median - current_peer_group_median) / current_peer_group_median) * 100
          current_warning <- current_bias >= tol
          warning_data <- data.table("Measured At" = current_date,
                                     "Monthly Start Date" = month_start_date,
                                     "Monthly Median" = current_month_median,
                                     "Monthly Number of Results" = current_month_results,
                                     "Peer Group Monthly Median" = current_peer_group_median,
                                     "Peer Group Start Date" = peer_group_start_date,
                                     "Peer Group Monthly Number of Results" = current_peer_group_results,
                                     "Bias Percentage" = round(current_bias, 3L),
                                     "Is Warning" = current_warning)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Monthly Start Date" = as.IDate("2000-01-01"),
                                     "Monthly Median" = NA_real_,
                                     "Monthly Number of Results" = NA_integer_,
                                     "Peer Group Monthly Median" = NA_real_,
                                     "Peer Group Start Date" = as.IDate("2000-01-01"),
                                     "Peer Group Monthly Number of Results" = NA_integer_,
                                     "Bias Percentage" = NA_real_,
                                     "Is Warning" = FALSE)
          return(warning_data)
        }
      }
      else if(measure == "hypo"){
        if(valid_windows_month[i] & valid_windows_month_peer_group[i]){

          current_date <- sub_data_month[[i]]$`Measured At`[which.max(sub_data_month[[i]]$`Measured At`)]
          month_start_date <- sub_data_month[[i]]$`Measured At`[which.min(sub_data_month[[i]]$`Measured At`)]
          current_month_median <- median(sub_data_month[[i]]$`Hypo Percentage`, na.rm = TRUE)
          current_month_results <- sum(!is.na(sub_data_month[[i]]$`Hypo Percentage`))
          peer_group_start_date <- sub_data_month_peer_group[[i]]$`Measured At`[which.min(sub_data_month_peer_group[[i]]$`Measured At`)][1]
          current_peer_group_median <- median(sub_data_month_peer_group[[i]]$`Hypo Percentage`, na.rm = TRUE)
          current_peer_group_results <- sum(!is.na(sub_data_month_peer_group[[i]]$`Hypo Percentage`))
          current_bias <- abs((current_month_median - current_peer_group_median) / current_peer_group_median) * 100
          current_warning <- current_bias >= tol
          warning_data <- data.table("Measured At" = current_date,
                                     "Monthly Start Date" = month_start_date,
                                     "Monthly Median" = current_month_median,
                                     "Monthly Number of Results" = current_month_results,
                                     "Peer Group Monthly Median" = current_peer_group_median,
                                     "Peer Group Start Date" = peer_group_start_date,
                                     "Peer Group Monthly Number of Results" = current_peer_group_results,
                                     "Bias Percentage" = round(current_bias, 3L),
                                     "Is Warning" = current_warning)
          return(warning_data)
        }
        else{
          warning_data <- data.table("Measured At" = as.IDate("2000-01-01"),
                                     "Monthly Start Date" = as.IDate("2000-01-01"),
                                     "Monthly Median" = NA_real_,
                                     "Monthly Number of Results" = NA_integer_,
                                     "Peer Group Monthly Median" = NA_real_,
                                     "Peer Group Start Date" = as.IDate("2000-01-01"),
                                     "Peer Group Monthly Number of Results" = NA_integer_,
                                     "Bias Percentage" = NA_real_,
                                     "Is Warning" = FALSE)
          return(warning_data)
        }
      }
    })
    output <- rbindlist(output)
    output$`Is Warning` <- reduce_bias_warnings(warnings = output$`Is Warning`, snooze = snooze)
    return(output)

  }

}

#' Alternative Warning Process for Multiple Instruments Within Multiple Laboratories
#'
#' @param data A \code{data.table} containing measurements.
#' @param by A \code{character} vector containing the names of the grouping variables.
#' @param from A \code{IDate} value. The date before the start of the warning process.
#' @param to A \code{IDate} value. The date before the end of the warning process.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param method A \code{character} signifying which smoothing method that is used. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing). Only relevant if \code{warning} is set to \code{slope}.
#' @param warning A \code{character} signifying which warning process that is of interest to simulate.
#' @param tol A \code{double} value. What is the warning threshold for relevant \code{warning}?
#' @param dur A \code{integer} value. How many days must the gradient be above the threshold given in \code{tol} in order to trigger a warning. Only relevant if \code{warning} is set to \code{slope}.
#' @param bw1 A \code{double}. The bandwidth used for the kernel smoothing. Typically set to 11.
#' @param bw2 A \code{double}. The bandwidth used for monthly moving median. Typically set to 30.
#' @param bw3 A \code{double}. The bandwidth used for yearly moving median. Typically set to 365.
#' @param snooze A \code{integer} signifying the number of days recurring warnings are snoozed. Only relevant if \code{warning} is \code{bias} or \code{peer_group}. Set to \code{0} to disable snoozing.
#' @param stringent_slope_warning A \code{logical} value. If set to \code{FALSE}, all triggered warnings are displayed. If set to \code{TRUE}, we avoid excessive number of warnings. See \code{?reduce_slope_warnings()} for more information. Only relevant if \code{warning} is set to \code{slope}.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#' @param approximate A \code{logical} value. If set to \code{TRUE}, second differences are used to estimate gradients. If \code{FALSE} other methods are utilized. Only relevant if \code{warning} is set to \code{slope}.
#' @param only_warnings A \code{logical} value. If set to \code{TRUE}, only the dates with warnings generated are outputted. If set to \code{FALSE}, the whole process it outputted.
#'
#' @return A \code{data.table} object containing the warning process information
#' @export
#'
#' @examples print(1)
xwarning_process <- function(data, by, from = "2021-12-31", to = "2022-12-30", measure = c("median", "hyper", "hypo"), method = c("lc", "ll"), warning = c("slope", "bias", "peer_group"), tol = 4, dur = 3, bw1 = 11, bw2 = 30, bw3 = 365, snooze = 0L, stringent_slope_warning = FALSE, attach = TRUE, approximate = TRUE, only_warnings = TRUE){

  `Is Warning` <- NULL

  measure <- measure[1]
  method <- method[1]
  warning <- warning[1]
  pg_data <- data
  output <- data[, xwarning_process0(data = .SD, pg_data = pg_data, from = from, to = to, measure = measure, method = method, warning = warning, tol = tol, dur = dur, bw1 = bw1, bw2 = bw2, bw3 = bw3, snooze = snooze, stringent_slope_warning = stringent_slope_warning, attach = attach, approximate = approximate),
                 by = by]
  if(only_warnings){
    return(output[`Is Warning` == TRUE])
  }
  else{
    return(output)
  }
}

#' Alternative Generation of All Warnings for All Instruments Within Each Laboratory
#'
#' @param data A \code{data.table} containing measurements.
#' @param from A \code{IDate} value. The start of the warning process, i.e., the first date where a warning is calculated.
#' @param to A \code{IDate} value. The end of the warning process, i.e., the last date where a warning is calculated.
#' @param method A \code{character} signifying which smoothing method that is used for the slope warnings. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing).
#' @param tol_median_slope A \code{double} value. What is the warning threshold for the median slope warning?
#' @param tol_hyper_slope A \code{double} value. What is the warning threshold for the hyper slope warning?
#' @param tol_hypo_slope A \code{double} value. What is the warning threshold for the hypo slope warning?
#' @param tol_median_bias A \code{double} value. What is the warning threshold for the median month-year bias warning?
#' @param tol_hyper_bias A \code{double} value. What is the warning threshold for the hyper month-year bias warning?
#' @param tol_hypo_bias A \code{double} value. What is the warning threshold for the hypo month-year bias warning?
#' @param tol_pg  A \code{double} value. What is the warning threshold for the median peer group bias warning?
#' @param dur A \code{integer} value. How many days must the gradient be above the threshold given in \code{tol} in order to trigger a slope warning. Will be the same for all measures for the slope warning.
#' @param bw1 A \code{double}. The bandwidth used for the kernel smoothing. Typically set to 11.
#' @param bw2 A \code{double}. The bandwidth used for monthly moving median. Typically set to 30.
#' @param bw3 A \code{double}. The bandwidth used for yearly moving median. Typically set to 365.
#' @param only_warnings A \code{logical} value. If set to \code{TRUE}, only the dates with warnings generated are outputted. If set to \code{FALSE}, the whole process it outputted.
#'
#' @return A \code{data.table} object containing the warning process information
#' @export
#'
#' @examples print(1)
xgenerate_warnings <- function(data, from = "2022-01-02", to = "2023-01-01", method = c("lc", "ll"), tol_median_slope = 2, tol_hyper_slope = 8, tol_hypo_slope = 8, tol_median_bias = 10, tol_hyper_bias = 10, tol_hypo_bias = 10, tol_pg = 10, dur = 3, bw1 = 11, bw2 = 30, bw3 = 365, only_warnings = FALSE){

  `Measured At` <- `Analyte Name` <- `Laboratory Code` <- `Instrument Code` <- `Warning Type` <- NULL
  method <- method[1]
  actual_from <- as.IDate(from) - 1
  actual_to <- as.IDate(to) - 1

  slope_median_warnings <- xwarning_process(data = data, by = c("Analyte Name", "Laboratory Code", "Instrument Code"), from = actual_from, to = actual_to, measure = "median", method = method, tol = tol_median_slope, warning = "slope", dur = dur, bw1 = bw1, stringent_slope_warning = FALSE, attach = FALSE, approximate = TRUE, only_warnings = only_warnings)
  slope_hyper_warnings <- xwarning_process(data = data, by = c("Analyte Name", "Laboratory Code", "Instrument Code"), from = actual_from, to = actual_to, measure = "hyper", method = method, tol = tol_hyper_slope, warning = "slope", dur = dur, bw1 = bw1, stringent_slope_warning = FALSE, attach = FALSE, approximate = TRUE, only_warnings = only_warnings)
  slope_hypo_warnings <- xwarning_process(data = data, by = c("Analyte Name", "Laboratory Code", "Instrument Code"), from = actual_from, to = actual_to, measure = "hypo", method = method, tol = tol_hypo_slope, warning = "slope", dur = dur, bw1 = bw1, stringent_slope_warning = FALSE, attach = FALSE, approximate = TRUE, only_warnings = only_warnings)

  slope_median_warnings$`Warning Type` <- "SMOOTHED_MEDIAN_DEVIATES_FOR_DAYS"
  slope_median_warnings$`Threshold Value` <- tol_median_slope
  slope_hyper_warnings$`Warning Type` <- "SMOOTHED_HYPER_DEVIATES_FOR_DAYS"
  slope_hyper_warnings$`Threshold Value` <- tol_hyper_slope
  slope_hypo_warnings$`Warning Type` <- "SMOOTHED_HYPO_DEVIATES_FOR_DAYS"
  slope_hypo_warnings$`Threshold Value` <- tol_hypo_slope

  slope_warnings <- rbind(slope_median_warnings,
                          slope_hyper_warnings,
                          slope_hypo_warnings, fill = TRUE)
  slope_warnings$`Warning At` <- slope_warnings$`Measured At` + 1

  setcolorder(slope_warnings, c("Warning At", "Measured At", "Analyte Name", "Instrument Code", "Laboratory Code", "Warning Type", "Threshold Value", "Evaluated Slopes", "Evaluated Days", "Evaluated Average", "Evaluated Standard Deviation", "Condition Met", "Is Warning"))

  slope_warnings$`Condition Met` <- NULL

  bias_median_warnings <- xwarning_process(data = data, by = c("Analyte Name", "Laboratory Code", "Instrument Code"), from = actual_from, to = actual_to, measure = "median", method = method, tol = tol_median_bias, warning = "bias", bw2 = bw2, bw3 = bw3, snooze = 0L, attach = FALSE, only_warnings = only_warnings)
  bias_hyper_warnings <- xwarning_process(data = data, by = c("Analyte Name", "Laboratory Code", "Instrument Code"), from = actual_from, to = actual_to, measure = "hyper", method = method, tol = tol_hyper_bias, warning = "bias", bw2 = bw2, bw3 = bw3, snooze = 0L, attach = FALSE, only_warnings = only_warnings)
  bias_hypo_warnings <- xwarning_process(data = data, by = c("Analyte Name", "Laboratory Code", "Instrument Code"), from = actual_from, to = actual_to, measure = "hypo", method = method, tol = tol_hypo_bias, warning = "bias", bw2 = bw2, bw3 = bw3, snooze = 0L, attach = FALSE, only_warnings = only_warnings)

  bias_median_warnings$`Warning Type` <- "MONTHLY_MEDIAN_DEVIATES_FROM_PREVIOUS_YEAR"
  bias_median_warnings$`Threshold Value` <- tol_median_bias
  bias_hyper_warnings$`Warning Type` <- "MONTHLY_HYPER_DEVIATES_FROM_PREVIOUS_YEAR"
  bias_hyper_warnings$`Threshold Value` <- tol_hyper_bias
  bias_hypo_warnings$`Warning Type` <- "MONTHLY_HYPO_DEVIATES_FROM_PREVIOUS_YEAR"
  bias_hypo_warnings$`Threshold Value` <- tol_hypo_bias

  bias_warnings <- rbind(bias_median_warnings,
                         bias_hyper_warnings,
                         bias_hypo_warnings, fill = TRUE)

  bias_warnings$`Warning At` <- bias_warnings$`Measured At` + 1

  setcolorder(bias_warnings, c("Warning At", "Measured At", "Analyte Name", "Instrument Code", "Laboratory Code", "Warning Type", "Threshold Value", "Monthly Median", "Monthly Start Date", "Monthly Number of Results", "Yearly Median", "Yearly Start Date", "Yearly Number of Results", "Bias Percentage", "Is Warning"))

  peer_group_warnings <- xwarning_process(data = data, by = c("Analyte Name", "Laboratory Code", "Instrument Code"), from = actual_from, to = actual_to, measure = "median", method = method, tol = tol_pg, warning = "peer_group", bw2 = bw2, snooze = 0L, attach = FALSE, only_warnings = only_warnings)

  peer_group_warnings$`Warning Type` <- "MONTHLY_MEDIAN_DEVIATES_FROM_PEER_GROUP"
  peer_group_warnings$`Warning At` <- peer_group_warnings$`Measured At` + 1
  peer_group_warnings$`Threshold Value` <- tol_pg

  setcolorder(peer_group_warnings, c("Warning At", "Measured At", "Analyte Name", "Instrument Code", "Laboratory Code", "Warning Type", "Threshold Value", "Monthly Median", "Monthly Start Date", "Monthly Number of Results", "Peer Group Monthly Median", "Peer Group Start Date", "Peer Group Monthly Number of Results", "Bias Percentage", "Is Warning"))

  slope_warnings$`Measured At` <- as.IDate(slope_warnings$`Measured At`)
  slope_warnings$`Warning At` <- as.IDate(slope_warnings$`Warning At`)

  all_warnings <- rbind(bias_warnings, peer_group_warnings, slope_warnings, fill = TRUE)

  setorder(all_warnings, `Measured At`, `Analyte Name`, `Laboratory Code`, `Instrument Code`, `Warning Type`)

  return(all_warnings)
}


