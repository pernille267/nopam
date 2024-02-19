#' Slope warning process
#'
#' @param data A \code{data.table} containing measurements.
#' @param from A \code{IDate} value. The start of the warning process.
#' @param end_dates A \code{IDate} vector. The endings of the warning process.
#' @param method A \code{character} signifying which smoothing method that is used. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing)
#' @param measure A \code{character} signifying which measure to smooth.
#' @param approximate A \code{logical} value. If set to \code{TRUE}, second differences are used to estimate gradients. If \code{FALSE} other methods are utilized.
#' @param tol A \code{double} value. What is the warning threshold for the smoothed curve?
#' @param dur A \code{integer} value. How many days must the gradient be above the threshold given in \code{tol} in order to trigger a warning
#' @param stringent_slope_warning A \code{logical} value. If set to \code{FALSE}, all triggered warnings are displayed. If set to \code{TRUE}, we avoid excessive number of warnings. See \code{?reduce_slope_warnings()} for more information.
#' @param bw A \code{double}. The bandwidth used for the smoothing.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#'
#' @return A \code{data.table} object containing the slope warning process information
#' @export
#'
#' @examples print(1)
slope_warning <- function(data, from, end_dates, method, measure, approximate = TRUE, tol = 2, dur = 3, stringent_slope_warning = TRUE, bw = 11, attach = TRUE){
  `Measured At` <- NULL
  smooth_list <- sapply(X = end_dates, FUN = function(x) kernel_smoothing0(data = data[`Measured At` >= from & `Measured At` <= x,], method = method, measure = measure, approximate = approximate, bw = 11, na_rm = TRUE), simplify = FALSE)
  last_smooth <- smooth_list[[length(smooth_list)]]
  if(approximate){
    smooth_list <- lapply(smooth_list, function(x) x[order(`Measured At`, decreasing = TRUE)][1:(min(dur + 1L, length(`Measured At`)))])
    warning_list <- lapply(smooth_list, FUN = function(x){
      if(nrow(x) >= 1){
        if(!all(is.na(x$`Measured At`))){
          number_of_weekend_days <- min(2, sum(weekdays(seq(from = min(x$`Measured At`, na.rm = TRUE), to = max(x$`Measured At`, na.rm = TRUE), by = "1 day")) %in% c("Saturday", "Sunday")))
        }
        else{
          number_of_weekend_days <- 0
        }
      }
      if(nrow(x) <= dur){
        measured_at <- max(x$`Measured At`, na.rm = TRUE)
        evaluated_days <- x$`Measured At`[-which.max(x$`Measured At`)]
        evaluated_slopes <- x$`Smoothed Gradient Degrees`[-which.max(x$`Measured At`)]
        raw_warnings <- data.table("Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                   "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                   "Measured At" = measured_at,
                                   "Mean Slope Magnitude" = ifelse(nrow(x) <= 2, 0, round(mean(evaluated_slopes, na.rm = TRUE), 3L)),
                                   "Condition Met" = FALSE)
      }
      else if(diff(range(x$`Measured At`)) >= dur + 1L + number_of_weekend_days){
        measured_at <- max(x$`Measured At`, na.rm = TRUE)
        evaluated_days <- x$`Measured At`[-which.max(x$`Measured At`)]
        evaluated_slopes <- x$`Smoothed Gradient Degrees`[-which.max(x$`Measured At`)]
        raw_warnings <- data.table("Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                   "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                   "Measured At" = measured_at,
                                   "Mean Slope Magnitude" = round(mean(evaluated_slopes, na.rm = TRUE), 3L),
                                   "Condition Met" = FALSE)
      }
      else{
        print(dur + 1L + number_of_weekend_days)
        measured_at <- max(x$`Measured At`, na.rm = TRUE)
        evaluated_days <- x$`Measured At`[-which.max(x$`Measured At`)]
        evaluated_slopes <- x$`Smoothed Gradient Degrees`[-which.max(x$`Measured At`)]
        all_above_threshold <- all(abs(evaluated_slopes) >= tol, na.rm = TRUE)
        all_same_direction <- all(sign(evaluated_slopes) == -1, na.rm = TRUE) | all(sign(evaluated_slopes) == 1, na.rm = TRUE)
        raw_warnings <- data.table("Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                   "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                   "Measured At" = measured_at,
                                   "Mean Slope Magnitude" = mean(evaluated_slopes, na.rm = TRUE),
                                   "Condition Met" = all_above_threshold & all_same_direction)
      }
      return(raw_warnings)
    })
  }
  else{
    smooth_list <- lapply(smooth_list, function(x) x[order(`Measured At`, decreasing = TRUE)][1:(min(dur, length(`Measured At`)))])
    warning_list <- lapply(smooth_list, FUN = function(x){
      if(nrow(x) >= 1){
        if(!all(is.na(x$`Measured At`))){
          number_of_weekend_days <- min(2, sum(weekdays(seq(from = min(x$`Measured At`, na.rm = TRUE), to = max(x$`Measured At`, na.rm = TRUE), by = "1 day")) %in% c("Saturday", "Sunday")))
        }
        else{
          number_of_weekend_days <- 0
        }
      }
      if(nrow(x) <= dur - 1L){
        measured_at <- max(x$`Measured At`, na.rm = TRUE)
        evaluated_days <- x$`Measured At`
        evaluated_slopes <- x$`Smoothed Gradient Degrees`
        raw_warnings <- data.table("Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                   "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                   "Measured At" = measured_at,
                                   "Mean Slope Magnitude" = ifelse(nrow(x) == 1, 0, round(mean(evaluated_slopes, na.rm = TRUE), 3L)),
                                   "Condition Met" = FALSE)
      }
      else if(diff(range(x$`Measured At`)) >= dur + number_of_weekend_days){
        measured_at <- max(x$`Measured At`, na.rm = TRUE)
        evaluated_days <- x$`Measured At`
        evaluated_slopes <- x$`Smoothed Gradient Degrees`
        raw_warnings <- data.table("Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                   "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                   "Measured At" = measured_at,
                                   "Mean Slope Magnitude" = round(mean(evaluated_slopes, na.rm = TRUE), 3L),
                                   "Condition Met" = FALSE)
      }
      else{
        measured_at <- max(x$`Measured At`, na.rm = TRUE)
        evaluated_days <- x$`Measured At`
        evaluated_slopes <- x$`Smoothed Gradient Degrees`
        all_above_threshold <- all(abs(evaluated_slopes) >= tol, na.rm = TRUE)
        all_same_direction <- all(sign(evaluated_slopes) == -1, na.rm = TRUE) | all(sign(evaluated_slopes) == 1, na.rm = TRUE)
        raw_warnings <- data.table("Evaluated Days" = paste(rev(evaluated_days), collapse = ", "),
                                   "Evaluated Slopes" = paste(rev(round(evaluated_slopes, 3L)), collapse = ", "),
                                   "Measured At" = measured_at,
                                   "Mean Slope Magnitude" = round(mean(evaluated_slopes, na.rm = TRUE), 3L),
                                   "Condition Met" = all_above_threshold & all_same_direction)
      }
      return(raw_warnings)
    })
  }
  extra_data <- data.table("Evaluated Days" = NA, "Evaluated Slopes" = NA, "Measured At" = as.IDate(from), "Mean Slope Magnitude" = 0, "Condition Met" = FALSE)
  warning_data <- rbindlist(warning_list)
  warning_data <- rbind(extra_data, warning_data)
  warning_data <- reduce_slope_warnings(warning_data, stringent = stringent_slope_warning) |> setDT()
  warning_data <- warning_data[duplicated(warning_data)==FALSE, ]
  if(attach){
    warning_data <- merge(last_smooth, warning_data, by = "Measured At")
  }
  return(warning_data)
}

#' Month-Year Bias Warning Process
#'
#' @param data A \code{data.table} containing measurements.
#' @param full_data A \code{data.table} containing all measurements for all dates.
#' @param from A \code{IDate} value. The start of the warning process.
#' @param end_dates A \code{IDate} vector. The endings of the warning process.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param tol A \code{double} value. What is the warning threshold for the month-year bias?
#' @param bw1 An \code{integer}. The window width used for the monthly moving median smoothing.
#' @param bw2 An \code{integer}. The window width used for the Yearly moving median smoothing.
#' @param snooze A \code{integer} signifying the number of days recurring warnings are snoozed.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#'
#' @return A \code{data.table} object containing the month-year bias warning process information
#' @export
#'
#' @examples print(1)
bias_warning <- function(data, full_data, from, end_dates, measure, tol, bw1 = 30, bw2 = 365, snooze = 10, attach = TRUE){
  `Measured At` <- NULL
  smooth_list_30 <- sapply(X = end_dates, FUN = function(x) moving_median_smoothing0(data = data[`Measured At` >= from & `Measured At` <= x,], measure = measure, bw = bw1, attach = "x"), simplify = FALSE)
  smooth_list_365 <- sapply(X = end_dates, FUN = function(x) moving_median_smoothing0(data = full_data[`Measured At` >= min(`Measured At`, na.rm = TRUE) & `Measured At` <= x,], measure = measure, bw = bw2, attach = "x"), simplify = FALSE)
  last_smooth <- smooth_list_30[[length(smooth_list_30)]]
  warning_list <- mapply(FUN = function(month, year){
    current_date <- max(month$`Measured At`, na.rm = TRUE)
    current_month_median <- month$m[which.max(month$`Measured At`)]
    current_year_median <- year$m[which.max(year$`Measured At`)]
    current_month_nor <- month$number_of_results[which.max(month$`Measured At`)]
    current_year_nor <- year$number_of_results[which.max(year$`Measured At`)]
    current_month_year_bias <- abs((current_month_median - current_year_median) / current_year_median) * 100
    current_enough_nor <- (current_month_nor >= 5) & (current_year_nor >= 5)
    return(data.table("Measured At" = current_date,
                      "Monthly Start Date" = max(current_date - 30 + 1, month$`Measured At`[which.min(month$`Measured At`)]),
                      "Monthly Median" = current_month_median,
                      "Monthly Number of Results" = current_month_nor,
                      "Yearly Median" = current_year_median,
                      "Yearly Start Date" = max(current_date - 365 + 1, year$`Measured At`[which.min(year$`Measured At`)]),
                      "Yearly Number of Results" = current_year_nor,
                      "Bias Percentage" = round(current_month_year_bias, 2L),
                      "Is Warning" = (current_month_year_bias > tol) & current_enough_nor))
  }, smooth_list_30, smooth_list_365, SIMPLIFY = FALSE)
  warning_data <- rbindlist(warning_list)
  warning_data$`Is Warning` <- reduce_bias_warnings(warnings = warning_data$`Is Warning`, snooze = snooze)
  if(attach){
    warning_data <- merge(last_smooth, warning_data, by = "Measured At")
  }
  return(warning_data)
}

#' Peer Group warning process
#'
#' @param data A \code{data.table} containing measurements.
#' @param full_data A \code{data.table} containing all measurements for all dates.
#' @param pg_data A \code{data.table} containing peer group measurements.
#' @param from A \code{IDate} value. The start of the warning process.
#' @param to A \code{IDate} value. The end date for the peer group.
#' @param end_dates A \code{IDate} vector. The endings of the warning process.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param tol A \code{double} value. What is the warning threshold for the peer group bias?
#' @param bw An \code{integer}. The window width used for the monthly moving median smoothing.
#' @param snooze snooze A \code{integer} signifying the number of days recurring warnings are snoozed.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#'
#' @return A \code{data.table} object containing the peer group bias warning process information
#' @export
#'
#' @examples print(1)
pg_warning <- function(data, full_data, pg_data, from, to, end_dates, measure, tol, bw = 30, snooze = 10, attach = TRUE){

  `Manufacturer Name` <- `Instrument Model Name` <- `Median` <- `Measured At` <- `Hypo Percentage` <- `Hyper Percentage` <- `Monthly Median` <- `Peer Group Monthly Median` <- `Monthly Number of Results` <- `Peer Group Monthly Number of Results` <- NULL

  # If 'pg_data' is not NULL, use it.
  if(!is.null(pg_data)){
    pg_data <- pg_data[`Manufacturer Name` %in% data$`Manufacturer Name` & `Instrument Model Name` %in% data$`Instrument Model Name`, ]

    # Taking date-wise medians
    if(measure == "median"){
      pg_data <- pg_data[, list(Median = ifelse(sum(!is.na(Median)) >= 1, median(Median, na.rm = TRUE), NA_real_)),
                         by = list(`Manufacturer Name`, `Instrument Model Name`, `Measured At`)]
    }
    else if(measure == "hypo"){
      pg_data <- pg_data[, list(`Hypo Percentage` = ifelse(sum(!is.na(`Hypo Percentage`)) >= 1, median(`Hypo Percentage`, na.rm = TRUE), NA_real_)),
                         by = list(`Manufacturer Name`, `Instrument Model Name`, `Measured At`)]
    }
    else if(measure == "hyper"){
      pg_data <- pg_data[, list(`Hyper Percentage` = ifelse(sum(!is.na(`Hyper Percentage`)) >= 1, median(`Hyper Percentage`, na.rm = TRUE), NA_real_)),
                         by = list(`Manufacturer Name`, `Instrument Model Name`, `Measured At`)]
    }
  }

  # If 'pg_data' is something else than NULL, use 'full_data'
  else{
    pg_data <- full_data
    # Taking date-wise medians
    if(measure == "median"){
      pg_data <- pg_data[, list(Median = ifelse(sum(!is.na(Median)) >= 1, median(Median, na.rm = TRUE), NA_real_)),
                         by = list(`Manufacturer Name`, `Instrument Model Name`, `Measured At`)]
    }
    else if(measure == "hypo"){
      pg_data <- pg_data[, list(`Hypo Percentage` = ifelse(sum(!is.na(`Hypo Percentage`)) >= 1, median(`Hypo Percentage`, na.rm = TRUE), NA_real_)),
                         by = list(`Manufacturer Name`, `Instrument Model Name`, `Measured At`)]
    }
    else if(measure == "hyper"){
      pg_data <- pg_data[, list(`Hyper Percentage` = ifelse(sum(!is.na(`Hyper Percentage`)) >= 1, median(`Hyper Percentage`, na.rm = TRUE), NA_real_)),
                         by = list(`Manufacturer Name`, `Instrument Model Name`, `Measured At`)]
    }
  }

  smooth_list_30 <- sapply(X = end_dates, FUN = function(x) moving_median_smoothing0(data = data[`Measured At` >= from & `Measured At` <= x,], measure = measure, bw = 30, attach = "x"), simplify = FALSE)
  last_smooth <- smooth_list_30[[length(smooth_list_30)]]

  # New end_dates bassed on 'pg_data'
  end_dates <- seq(from = min(pg_data$`Measured At`, na.rm = TRUE), to = as.IDate(to), by = "1 day")
  start_date <- min(pg_data$`Measured At`, na.rm = TRUE)
  if(length(end_dates) > 1){
    end_dates <- end_dates[which(end_dates %in% pg_data$`Measured At`)][-1]
    #end_dates <- end_dates[-1]
  }
  pg_smooth_list_30 <- sapply(X = end_dates, FUN = function(x) moving_median_smoothing0(data = pg_data[`Measured At` >= start_date & `Measured At` <= x,], measure = measure, bw = 30, attach = "x"), simplify = FALSE)
  pg_last_smooth <- pg_smooth_list_30[[length(pg_smooth_list_30)]]
  moving_median_instrument <- lapply(X = smooth_list_30, FUN = function(x){
    current_date <- max(x$`Measured At`, na.rm = TRUE)
    current_month_median <- x$m[which.max(x$`Measured At`)]
    current_month_nor <- x$number_of_results[which.max(x$`Measured At`)]
    return(data.table("Measured At" = current_date,
                      "Monthly Median" = current_month_median,
                      "Monthly Start Date" = max(current_date - 30 + 1, x$`Measured At`[which.min(x$`Measured At`)]),
                      "Monthly Number of Results" = current_month_nor))
  })
  moving_median_peer_group <- lapply(X = pg_smooth_list_30, FUN = function(x){
    current_date <- max(x$`Measured At`, na.rm = TRUE)
    current_month_median <- x$m[which.max(x$`Measured At`)]
    current_month_nor <- x$number_of_results[which.max(x$`Measured At`)]
    return(data.table("Measured At" = current_date,
                      "Peer Group Monthly Median" = x$m[which.max(x$`Measured At`)],
                      "Peer Group Start Date" = max(current_date - 30 + 1, x$`Measured At`[which.min(x$`Measured At`)]),
                      "Peer Group Monthly Number of Results" = current_month_nor))
  })
  moving_median_instrument <- rbindlist(moving_median_instrument)
  moving_median_peer_group <- rbindlist(moving_median_peer_group)
  warning_data <- merge(moving_median_instrument, moving_median_peer_group, by = "Measured At", allow.cartesian = TRUE)

  warning_data[, `:=` ("Bias Percentage" = round(abs((`Monthly Median` - `Peer Group Monthly Median`) / `Peer Group Monthly Median`) * 100, 2L),
                       "Is Warning" = abs((`Monthly Median` - `Peer Group Monthly Median`) / `Peer Group Monthly Median`) * 100 > tol & (`Monthly Number of Results` > 5) & (`Peer Group Monthly Number of Results` > 5))]
  warning_data$`Is Warning` <- reduce_bias_warnings(warnings = warning_data$`Is Warning`, snooze = snooze)
  if(attach){
    warning_data <- merge(last_smooth, warning_data, by = "Measured At")
  }
  return(warning_data)
}

#' Simulate A Warning Process for One Instrument Within One Laboratory
#'
#' @param data A \code{data.table} containing measurements.
#' @param pg_data A \code{data.table} containing peer group measurements.
#' @param from A \code{IDate} value. The start of the warning process.
#' @param to A \code{IDate} value. The end of the warning process.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param method A \code{character} signifying which smoothing method that is used. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing)
#' @param warning A \code{character} signifying which warning process that is of interest to simulate.
#' @param tol A \code{double} value. What is the warning threshold for relevant \code{warning}?
#' @param dur dur A \code{integer} value. How many days must the gradient be above the threshold given in \code{tol} in order to trigger a warning. Only relevant if \code{warning} is set to \code{slope}.
#' @param bw A \code{double}. The bandwidth used for the kernel smoothing.
#' @param snooze A \code{integer} signifying the number of days recurring warnings are snoozed. Only relevant if \code{warning} is \code{bias} or \code{peer_group}.
#' @param stringent_slope_warning A \code{logical} value. If set to \code{FALSE}, all triggered warnings are displayed. If set to \code{TRUE}, we avoid excessive number of warnings. See \code{?reduce_slope_warnings()} for more information. Only relevant if \code{warning} is set to \code{slope}.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#' @param approximate A \code{logical} value. If set to \code{TRUE}, second differences are used to estimate gradients. If \code{FALSE} other methods are utilized. Only relevant if \code{warning} is set to \code{slope}.
#'
#' @return A \code{data.table} object containing the warning process information
#' @export
#'
#' @examples print(1)
warning_process0 <- function(data, pg_data = NULL, from = "2021-01-01", to = "2022-12-31", measure = c("median", "hyper", "hypo"), method = c("lc", "ll"), warning = c("slope", "bias", "peer_group"), tol = 4, dur = 3, bw = 11, snooze = 10L, stringent_slope_warning = TRUE, attach = TRUE, approximate = TRUE){

  `average` <- `Median` <- `standard_deviation` <- `Measured At` <- NULL
  #cat("Laboratory + Instrument Code with ID ", round(runif(1, min = 1, 1e3L),0L), sample(c("A", "B", "C", "D", "E", "F"), size = 1), "\n", sep = "")
  measure <- measure[1]
  warning <- warning[1]
  method <- method[1]

  # If average and standard deviation does not exist, calculate them if slope warning
  if(warning == "slope" && (!all(c("average", "standard_deviation") %in% names(data)))){
    if(measure == "median"){
      data[, average := mean(Median, na.rm = TRUE)]
      data[, standard_deviation := sd(Median, na.rm = TRUE)]
    }
    else{
      data[, average := 0]
      data[, standard_deviation := 1]
    }
  }

  if(as.IDate(from) < min(data$`Measured At`, na.rm = TRUE)){
    from <- min(data$`Measured At`, na.rm = TRUE)
  }
  else if(as.IDate(from) > max(data$`Measured At`, na.rm = TRUE)){
    from <- min(data$`Measured At`, na.rm = TRUE)
  }
  if(as.IDate(to) > max(data$`Measured At`, na.rm = TRUE)){
    to <- max(data$`Measured At`, na.rm = TRUE)
  }
  else if(as.IDate(to) < min(data$`Measured At`, na.rm = TRUE)){
    to <- max(data$`Measured At`, na.rm = TRUE)
  }
  full_data <- copy(data)[`Measured At` <= to, ]
  data <- data[`Measured At` >= from & `Measured At` <= to, ]
  end_dates <- seq(from = as.IDate(from), to = as.IDate(to), by = "1 day")
  if(length(end_dates) > 1){
    end_dates <- end_dates[which(end_dates %in% data$`Measured At`)][-1]
  }

  if(warning == "slope"){
    warning_data <- slope_warning(data, from, end_dates, method, measure, approximate, tol, dur, stringent_slope_warning, bw, attach)
    return(warning_data)
  }
  else if(warning == "bias"){
    warning_data <- bias_warning(data, full_data, from, end_dates, measure, tol, 30, 365, snooze, attach)
    return(warning_data)
  }
  else if(warning == "peer_group"){
    warning_data <- pg_warning(data, full_data, pg_data, from, to, end_dates, measure, tol, 30, snooze, attach)
    return(warning_data)
  }
}

#' Simulate A Warning Process for Multiple Instruments Within Multiple Laboratories
#'
#' @param data A \code{data.table} containing measurements.
#' @param by A \code{character} vector containing the names of the grouping variables.
#' @param pg_data A \code{data.table} containing peer group measurements.
#' @param from A \code{IDate} value. The start of the warning process.
#' @param to A \code{IDate} value. The end of the warning process.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param method A \code{character} signifying which smoothing method that is used. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing)
#' @param warning A \code{character} signifying which warning process that is of interest to simulate.
#' @param tol A \code{double} value. What is the warning threshold for relevant \code{warning}?
#' @param dur A \code{integer} value. How many days must the gradient be above the threshold given in \code{tol} in order to trigger a warning. Only relevant if \code{warning} is set to \code{slope}.
#' @param bw A \code{double}. The bandwidth used for the kernel smoothing.
#' @param snooze A \code{integer} signifying the number of days recurring warnings are snoozed. Only relevant if \code{warning} is \code{bias} or \code{peer_group}.
#' @param stringent_slope_warning A \code{logical} value. If set to \code{FALSE}, all triggered warnings are displayed. If set to \code{TRUE}, we avoid excessive number of warnings. See \code{?reduce_slope_warnings()} for more information. Only relevant if \code{warning} is set to \code{slope}.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#' @param approximate A \code{logical} value. If set to \code{TRUE}, second differences are used to estimate gradients. If \code{FALSE} other methods are utilized. Only relevant if \code{warning} is set to \code{slope}.
#'
#' @return A \code{data.table} object containing the warning process information
#' @export
#'
#' @examples print(1)
warning_process <- function(data, by, pg_data = NULL, from = "2021-01-01", to = "2022-12-31", measure = c("median", "hyper", "hypo"), method = c("lc", "ll"), warning = c("slope", "bias", "peer_group"), tol = 4, dur = 3, bw = 11, snooze = 10L, stringent_slope_warning = TRUE, attach = TRUE, approximate = TRUE){

  `Measured At` <- `Median` <- NULL
  measure <- measure[1]
  ksc <- kernel_smoothing_check(data = data, by = by, measure = measure, only_acceptable = TRUE)
  if(!isFALSE(ksc)){
    indicator <- !logical(nrow(data))
    for(i in 1:length(by)){
      indicator <- indicator & (data[[which(names(data) == by[i])]] %in% ksc[[i]])
    }
    data <- data[which(indicator),]
  }
  else{
    stop("Not enough results to perform smoothing.")
  }

  method <- method[1]
  warning <- warning[1]
  if(measure == "median"){
    average_and_standard_deviation <- data[`Measured At` <= as.IDate(from),
                                           list(average = mean(`Median`, na.rm = TRUE),
                                                #standard_deviation = sqrt(var(diff(`Median`), na.rm = TRUE)/2)
                                                standard_deviation = sd(`Median`, na.rm = TRUE)),
                                           by = by]

    data <- merge(data, average_and_standard_deviation, by = by)
  }
  else{
    average_and_standard_deviation <- data[, list(average = 0,
                                                  standard_deviation = 1),
                                           by = by]

    data <- merge(data, average_and_standard_deviation, by = by)
  }
  if(warning == "peer_group" && is.null(pg_data)){
    pg_data <- copy(data)
  }

  output <- data[, warning_process0(data = .SD, from = from, to = to, measure = measure, method = method, warning = warning, tol = tol, attach = attach, dur = dur, bw = bw, approximate = approximate, pg_data = pg_data, snooze = snooze, stringent_slope_warning = stringent_slope_warning), by = by]
  return(output)
}
