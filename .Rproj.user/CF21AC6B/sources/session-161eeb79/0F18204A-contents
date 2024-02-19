#' Kernel Smoothing for one Instrument Within one Laboratory
#'
#' @param data A \code{data.table} containing measurements
#' @param method A \code{character} signifying which smoothing method that is used. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing)
#' @param measure A \code{character} signifying which measure to smooth.
#' @param bw A \code{double}. The bandwidth used for the smoothing.
#' @param average A \code{double}. For calculation of the smoothed median curve slopes. Should be equal to 0 for hypo and hyper.
#' @param standard_deviation A \code{double}. For calculation of the smoothed median curve slopes. Should be equal to 1 for hypo and hyper.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#' @param approximate A \code{logical} value. If set to \code{TRUE}, second differences are used to estimate gradients. If \code{FALSE} other methods are utilized.
#' @param matrix_approach A \code{logical} value. Should matrix calculus be used to smooth. Only relevant if \code{method} is \code{ll}. Matrix calculus is faster, but may not be stable in some cases.
#' @param diagnostics A \code{logical} value. Should diagnostic messages be printed where relevant?
#' @param na_rm A \code{logical} value. Should rows with missing measurements be stripped from the output?
#'
#' @return A \code{data.table} containing the smoothed values for one instrument within one laboratory
#' @export
#'
#' @examples print(1)

kernel_smoothing0 <- function(data, method = "lc", measure = c("median", "hyper", "hypo"), bw = 11, average = 0, standard_deviation = 1, attach = TRUE, approximate = TRUE, matrix_approach = FALSE, diagnostics = FALSE, na_rm = TRUE){
  `m` <- `Median` <- `Hyper Percentage` <- `Hypo Percentage` <- NULL
  # Get measure (can be median, hyper or hypo)
  measure <- measure[1]
  # Convert from "Date" class to "numeric" class.
  data$date <- as.numeric(data$`Measured At` - min(data$`Measured At`, na.rm = TRUE)) + 1
  # Detect NA-values in the dates, which are removed later
  date_NAs <- which(is.na(data$date))
  # Initialize output names
  smooth_name <- ""
  raw_name <- ""
  # Convert our measure to "m", so that is common for all measures
  if(measure == "median"){
    data[, m := Median]
    smooth_name <- "Smoothed Median"
    raw_name <- "Median"
  }
  else if(measure == "hyper"){
    data[, m := `Hyper Percentage`]
    smooth_name <- "Smoothed Hyper Percentage"
    raw_name <- "Hyper Percentage"
  }
  else if(measure == "hypo"){
    data[, m := `Hypo Percentage`]
    smooth_name <- "Smoothed Hypo Percentage"
    raw_name <- "Hypo Percentage"
  }
  else{
    stop("measure must be either median, hyper or hypo.")
  }

  if(length(date_NAs) > 0){
    sub_data <- data[-date_NAs,]
  }
  else{
    sub_data <- data
  }

  if(method == "lc"){
    if(is.null(bw)){
      bw <- 11
    }
    # Check if 'average' and 'standard_deviation' is in 'sub_data' (should be the case for slope warnings!)
    if(all(c("average", "standard_deviation") %in% names(sub_data))){
      kernel_smoothed <- laks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = sub_data$average[1], standard_deviation = sub_data$standard_deviation[1], approximate = approximate)
    }
    # Try to handle the case when at leat one of 'average' and 'standard_deviation' cannot be found in 'sub_data'
    else{
      # If both 'average' >= 0 and 'standard_deviation' > 0 (spans the default settings), we use these values!
      if(isTRUE(average >= 0 && standard_deviation > 0)){
        kernel_smoothed <- laks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = average, standard_deviation = standard_deviation, approximate = approximate)
      }
      # If at least one of 'average' and 'standard_deviation' < 0, 'average' and 'standard_deviation' should be calculated in the laks() function! Negative values for these in the laks() function trigger their calculations.
      else if(isTRUE(average < 0 || standard_deviation < 0)){
        # If we consider hyper / hypo values, we do not wish to calculate 'average' and 'standard_deviation' inside laks(), as they will always be 0 and 1 in these cases.
        if(measure == "hyper" || measure == "hypo"){
          kernel_smoothed <- laks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = 0, standard_deviation = 1, approximate = approximate)
        }
        # Last resort!
        else{
          if(diagnostics){
            cat("Average and standard deviation not found in data, and are calculated via the laks() function!", "\n")
          }
          kernel_smoothed <- laks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = -1, standard_deviation = -1, approximate = approximate)
        }

      }
      else{
        if(diagnostics){
          cat("Average and standard deviation not found in data. If this is not intentional, you should have a closer look", "\n")
        }
        kernel_smoothed <- laks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = 0, standard_deviation = 1, approximate = approximate)
      }

    }
    number_of_valid_results <- length(sub_data$m[!is.na(sub_data$m)])
    if(diagnostics && number_of_valid_results <= 10){
      cat("There are very few valid observations in this dataset (", number_of_valid_results, "). The output will not be particularly informative...", "\n", sep = "")
    }
    mean_absolute_residual_error <- mean(abs((sub_data$m[!is.na(sub_data$m)] - kernel_smoothed$smoothed_median)/kernel_smoothed$smoothed_median) * 100, na.rm = TRUE)
    setDT(kernel_smoothed)
  }
  else if(method == "ll"){
    if(is.null(bw)){
      bw <- 11
    }
    if(nrow(sub_data) >= 1){
      if(all(c("average", "standard_deviation") %in% names(sub_data))){
        kernel_smoothed <- llks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = sub_data$average[1], standard_deviation = sub_data$standard_deviation[1], approximate = approximate, matrix_approach = matrix_approach)
      }
      # Try to handle the case when at leat one of 'average' and 'standard_deviation' cannot be found in 'sub_data'
      else{
        if(isTRUE(average >= 0 && standard_deviation > 0)){
          kernel_smoothed <- llks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = average, standard_deviation = standard_deviation, approximate = approximate, matrix_approach = matrix_approach)
        }
        else if(isTRUE(average < 0 || standard_deviation < 0)){
          if(measure == "hyper" || measure == "hypo"){
            kernel_smoothed <- llks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = 0, standard_deviation = 1, approximate = approximate, matrix_approach = matrix_approach)
          }
          else{
            if(diagnostics){
              cat("Average and standard deviation not found in data, and are calculated via the llks() function!", "\n")
            }
            kernel_smoothed <- llks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = -1, standard_deviation = -1, approximate = approximate, matrix_approach = matrix_approach)
          }
        }
        else{
          if(diagnostics){
            cat("Average and standard deviation not found in data. If this is not intentional, you should have a closer look", "\n")
          }
          kernel_smoothed <- llks(date = sub_data$date, median = sub_data$m, bandwidth = bw, average = 0, standard_deviation = 1, approximate = approximate, matrix_approach = matrix_approach)
        }

      }
      number_of_valid_results <- length(sub_data$m[!is.na(sub_data$m)])
      if(diagnostics && number_of_valid_results <= 10){
        cat("There are very few valid observations in this dataset (", number_of_valid_results, "). The output will not be particularly informative...", "\n", sep = "")
      }
      mean_absolute_residual_error <- mean(abs((sub_data$m[!is.na(sub_data$m)] - kernel_smoothed$smoothed_median)/kernel_smoothed$smoothed_median) * 100, na.rm = TRUE)
    }
    else{
      if(diagnostics){
        cat("There are very few valid observations in this dataset (", 0, "). Cannot generate meaningful output...", "\n", sep = "")
      }
      kernel_smoothed <- list("date" = data$date, "smoothed_median" = rep(NA_real_, length(data$date)), "gradients" = rep(NA_real_, length(data$date)))
      number_of_valid_results <- 0
      mean_absolute_residual_error <- NA_real_
    }
    setDT(kernel_smoothed)
  }
  if(nrow(sub_data) >= 1){
    output <- merge(sub_data, kernel_smoothed, by = "date", all.x = !na_rm)
  }
  else{
    output <- merge(data, kernel_smoothed, by = "date", all.x = !na_rm)
  }
  setnames(output, old = c("smoothed_median", "gradients"), new = c(smooth_name, "Smoothed Gradient Degrees"))
  if(isTRUE(attach)){
    output$bw <- bw
    output$mean_absolute_residual_error <- mean_absolute_residual_error
    output$date <- NULL
    output$m <- NULL
    return(output)
  }
  else if(attach == "x"){
    output$date <- NULL
    if(measure == "median"){
      output$m <- output$`Smoothed Median`
    }
    else if(measure == "hyper"){
      output$m <- output$`Smoothed Hyper Percentage`
    }
    else if(measure == "hypo"){
      output$m <- output$`Smoothed Hypo Percentage`
    }
  }
  else if(attach == "y"){
    output$date <- NULL
    if(measure == "median"){
      output$sm <- output$`Smoothed Median`
      output$m <- output$`Median`
    }
    else if(measure == "hyper"){
      output$sm <- output$`Smoothed Hyper Percentage`
      output$m <- output$`Hyper Percentage`
    }
    else if(measure == "hypo"){
      output$sm <- output$`Smoothed Hypo Percentage`
      output$m <- output$`Hypo Percentage`
    }
    output <- output[, c("Measured At", raw_name, smooth_name, "Smoothed Gradient Degrees", "m", "sm"), with = FALSE]
    return(output)
  }
  else{
    output <- output[, c("Measured At", raw_name, smooth_name, "Smoothed Gradient Degrees"), with = FALSE]
    return(output)
  }
  return(output)
}

#' Kernel Smoothing for Multiple Instruments Within Multiple Laboratories
#'
#' @param data A \code{data.table} containing measurements.
#' @param by A \code{character} vector containing the names of the grouping variables.
#' @param method A \code{character} signifying which smoothing method that is used. Valid options include \code{lc} (Local-Average kernel smoothing) and \code{ll} (Local-Linear kernel smoothing)
#' @param measure A \code{character} signifying which measure to smooth.
#' @param bw A \code{double}. The bandwidth used for the smoothing.
#' @param average A \code{double}. For calculation of the smoothed median curve slopes. Should be equal to 0 for hypo and hyper.
#' @param standard_deviation A \code{double}. For calculation of the smoothed median curve slopes. Should be equal to 1 for hypo and hyper.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#' @param approximate A \code{logical} value. If set to \code{TRUE}, second differences are used to estimate gradients. If \code{FALSE} other methods are utilized.
#' @param matrix_approach A \code{logical} value. Should matrix calculus be used to smooth. Only relevant if \code{method} is \code{ll}. Matrix calculus is faster, but may not be stable in some cases.
#'
#' @return A \code{data.table} containing the smoothed values for multiple instruments within multiple laboratories
#' @export
#'
#' @examples print(1)
kernel_smoothing <- function(data, by, method = "lc", measure = c("median", "hyper", "hypo"), bw = 11, average = 0, standard_deviation = 1, attach = TRUE, approximate = TRUE, matrix_approach = TRUE){
  `Median` <- NULL
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
  if(measure == "median"){
    average_and_standard_deviation <- data[, list(average = mean(`Median`, na.rm = TRUE),
                                                  standard_deviation = ifelse(length(`Median`)<=2, sd(`Median`, na.rm = TRUE), sqrt(var(diff(`Median`), na.rm = TRUE)/2))),
                                           by = by]

    data <- merge(data, average_and_standard_deviation, by = by)
  }
  else{
    average_and_standard_deviation <- data[, list(average = 0,
                                                  standard_deviation = 1),
                                           by = by]

    data <- merge(data, average_and_standard_deviation, by = by)
  }

  output <- data[, kernel_smoothing0(.SD, method = method, measure = measure, bw = bw, attach = attach, average = average, standard_deviation = standard_deviation, approximate = approximate, matrix_approach = matrix_approach), by = by]
  return(output)
}
