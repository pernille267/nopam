#' Moving Median Curve for One Instrument Within One Laboratory
#'
#' @param data A \code{data.table} containing measurements.
#' @param measure  A \code{character} signifying which measure to smooth.
#' @param bw An \code{integer}. The window width used for the moving median smoothing.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#'
#' @return A \code{data.table} containing the smoothed values for one instrument within one laboratory
#' @export
#'
#' @examples print(1)
moving_median_smoothing0 <- function(data, measure = c("median", "hyper", "hypo"), bw = 30, attach = TRUE){
  `m` <- `Median` <- `Hyper Percentage` <- `Hypo Percentage` <- NULL
  measure <- measure[1]
  data$date <- as.numeric(data$`Measured At` - min(data$`Measured At`, na.rm = TRUE)) + 1
  date_NAs <- which(is.na(data$date))
  smooth_name <- ""
  raw_name <- ""
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

  moving_median_smoothed <- moving_median(date = sub_data$date, median = sub_data$m, bandwidth = bw)
  output <- merge(sub_data, moving_median_smoothed, by = "date")
  setnames(output, old = c("moving_median"), new = c(smooth_name))

  if(isTRUE(attach)){
    output$m <- NULL
    output$date <- NULL
    return(output)
  }
  else if(attach == "x"){
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
  else{
    output <- output[, c("Measured At", raw_name, smooth_name), with = FALSE]
    return(output)
  }
  return(output)

}

#' Moving Median Curve for Multiple Instruments Within Multiple Laboratories
#'
#' @param data A \code{data.table} containing measurements.
#' @param by A \code{character} vector containing the names of the grouping variables.
#' @param measure  A \code{character} signifying which measure to smooth.
#' @param bw An \code{integer}. The window width used for the moving median smoothing.
#' @param attach A \code{logical} value. Should we attach the information found in \code{data} to the output?
#'
#' @return A \code{data.table} containing the smoothed values for multiple instruments within multiple laboratories
#' @export
#'
#' @examples print(1)
moving_median_smoothing <- function(data, by, measure = c("median", "hyper", "hypo"), bw = 30, attach = TRUE){
  measure <- measure[1]
  output <- data[, moving_median_smoothing0(data = .SD, measure = measure, bw = bw, attach = attach), by = by]
  return(output)
}
