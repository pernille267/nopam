#' Check if we have enough observations to do kernel smoothing correctly
#'
#' @param data A \code{data.table} containing measurements.
#' @param by A \code{character} vector containing the names of the grouping variables.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param tol A \code{integer}. The minimum acceptable number of results. Only groups with more than this value will be acceptable if \code{only_acceptable} is set to \code{TRUE}.
#' @param only_acceptable A \code{logical} value. If set to \code{TRUE}, only the grouping variables with enough observations are presented. Otherwise \code{FALSE} is returned by the function. If set to \code{FALSE}, a \code{data.table} object containing relevant information is outputted.
#'
#' @return A \code{logical} value or a \code{data.table} object depending on \code{data}.
#' @export
#'
#' @examples print(1)

kernel_smoothing_check <- function(data, by, measure = c("median", "hyper", "hypo"), tol = 2, only_acceptable = TRUE){

  `m` <- `Median` <- `Hyper Percentage` <- `Hypo Percentage` <- `Measured At` <- `Laboratory Code` <- `Number of results` <- NULL

  measure <- measure[1]
  if(measure == "median"){
    data[, m := Median,]
  }
  else if(measure == "hyper"){
    data[, m := `Hyper Percentage`]
  }
  else if(measure == "hypo"){
    data[, m := `Hypo Percentage`]
  }
  output <- data[, list(`Number of results` = sum(!is.na(m)),
                        `Number of missing results` = sum(is.na(m)),
                        `Number of days in period` = diff(range(`Measured At`, na.rm = TRUE)) + 1,
                        `Start Date` = min(`Measured At`, na.rm = TRUE),
                        `End Date` = max(`Measured At`, na.rm = TRUE),
                        `Date Dispersion` = diff(range(`Measured At`, na.rm = TRUE)) / max(sum(!is.na(m)), 1)),
                 by = by]
  if("Laboratory Code" %in% names(output)){
    setorder(output, `Laboratory Code`)
  }
  if(only_acceptable){
    output <- output[`Number of results` >= tol,]
    if(nrow(output) < 1){
      return(FALSE)
    }
    return(output)
  }
  else{
    return(output)
  }

}
