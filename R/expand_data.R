#' Internal Function for Data Conversion
#'
#' This function is used internally for converting different data types to data.table.
#' @keywords internal
convert_to_data_table <- function(data){
  if(is.list(data)){
    setDT(data)
  }
  else if(is.data.frame(data)){
    as.data.table(data)
  }
  else{
    stop("'data' is expected to be a data.table object but is a ", class(data))
  }
}


#' Expand Data to Include Entries Inbetween Missing Dates
#'
#' @param data The dataset to be expanded. Must be a \code{data.table} object. Are attempted converted if it is not a \code{data.table} object.
#' @param id_cols A \code{character} vector of column names to group by. All column names included must be part of \code{data}.
#' @param date_col A \code{character} value corresponding to the name of the date column. Must be part of \code{data}.
#' @param measurement_cols A \code{character} vector corresponding to the names of the measurement columns. All column names included must be part of \code{data}. Currently not used.
#'
#' @return A \code{data.table} object that is the expanded version of the input \code{data}
#' @export
#'
#' @examples print(1)
expand_data <- function(data, id_cols = NULL, date_col = NULL, measurement_cols = NULL){

  # Checks if 'data' is of correct class. Convert if possible, error thrown if not.
  data <- convert_to_data_table(data)

  # Handle NULLs
  if(is.null(id_cols)){
    id_cols <- c("Analyte Name", "Country", "Manufacturer Name", "Instrument Model Name", "Laboratory Code", "Instrument Id", "Instrument Code")
  }
  if(is.null(date_col)){
    date_col <- "Measured At"
  }
  if(is.null(measurement_cols)){
    measurement_cols <- c("Median", "Hyper Percentage", "Hypo Percentage")
  }

  # Get column names of 'data'
  all_cols <- names(data)

  # Matching
  id_cols_ids <- id_cols %in% all_cols
  date_col_id <- date_col %in% all_cols
  measurement_cols_ids <- measurement_cols %in% all_cols

  # Checks
  id_cols_exist <- all(id_cols_ids)
  date_col_exists <- any(date_col_id)
  measurement_cols_exist <- all(measurement_cols_ids)

  # Throw an error if some of the given columns do not exist
  if(isFALSE(id_cols_exist)){
    stop("The following ID columns do not exist in 'data': ", id_cols[!id_cols_ids], ".")
  }
  if(isFALSE(date_col_exists)){
    stop("The following date column does not exist in 'data': ", date_col[!date_col_id], ".")
  }
  if(isFALSE(measurement_cols_exist)){
    stop("The following measurement columns do not exist in 'data': ", measurement_cols[!measurement_cols_ids], ".")
  }
  if(sum(names(data) %in% date_col) > 1){
    warning("Several matching colums are found for ", date_col, ".", " Only the first one is used.")
  }

  # Rename date_col name to start name, a.k.a, 'Measured At'
  names(data)[which(names(data) %in% date_col)[1]] <- "Measured At"
  data_list <- split(x = data, by = id_cols, keep.by = TRUE, sorted = FALSE)

  # Inner function for expansion
  expand_group_data <- function(data, id_cols){
    all_dates <- seq(min(data[["Measured At"]], na.rm = TRUE), max(data[["Measured At"]], na.rm = TRUE), by = "1 day")
    n_rows <- length(all_dates)
    expanded_group_data <- CJ("Measured At" = all_dates)
    for(id in id_cols){
      expanded_group_data[[id]] <- rep(data[[id]][1], n_rows)
    }

    out <- merge(expanded_group_data, data, by = c("Measured At", id_cols), all.x = TRUE)
    out$`Is Weekend` <- ifelse(weekdays(out$`Measured At`) == "Saturday" | weekdays(out$`Measured At`) == "Sunday", TRUE, FALSE)
    out$`Weekday` <- weekdays(out$`Measured At`)
    return(out)
  }

  result <- lapply(data_list, expand_group_data, id_cols = id_cols) |> rbindlist(idcol = NULL)
  return(result)

}
