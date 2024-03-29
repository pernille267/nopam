% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_warning_process_p.R
\name{bias_warning}
\alias{bias_warning}
\title{Month-Year Bias Warning Process}
\usage{
bias_warning(
  data,
  full_data,
  from,
  end_dates,
  measure,
  tol,
  bw1 = 30,
  bw2 = 365,
  snooze = 10,
  attach = TRUE
)
}
\arguments{
\item{data}{A \code{data.table} containing measurements.}

\item{full_data}{A \code{data.table} containing all measurements for all dates.}

\item{from}{A \code{IDate} value. The start of the warning process.}

\item{end_dates}{A \code{IDate} vector. The endings of the warning process.}

\item{measure}{A \code{character} signifying which measure to smooth.}

\item{tol}{A \code{double} value. What is the warning threshold for the month-year bias?}

\item{bw1}{An \code{integer}. The window width used for the monthly moving median smoothing.}

\item{bw2}{An \code{integer}. The window width used for the Yearly moving median smoothing.}

\item{snooze}{A \code{integer} signifying the number of days recurring warnings are snoozed.}

\item{attach}{A \code{logical} value. Should we attach the information found in \code{data} to the output?}
}
\value{
A \code{data.table} object containing the month-year bias warning process information
}
\description{
Month-Year Bias Warning Process
}
\examples{
print(1)
}
