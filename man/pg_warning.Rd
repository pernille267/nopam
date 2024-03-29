% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_warning_process_p.R
\name{pg_warning}
\alias{pg_warning}
\title{Peer Group warning process}
\usage{
pg_warning(
  data,
  full_data,
  pg_data,
  from,
  to,
  end_dates,
  measure,
  tol,
  bw = 30,
  snooze = 10,
  attach = TRUE
)
}
\arguments{
\item{data}{A \code{data.table} containing measurements.}

\item{full_data}{A \code{data.table} containing all measurements for all dates.}

\item{pg_data}{A \code{data.table} containing peer group measurements.}

\item{from}{A \code{IDate} value. The start of the warning process.}

\item{to}{A \code{IDate} value. The end date for the peer group.}

\item{end_dates}{A \code{IDate} vector. The endings of the warning process.}

\item{measure}{A \code{character} signifying which measure to smooth.}

\item{tol}{A \code{double} value. What is the warning threshold for the peer group bias?}

\item{bw}{An \code{integer}. The window width used for the monthly moving median smoothing.}

\item{snooze}{snooze A \code{integer} signifying the number of days recurring warnings are snoozed.}

\item{attach}{A \code{logical} value. Should we attach the information found in \code{data} to the output?}
}
\value{
A \code{data.table} object containing the peer group bias warning process information
}
\description{
Peer Group warning process
}
\examples{
print(1)
}
