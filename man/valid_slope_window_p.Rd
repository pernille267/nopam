% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/valid_slope_window.R
\name{valid_slope_window_p}
\alias{valid_slope_window_p}
\title{Checks Whether a Smoothing Window is Valid - P}
\usage{
valid_slope_window_p(window, measure, dur = 3L)
}
\arguments{
\item{window}{A \code{data.table} object containing the relevant smoothing data for a particular window.}

\item{measure}{A \code{character} value referring to the measure we want to smooth. Available choices include \code{median}, \code{hyper} and \code{hypo}. Nothihng else is accepted.}

\item{dur}{A \code{integer} value. Signifying the number of days the slope must exceed the threshold value to return a warning.}
}
\value{
A \code{logical} value. If the requirements are satisfied, returns \code{TRUE}.
}
\description{
Checks Whether a Smoothing Window is Valid - P
}
\examples{
print(1)
}
