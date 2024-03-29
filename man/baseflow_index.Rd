% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signatures.R
\name{baseflow_index}
\alias{baseflow_index}
\title{Compute the baseflow index}
\usage{
baseflow_index(x, method = "LH", min_segment_length = 30, ...)
}
\arguments{
\item{x}{A tibble.}

\item{method}{A string specifying the method to use.}

\item{min_segment_length}{An integer specifying the
minimum number of consecutive time points for which to
perform baseflow separation.}

\item{...}{Additional arguments passed to \code{baseflow_sep}}
}
\value{
A tibble
}
\description{
Compute the baseflow index
}
\examples{
\dontrun{
print("Hello, world")
}
}
