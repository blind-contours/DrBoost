% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DrBoost.R
\name{build_box_string}
\alias{build_box_string}
\title{Build Human-Readable Box String}
\usage{
build_box_string(lb, ub, feat_names, digits = 3)
}
\arguments{
\item{lb}{numeric vector of lower bounds}

\item{ub}{numeric vector of upper bounds}

\item{feat_names}{character vector of length p (feature names)}

\item{digits}{how many decimal places to round}
}
\value{
A single character string describing the bounding box constraints
}
\description{
Creates a short text rule like
"Petal.Length >= 2.3 & Petal.Length <= 5 & Petal.Width <= 1.0"
from numeric lb/ub and feature names.
}
\keyword{internal}
