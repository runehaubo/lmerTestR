% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contrast_utils.R
\name{ensure_full_rank}
\alias{ensure_full_rank}
\title{Ensure a Design Matrix has Full (Column) Rank}
\usage{
ensure_full_rank(X, tol = 1e-07, silent = FALSE, test.ans = FALSE)
}
\arguments{
\item{X}{a design matrix as produced by \code{model.matrix}.}

\item{tol}{\code{qr} tolerance.}

\item{silent}{throw message if columns are dropped from \code{X}? Default
is \code{FALSE}.}

\item{test.ans}{Test if the resulting/returned matrix has full rank? Default
is \code{FALSE}.}
}
\value{
A design matrix in which redundant columns are dropped
}
\description{
Determine and drop redundant columns using the \code{\link{qr}}
decomposition.
}
\keyword{internal}
