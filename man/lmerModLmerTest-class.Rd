% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmer.R
\docType{class}
\name{lmerModLmerTest-class}
\alias{lmerModLmerTest-class}
\alias{lmerModLmerTest}
\title{Represent Linear Mixed-Effects Models}
\value{
An object of class \code{lmerModLmerTest} with slots as in
\code{lmerMod} objects (see \code{\link[lme4]{merMod}}) and a few
additional slots as described in the slots section.
}
\description{
The \code{lmerModLmerTest} class extends \code{lmerMod} (which extends
\code{merMod}) from the \pkg{lme4}-package.
}
\section{Slots}{

\describe{
\item{\code{vcov_varpar}}{a numeric matrix holding the asymptotic variance-covariance
matrix of the variance parameters (including sigma).}

\item{\code{Jac_list}}{a list of gradient matrices (Jacobians) for the gradient of
the variance-covariance of beta with respect to the variance parameters,
where beta are the mean-value parameters available in \code{fixef(object)}.}

\item{\code{vcov_beta}}{a numeric matrix holding the asymptotic variance-covariance
matrix of the fixed-effect regression parameters (beta).}

\item{\code{sigma}}{the residual standard deviation.}
}}

\seealso{
\code{\link[lme4]{lmer}} and \code{\link[lme4]{merMod}}
}
\author{
Rune Haubo B. Christensen
}
