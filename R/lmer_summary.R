#############################################################################
#    Copyright (c) 2013-2018 Alexandra Kuznetsova, Per Bruun Brockhoff, and
#    Rune Haubo Bojesen Christensen
#
#    This file is part of the lmerTest package for R (*lmerTest*)
#
#    *lmerTest* is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    *lmerTest* is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at
#    <https://www.r-project.org/Licenses/> and/or
#    <http://www.gnu.org/licenses/>.
#############################################################################
#
# lmer_summary.R - summary method for lmerModLmerTest objects

# ------- Contents: --------
#
# summary.lmerModLmerTest
#
# --- utility functions: ---
#
# get_coefmat
#

#' @include lmer.R
NULL

##############################################
######## summary method for lmerModLmerTest
##############################################
#' Summary Method for Linear Mixed Models
#'
#' Summaries of Linear Mixed Models with coefficient tables including t-tests
#' and p-values using Satterthwaites's or Kenward-Roger's methods for
#' degrees-of-freedom and t-statistics.
#'
#' The returned object is of class
#' \code{c("summary.lmerModLmerTest", "summary.merMod")} utilizing \code{print},
#' \code{coef} and other methods defined for \code{summary.merMod} objects.
#' The \code{"Kenward-Roger"} method use methods from the \pkg{pbkrtest} package internally
#' to compute t-statistics and associated degrees-of-freedom.
#'
#' @param object an lmerModLmerTest object.
#' @param ddf the method for computing the degrees of freedom and
#' t-statistics. \code{ddf="Satterthwaite"} (default) uses Satterthwaite's method;
#' \code{ddf="Kenward-Roger"} uses Kenward-Roger's method,
#' \code{ddf = "lme4"} returns the lme4-summary i.e., using the summary
#' method for \code{lmerMod} objects as defined in the \pkg{lme4}-package and
#' ignores the \code{type} argument. Partial matching is allowed.
#' @param ... additional arguments passed on to \code{lme4::summary.merMod}
#'
#' @return A summary object with a coefficient table (a \code{matrix}) including
#' t-values and p-values. The coefficient table can be extracted with
#' \code{coef(summary(<my-model>))}.
#'
#' @seealso \code{\link{contest1D}} for one degree-of-freedom contrast tests
#' and \code{\link[pbkrtest]{KRmodcomp}} for Kenward-Roger F-tests.
#' @author Rune Haubo B. Christensen and Alexandra Kuznetsova
#' @export
#' @importFrom methods as signature
#'
#' @examples
#'
#' # Fit example model:
#' data("sleepstudy", package="lme4")
#' fm <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
#'
#' # Get model summary:
#' summary(fm) # Satterthwaite df and t-tests
#'
#' # Extract coefficient table:
#' coef(summary(fm))
#'
#' # Use the Kenward-Roger method
#' if(requireNamespace("pbkrtest", quietly = TRUE))
#'   summary(fm, ddf="Kenward-Roger")
#'
#' # The lme4-summary table:
#' summary(fm, ddf="lme4") # same as summary(as(fm, "lmerMod"))
#'
#' \dontshow{
#'   # Check that summaries are as expected:
#'   summ_fm <- coef(summary(fm))
#'   summ_fm_lme4 <- coef(summary(fm, ddf="lme4"))
#'   stopifnot(
#'     all(colnames(summ_fm) == c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")),
#'     all(colnames(summ_fm_lme4) == c("Estimate", "Std. Error", "t value")),
#'     all(!(is.na(summ_fm))),
#'     all(!(is.na(summ_fm_lme4)))
#'   )
#'  if(requireNamespace("pbkrtest", quietly = TRUE) && getRversion() >= "3.3.3") {
#'     summ_fm_kr <- coef(summary(fm, ddf="Kenward-Roger"))
#'      stopifnot(
#'        all(colnames(summ_fm_kr) == c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")),
#'        all(!(is.na(summ_fm_kr)))
#'     )
#'  }
#' }
summary.lmerModLmerTest <- function(object, ...,
                                    ddf=c("Satterthwaite", "Kenward-Roger", "lme4")) {
  ddf <- match.arg(ddf)
  if(!inherits(object, "lmerModLmerTest") && !inherits(object, "lmerMod")) {
    stop("Cannot compute summary for objects of class: ",
         paste(class(object), collapse = ", "))
  }
  if(!inherits(object, "lmerModLmerTest") && inherits(object, "lmerMod")) {
    message("Coercing object to class 'lmerModLmerTest'")
    object <- as_lmerModLmerTest(object)
    if(!inherits(object, "lmerModLmerTest")) {
      warning("Failed to coerce object to class 'lmerModLmerTest'")
      return(summary(object))
    }
  }
  summ <- summary(as(object, "lmerMod"), ...)
  if(ddf == "lme4") return(summ)
  summ$coefficients <- get_coefmat(object, ddf=ddf)
  ddf_nm <- switch(ddf, "Satterthwaite" = "Satterthwaite's",
                   "Kenward-Roger" = "Kenward-Roger's")
  summ$objClass <- class(object) # Used by lme4:::print.summary.lmerMod
  summ$methTitle <- paste0(summ$methTitle, ". t-tests use ", ddf_nm, " method")
  class(summ) <- c("summary.lmerModLmerTest", class(summ))
  summ
}


##############################################
######## get_coefmat
##############################################
#' @importFrom lme4 fixef
get_coefmat <- function(model, ddf=c("Satterthwaite", "Kenward-Roger")) {
  ddf <- match.arg(ddf)
  p <- length(fixef(model))
  if(p < 1)
    return(as.matrix(contest1D(model, numeric(0L), ddf=ddf)))
  Lmat <- diag(p)
  tab <- rbindall(lapply(1:p, function(i) contest1D(model, Lmat[i, ], ddf=ddf)))
  rownames(tab) <- names(fixef(model))
  as.matrix(tab)
}
