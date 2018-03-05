# lmer_summary.R - summary method for lmerModLmerTest objects

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
#' summary(fm, ddf="Kenward-Roger")
#'
#' # The lme4-summary table:
#' summary(fm, ddf="lme4") # same as summary(as(fm, "lmerMod"))
#'
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

#' @importFrom lme4 fixef
get_coefmat <- function(model, ddf=c("Satterthwaite", "Kenward-Roger")) {
  ddf <- match.arg(ddf)
  p <- length(fixef(model))
  if(p < 1)
    return(as.matrix(contest1D(numeric(0L), model=model, ddf=ddf)))
  Lmat <- diag(p)
  tab <- rbindall(lapply(1:p, function(i)
    contest1D(Lmat[i, ], model=model, ddf=ddf)))
  rownames(tab) <- names(fixef(model))
  as.matrix(tab)
}
