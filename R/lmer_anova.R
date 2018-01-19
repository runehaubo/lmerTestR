# lmer_anova.R - anova method for lmerModLmerTest objects

#' @include lmer.R
NULL

##############################################
######## anova method for lmerModLmerTest
##############################################
#' ANOVA Tables for Linear Mixed Models
#'
#' ANOVA table with F-tests and p-values using Satterthwaite's or
#' Kenward-Roger's method for denominator degrees-of-freedom and F-statistic.
#' Models should be fitted with
#' \code{\link{lmer}} from the \pkg{lmerTestR}-package.
#'
#' The \code{"KR"} method calls \code{pbkrtest::KRmodcomp} internally and
#' reports scaled F-statistics and associated denominator degrees-of-freedom.
#'
#' @param object an \code{lmerModLmerTest} object; the result of \code{lmer()}
#' after loading the \pkg{lmerTestR}-package.
#' @param ... potentially additional \code{lmer} or \code{lm} model objects for
#' comparison of models in which case \code{type} and \code{ddf} arguments are
#' ignored.
#' @param type the type of ANOVA table requested (using SAS terminology)
#' with Type I being the familiar sequential ANOVA table.
#' @param ddf the method for computing the denominator degrees of freedom and
#' F-statistics. \code{ddf="Satterthwaite"} (default) uses Satterthwaite's method;
#' \code{ddf="KR"} uses Kenward-Roger's method,
#' \code{ddf = "lme4"} returns the lme4-anova table, i.e., using the anova
#' method for \code{lmerMod} objects as defined in the \pkg{lme4}-package and
#' ignores the \code{type} argument. Partial matching is allowed.
#'
#' @return an ANOVA table
#' @seealso \code{\link{contestMD}} for multi degree-of-freedom contrast tests
#' and \code{\link[pbkrtest]{KRmodcomp}} for the \code{"KR"} method.
#' @author Rune Haubo B. Christensen
#' @importFrom methods is callNextMethod
#' @export
#'
#' @examples
#'
#' data("sleepstudy", package="lme4")
#' m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' anova(m) # with p-values from F-tests using Satterthwaite's denominator df
#' anova(m, ddf="lme4")
#'
setMethod("anova",
          signature(object="lmerModLmerTest"),
          function(object, ..., type = c("I", "II", "III", "1", "2", "3"),
                   ddf=c("Satterthwaite", "KR", "lme4")) {
            dots <- list(...)
            models <- if(length(dots))
              sapply(dots, is, "merModLmerTest") | sapply(dots, is, "merMod") |
              sapply(dots, is, "lm") else logical(0)
            if(any(models)) return(callNextMethod())
            ddf <- match.arg(ddf)
            # type <- match.arg(type) # not actually needed
            if(ddf=="lme4") return(callNextMethod())
            single_anova(object=object, type=type, ddf=ddf)
          })


##############################################
######## anova function for single models
##############################################
#' ANOVA Tables for Linear Mixed Models
#'
#' @param object an \code{lmerModLmerTest} object; the result of \code{lmer()}
#' after loading the \pkg{lmerTestR}-package.
#' @param type the type of ANOVA table requested (using the SAS terminology for
#' these) with Type I being the familiar sequential ANOVA table.
#' @param ddf method for computing denominator degrees of freedom.
#'
#' @return an ANOVA table
#' @importFrom utils as.roman
#' @importFrom stats model.matrix terms formula
#' @author Rune Haubo B. Christensen
#'
#' @keywords internal
single_anova <- function(object, type = c("I", "II", "III", "1", "2", "3"),
                         ddf=c("Satterthwaite", "KR")) {
  if(!inherits(object, "lmerModLmerTest"))
    warning("calling anova(<fake-lmerModLmerTest-object>) ...")
  if(!is.character(type)) type <- as.character(type)
  type <- as.integer(as.roman(match.arg(type)))
  if(type > 1L) {
    warning("Type II and III anova tables are not yet implemented; returning type I")
    type <- 1L
  }
  ddf <- match.arg(ddf)
  # Get list of contrast matrices (L) - one for each model term:
  L_list <- get_contrasts_type1(model.matrix(object), terms(object))
  # Get F-test for each term and collect in table:
  table <- rbindall(lapply(L_list, contestMD, model=object, ddf=ddf))
  # Format ANOVA table and return:
  rownames(table) <- names(L_list)
  response_name <- deparse(formula(object)[[2L]], width.cutoff = 500L)
  method <- switch(ddf, "Satterthwaite" = "Satterthwaite's",
                   "KR" = "Kenward-Roger's")
  attr(table, "heading") <-
    paste("Type", as.roman(type), "Analysis of Variance Table",
          "with", method, "method")
  class(table) <- c("anova", "data.frame")
  table
}

###############################################################################
######## Type I anova functions below
###############################################################################

##############################################
######## get_contrasts_type1
##############################################

#' Type I ANOVA table contrasts
#'
#' @param X a design matrix - usually from \code{model.matrix}. \code{X} must
#' have an \code{assign} attribute.
#' @param terms a terms object.
#' @param keep_intercept defaults to \code{FALSE}. If \code{TRUE} a contrast
#' for the intercept is included in the returned list of contrast matrices.
#'
#' @return List of contrast matrices - one contrast matrix for each model term.
#' @importFrom stats setNames
#' @author Rune Haubo B. Christensen
#'
#' @keywords internal
get_contrasts_type1 <- function(X, terms, keep_intercept = FALSE) {
  p <- ncol(X)
  if(p == 0L) return(list(matrix(numeric(0L), nrow=0L))) # no fixef
  if(p == 1L && attr(terms, "intercept")) # intercept-only model
    return(list(matrix(numeric(0L), ncol=1L)))
  # Compute 'normalized' doolittle factorization of XtX:
  # L <- t(doolittle(crossprod(X)))
  L <- if(p == 1L) matrix(1L) else normalized_doolittle(crossprod(X))
  # Determine which rows of L belong to which term:
  asgn <- attr(X, "assign")
  stopifnot(!is.null(asgn))
  term_labels <- attr(terms, "term.labels")
  term_names <- c("(Intercept)", term_labels)
  term_names <- term_names[1 + unique(asgn)] # order appropriately
  # Compute list of row indicators for L matrix:
  ind.list <- setNames(split(1L:p, asgn), nm=term_names)
  ind.list <- ind.list[term_labels] # rm intercept if present
  # if(length(ind.list) == 0L) return(list(matrix(numeric(0L), nrow=0L)))
  lapply(ind.list, function(rows) L[rows, , drop=FALSE])
}

