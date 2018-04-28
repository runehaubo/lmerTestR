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
# lmer_anova.R - anova method for lmerModLmerTest objects

# ------- Contents: --------
#
# --- Generics: ---
#
# show_tests
#
# --- methods: ---
#
# anova.lmerModLmerTest
#
# show_tests.default
# show_tests.anova
#
# --- other exported function: ---
#
# show_contrasts
#
# --- utility functions: ---
#
# single_anova
#

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
#' \code{\link{lmer}} from the \pkg{lmerTest}-package.
#'
#' The \code{"Kenward-Roger"} method calls \code{pbkrtest::KRmodcomp} internally and
#' reports scaled F-statistics and associated denominator degrees-of-freedom.
#'
#' @param object an \code{lmerModLmerTest} object; the result of \code{lmer()}
#' after loading the \pkg{lmerTest}-package.
#' @param ... potentially additional \code{lmer} or \code{lm} model objects for
#' comparison of models in which case \code{type} and \code{ddf} arguments are
#' ignored.
#' @param type the type of ANOVA table requested (using SAS terminology)
#' with Type I being the familiar sequential ANOVA table.
#' @param ddf the method for computing the denominator degrees of freedom and
#' F-statistics. \code{ddf="Satterthwaite"} (default) uses Satterthwaite's method;
#' \code{ddf="Kenward-Roger"} uses Kenward-Roger's method,
#' \code{ddf = "lme4"} returns the lme4-anova table, i.e., using the anova
#' method for \code{lmerMod} objects as defined in the \pkg{lme4}-package and
#' ignores the \code{type} argument. Partial matching is allowed.
#'
#' @return an ANOVA table
#' @seealso \code{\link{contestMD}} for multi degree-of-freedom contrast tests
#' and \code{\link[pbkrtest]{KRmodcomp}} for the \code{"Kenward-Roger"} method.
#' @author Rune Haubo B. Christensen and Alexandra Kuznetsova
#' @importFrom methods is callNextMethod
#' @importFrom stats anova
#' @export
#'
#' @examples
#'
#' data("sleepstudy", package="lme4")
#' m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' anova(m) # with p-values from F-tests using Satterthwaite's denominator df
#' anova(m, ddf="lme4") # no p-values
#'
#' # Use the Kenward-Roger method
#' if(requireNamespace("pbkrtest", quietly = TRUE))
#'   anova(m, ddf="Kenward-Roger")
#'
#' \dontshow{
#'   an1 <- anova(m) # with p-values from F-tests using Satterthwaite's denominator df
#'   an2 <- anova(m, ddf="lme4")
#'   stopifnot(
#'     all(colnames(an1) == c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F value", "Pr(>F)")),
#'     !"Pr(>F)" %in% colnames(an2),
#'     all(!is.na(an1)),
#'     all(!is.na(an2))
#'   )
#' }
anova.lmerModLmerTest <- function(object, ..., type = c("III", "II", "I", "3", "2", "1"),
                                  ddf=c("Satterthwaite", "Kenward-Roger", "lme4")) {
  if(!inherits(object, "lmerModLmerTest") && !inherits(object, "lmerMod")) {
    stop("'object' of class: ", paste(class(object), collapse = ", "),
         ". Expecting object of class 'lmerModLmerTest'")
  }
  if(!inherits(object, "lmerModLmerTest") && inherits(object, "lmerMod")) {
    message("Coercing object to class 'lmerModLmerTest'")
    object <- as_lmerModLmerTest(object)
    if(!inherits(object, "lmerModLmerTest")) {
      warning("Failed to coerce object to class 'lmerModLmerTest'")
      return(NextMethod())
    }
  }
  dots <- list(...)
  models <- if(length(dots))
    sapply(dots, is, "lmerModLmerTest") | sapply(dots, is, "merMod") |
    sapply(dots, is, "lm") else logical(0)
  if(any(models)) return(NextMethod()) # return(anova(as(object, "lmerMod"), ...))
  # Note: Need 'NextMethod' here to get printing from anova.merMod right.
  ddf <- match.arg(ddf)
  # Commented since we need to pass 'hidden' type options to single_anova
  # type <- match.arg(type)
  if(ddf=="lme4") return(anova(as(object, "lmerMod"), ...)) # return(NextMethod())
  # FIXME: Warn that 'type' is ignored when ddf="lme4"?
  single_anova(object=object, type=type, ddf=ddf)
}


# #' @export
# #' @keywords internal
# anova <- function(object, ...) UseMethod("anova")


##############################################
######## single_anova()
##############################################
#' ANOVA Tables for Linear Mixed Models
#'
#' @param object an \code{lmerModLmerTest} object; the result of \code{lmer()}
#' after loading the \pkg{lmerTest}-package.
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
single_anova <- function(object,
                         type = c("III", "II", "I", "3", "2", "1", "yates", "marginal", "2b"),
                         ddf=c("Satterthwaite", "Kenward-Roger")) {
  if(!inherits(object, "lmerModLmerTest"))
    warning("calling single_anova(<fake-lmerModLmerTest-object>) ...")
  type <- type[1L]
  if(!is.character(type)) type <- as.character(type)
  type <- match.arg(type)
  if(type %in% c("I", "II", "III"))
    type <- as.character(as.integer(as.roman(type)))
  ddf <- match.arg(ddf)
  # Get list of contrast matrices (L) - one for each model term:
  L_list <- if(type == "1") {
    get_contrasts_type1(object)
  } else if(type == "2") {
    get_contrasts_type2_unfolded(object)
  } else if(type == "2b") {
    get_contrasts_type2(object)
  } else if(type == "3") {
    get_contrasts_type3(object)
  } else if(type == "yates") {
    get_contrasts_yates(object)
  } else if(type == "marginal") {
    get_contrasts_marginal(object)
  } else {
    stop("'type' not recognized")
  }
  # Get F-test for each term and collect in table:
  table <- rbindall(lapply(L_list, function(L) contestMD(object, L, ddf=ddf)))
  # Format ANOVA table and return:
  if(length(nm <- setdiff(names(L_list), rownames(table)))) {
    tab <- array(NA_real_, dim=c(length(nm), 6L),
                 dimnames = list(nm, colnames(table)))
    table <- rbind(table, tab)
  }
  method <- switch(ddf, "Satterthwaite" = "Satterthwaite's",
                   "Kenward-Roger" = "Kenward-Roger's")
  # Format 'type':
  type <- if(type == "marginal") {
    "Marginal"
  } else if (type == "yates" || type == "3b") {
    "Yates"
  } else if(grepl("b|c", type)) {
    alph <- gsub("[0-9]", "", type)
    paste0("Type ", as.roman(as.integer(gsub("b|c", "", type))), alph)
  } else paste("Type", as.roman(as.integer(type)))
  attr(table, "heading") <-
    paste(type, "Analysis of Variance Table", "with", method, "method")
  attr(table, "hypotheses") <- L_list
  class(table) <- c("anova", "data.frame")
  table
}

##############################################
######## show_tests.anova()
##############################################
#' Show Hypothesis Tests in ANOVA Tables
#'
#' Extracts hypothesis matrices for terms in ANOVA tables detailing exactly which
#' functions of the parameters are being tested in anova tables.
#'
#' @param object an anova table with a \code{"hypotheses"} attribute.
#' @param fractions display entries in the hypothesis matrices as fractions?
#' @param names if \code{FALSE} column and row names of the hypothesis matrices
#' are suppressed.
#' @param ... currently not used.
#'
#' @return a list of hypothesis matrices.
#' @importFrom MASS fractions
#' @author Rune Haubo B. Christensen
#' @seealso \code{\link[=show_tests.ls_means]{show_tests}} for \code{ls_means}
#' objects.
#' @export
#'
#' @examples
#'
#' # Fit basic model to the 'cake' data:
#' data("cake", package="lme4")
#' fm1 <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
#'
#' # Type 3 anova table:
#' (an <- anova(fm1, type="3"))
#'
#' # Display tests/hypotheses for type 1, 2, and 3 ANOVA tables:
#' # (and illustrate effects of 'fractions' and 'names' arguments)
#' show_tests(anova(fm1, type="1"))
#' show_tests(anova(fm1, type="2"), fractions=TRUE, names=FALSE)
#' show_tests(an, fractions=TRUE)
#'
show_tests.anova <- function(object, fractions=FALSE, names=TRUE, ...)
  NextMethod() # use default method

##############################################
######## show_tests()
##############################################
#' Show Tests Generic Function and Default Method
#'
#' @param object a suitable object with an \code{"hypotheses"} attribute, e.g. an
#' anova table or an \code{ls_means} table as defined in \pkg{lmerTest}.
#' @param ... parsed on to methods; currently not used in the default method.
#'
#' @export
#' @author Rune Haubo B. Christensen
#' @seealso \code{\link{show_tests.anova}} and \code{\link{show_tests.ls_means}}
#' @keywords internal
show_tests <- function(object, ...) UseMethod("show_tests")


##############################################
######## show_tests.default()
##############################################
#' @rdname show_tests
#'
#' @param fractions display entries in the hypothesis matrices as fractions?
#' @param names if \code{FALSE} column and row names of the hypothesis matrices
#' are suppressed.
#' @export
#' @author Rune Haubo B. Christensen
#' @keywords internal
show_tests.default <- function(object, fractions=FALSE, names=TRUE, ...) {
  tests <- attr(object, "hypotheses")
  # FIXME: Maybe this should be a generic with a method for anova objects?
  if(is.null(tests))
    stop("'object' does not have an 'hypotheses' attribute")
  if(fractions) tests <- lapply(tests, MASS::fractions)
  if(names) tests else lapply(tests, unname)
}
