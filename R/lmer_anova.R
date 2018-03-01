# lmer_anova.R - anova method for lmerModLmerTest objects

# anova - lmerModLmerTest
# single_anova
# show_tests

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
                   ddf=c("Satterthwaite", "Kenward-Roger", "lme4")) {
            dots <- list(...)
            models <- if(length(dots))
              sapply(dots, is, "lmerModLmerTest") | sapply(dots, is, "merMod") |
              sapply(dots, is, "lm") else logical(0)
            if(any(models)) return(callNextMethod()) # return(anova(as(object, "lmerMod"), ...))
            # Note: Need 'callNextMethod' here to get printing from anova.merMod right.
            ddf <- match.arg(ddf)
            # type <- match.arg(type) # not actually needed
            if(ddf=="lme4") return(anova(as(object, "lmerMod"), ...)) # return(callNextMethod())
            # FIXME: Warn that 'type' is ignored when ddf="lme4"
            single_anova(object=object, type=type, ddf=ddf)
          })


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
                         type = c("I", "II", "III", "1", "2", "3", "yates", "marginal", "2b", "3b", "3c"),
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
  } else if(type == "3b" || type == "yates") {
    get_contrasts_type3b(object)
  } else if(type == "3c") {
    get_contrasts_type3old(object)
  } else if(type == "marginal") {
    get_contrasts_marginal(object)
  } else {
    stop("'type' not recognized")
  }
  # Get F-test for each term and collect in table:
  table <- rbindall(lapply(L_list, contestMD, model=object, ddf=ddf))
  # Format ANOVA table and return:
  rownames(table) <- names(L_list)
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
######## show_tests()
##############################################
#' Show Hypothesis Tests in ANOVA Tables
#'
#' Extracts hypothesis matrices for terms in ANOVA tables detailing exactly which
#' functions of the parameters are being tested in anova tables.
#'
#' @param object an anova table with a \code{"hypotheses"} attribute.
#' @param fractions Display entries in the hypothesis matrices as fractions?
#' @param names if \code{FALSE} column and row names of the hypothesis matrices
#' are suppressed.
#'
#' @return a list of hypothesis matrices.
#' @importFrom MASS fractions
#' @author Rune Haubo B. Christensen
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
show_tests <- function(object, fractions=FALSE, names=TRUE) {
  tests <- attr(object, "hypotheses")
  # FIXME: Maybe this should be a generic with a method for anova objects?
  if(is.null(tests))
    stop("'object' does not have an 'hypotheses' attribute")
  if(fractions) tests <- lapply(tests, MASS::fractions)
  if(names) tests else lapply(tests, unname)
}
