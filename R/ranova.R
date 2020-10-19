#############################################################################
#    Copyright (c) 2013-2020 Alexandra Kuznetsova, Per Bruun Brockhoff, and
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
# ranova.R - random effects ANOVA table

# ------- Contents: --------
#
# --- exported function: ---
#
# ranova
# rand
#
# --- utility functions: ---
#
# rm_complete_terms
# get_lm_call
# get_newforms
# get_logLik
# mk_LRtab
# has_ranef
# has_terms
# get_lhs
# get_rhs
#

##############################################
######## ranova(); rand()
##############################################
#' ANOVA-Like Table for Random-Effects
#'
#' Compute an ANOVA-like table with tests of random-effect terms in the model.
#' Each random-effect term is reduced or removed and likelihood ratio tests of
#' model reductions are presented in a form similar to that of
#' \code{\link[=drop1.lmerModLmerTest]{drop1}}.
#' \code{rand} is an alias for \code{ranova}.
#'
#' If the model is fitted with REML the tests are REML-likelihood ratio tests.
#'
#' A random-effect term of the form \code{(f1 + f2 | gr)} is reduced to
#' terms of the form \code{(f2 | gr)} and \code{(f1 | gr)} and these reduced
#' models are compared to the original model.
#' If \code{reduce.terms} is \code{FALSE} \code{(f1 + f2 | gr)} is removed
#' instead.
#'
#' A random-effect term of the form \code{(f1 | gr)} is reduced to \code{(1 | gr)}
#' (unless \code{reduce.terms} is \code{FALSE}).
#'
#' A random-effect term of the form \code{(1 | gr)} is not reduced but
#' simply removed.
#'
#' A random-effect term of the form \code{(0 + f1 | gr)} or \code{(-1 + f1 | gr)}
#' is reduced (if \code{reduce.terms = TRUE}) to \code{(1 | gr)}.
#'
#' A random-effect term of the form \code{(1 | gr1/gr2)} is automatically
#' expanded to two terms: \code{(1 | gr2:gr1)} and \code{(1 | gr1)} using
#' \code{\link[lme4]{findbars}}.
#'
#' In this exposition it is immaterial whether \code{f1} and \code{f2} are
#' factors or continuous variables.
#'
#' @note Note that \code{anova} can be used to compare two models and will often
#' be able to produce the same tests as \code{ranova}. This is, however, not always the
#' case as illustrated in the examples.
#'
#' @section Warning:
#' In certain cases tests of non-nested models may be generated. An example
#' is when \code{(0 + poly(x, 2) | gr)} is reduced (the default) to \code{(1 | gr)}.
#' To our best knowledge non-nested model comparisons are only generated in
#' cases which are statistical nonsense anyway (such as in this example where
#' the random intercept is suppressed).
#'
#'
#' @param model a linear mixed effect model fitted with \code{lmer()}
#' (inheriting from class \code{lmerMod}).
#' @param reduce.terms if \code{TRUE} (default) random-effect terms are
#' reduced (if possible). If \code{FALSE} random-effect terms are simply
#' removed.
#' @param ... currently ignored
#'
#' @return an ANOVA-like table with single term deletions of random-effects
#' inheriting from class \code{anova} and \code{data.frame} with the columns:
#' \item{npar}{number of model parameters.}
#' \item{logLik}{the log-likelihood for the model. Note that this is the
#' REML-logLik if the model is fitted with REML.}
#' \item{AIC}{the AIC for the model evaluated as \code{-2*(logLik - npar)}.
#' Smaller is better.}
#' \item{LRT}{the likelihood ratio test statistic; twice the difference in
#' log-likelihood, which is asymptotically chi-square distributed.}
#' \item{Df}{degrees of freedom for the likelihood ratio test: the difference in
#' number of model parameters.}
#' \item{Pr(>Chisq)}{the p-value.}
#' @export
#' @author Rune Haubo B. Christensen and Alexandra Kuznetsova
#'
#' @seealso \code{\link[=drop1.lmerModLmerTest]{drop1}} for tests of marginal
#' fixed-effect terms and
#' \code{\link{anova}} for usual anova tables for fixed-effect terms.
#' @importFrom stats formula nobs update
#' @importFrom lme4 getME findbars nobars
#'
#' @examples
#'
#' # Test reduction of (Days | Subject) to (1 | Subject):
#' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#' ranova(fm1) # 2 df test
#'
#' # This test can also be achieved with anova():
#' fm2 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
#' anova(fm1, fm2, refit=FALSE)
#'
#' # Illustrate reduce.test argument:
#' # Test removal of (Days | Subject):
#' ranova(fm1, reduce.terms = FALSE) # 3 df test
#'
#' # The likelihood ratio test statistic is in this case:
#' fm3 <- lm(Reaction ~ Days, sleepstudy)
#' 2*c(logLik(fm1, REML=TRUE) - logLik(fm3, REML=TRUE)) # LRT
#'
#' # anova() is not always able to perform the same tests as ranova(),
#' # for example:
#' anova(fm1, fm3, refit=FALSE) # compares REML with ML and should not be used
#' anova(fm1, fm3, refit=TRUE) # is a test of ML fits and not what we seek
#'
#' # Also note that the lmer-fit needs to come first - not an lm-fit:
#' # anova(fm3, fm1) # does not work and gives an error
#'
#' # ranova() may not generate all relevant test:
#' # For the following model ranova() indicates that we should not reduce
#' # (TVset | Assessor):
#' fm <- lmer(Coloursaturation ~ TVset * Picture + (TVset | Assessor), data=TVbo)
#' ranova(fm)
#' # However, a more appropriate model is:
#' fm2 <- lmer(Coloursaturation ~ TVset * Picture + (1 | TVset:Assessor), data=TVbo)
#' anova(fm, fm2, refit=FALSE)
#' # fm and fm2 has essentially the same fit to data but fm uses 5 parameters
#' # more than fm.
#'
ranova <- function(model, reduce.terms=TRUE, ...) {
  if(!inherits(model, "lmerMod"))
    stop("'model' should be an lmer-fit: \"inherits(model, 'lmerMod')\" is not TRUE")
  isREML <- getME(model, "is_REML")
  nobs_model <- nobs(model)
  orig_form <- formula(model)
  orig_rhs <- orig_form[[length(orig_form)]]
  if(!has_ranef(orig_rhs))
    stop("Model should have at least one random-effects term")

  # Reconstruct formula - needed for terms like (1 | g1 / g2):
  fe_rhs <- deparse2(nobars(orig_rhs))
  reforms <- lapply(findbars(orig_rhs), deparse2) # random-effect forms
  re_rhs <- lapply(reforms, function(rf) paste0("(", rf, ")"))
  full_rhs <- paste(c(list(fe_rhs), re_rhs), collapse=" + ")
  full_form <- update(orig_form, paste0(". ~", full_rhs))

  # Compute new model formulae with reduced ranef formulae:
  new_forms <- if(!reduce.terms)
    rm_complete_terms(reforms, full_form) else
      unlist(lapply(reforms, get_newforms, full_formula=full_form))
  ll <- get_logLik(model) # store df and logLik

  for(nform in new_forms) { # For each new formula. nform <- new_forms[[1]]
    newfit <- if(!has_ranef(nform)) { # If no random effects: fit with lm
      lm_call <- get_lm_call(model, nform)
      eval.parent(as.call(lm_call))
    } else eval.parent(update(model, formula=nform))
  # } else eval.parent(update(model, formula=nform, ...))
    # Check that models were fit to the same number of observations:
    nobs_newfit <- nobs(newfit)
    if(all(is.finite(c(nobs_model, nobs_newfit))) && nobs_newfit != nobs_model)
      stop("number of rows in use has changed: remove missing values?")
    ll <- rbind(ll, get_logLik(newfit, REML=isREML)) # store df and logLik
  }

  # Collect information in ANOVA table and return:
  aov <- mk_LRtab(ll)
  rownames(aov) <- c("<none>", names(new_forms))
  head <- c("ANOVA-like table for random-effects: Single term deletions",
            "\nModel:", deparse2(full_form))
  attr(aov, "formulae") <- new_forms
  structure(aov, heading = head, class = c("anova", "data.frame"))
}

#' @rdname ranova
#' @export
rand <- ranova

##############################################
######## ranova utility functions below
##############################################

#' Remove Terms from Formula
#'
#' Remove fixef or ranef terms from formula, return a list of modified formulae
#' with environment restored to that of the original formula.
#'
#' @param terms character vector (or list) of terms to remove from
#' \code{full_formula}
#' @param full_formula formula
#' @param random if \code{TRUE} names of the return list have parentheses around
#' them.
#'
#' @importFrom stats update.formula
#' @keywords internal
rm_complete_terms <- function(terms, full_formula, random=TRUE) {
  # Remove random-effect formula terms from original model formula (full_formula)
  forms <- lapply(terms, function(reform) {
    form <- update.formula(full_formula, paste0("~.- (", reform, ")"))
    environment(form) <- environment(full_formula)
    form
  })
  names(forms) <- if(!random) terms else
    sapply(terms, function(form) paste0("(", form, ")"))
  forms
}


#' @importFrom stats getCall
get_lm_call <- function(object, formula) {
  # object: lmerMod object
  # formula: model formula without random effects
  Call <- as.list(getCall(object))
  notkeep <- c("control", "start", "verbose", "devFunOnly", "REML")
  Call <- Call[!names(Call) %in% notkeep]
  Call$formula <- formula
  Call[[1]] <- as.name("lm")
  Call
}

#' @importFrom stats update.formula drop.scope
get_newforms <- function(form, full_formula) {
  # Update full_formula by reducing the random-effect structure of 'form'
  #
  # form: a deparse'd random-effect formula term
  # full_formula: the original model formula with lhs, fixed and random terms
  #
  rhs <- get_rhs(form) # rhs of random term: (lhs | rhs)
  lhs <- get_lhs(form) # lhs of random term: (lhs | rhs)
  scope <- drop.scope(lhs) # Detemine terms to drop from lhs
  # Determine list of updates to 'form'
  update_forms <- if(!has_terms(lhs) || length(scope) == 0L) {# length(scope) >= 1
    # Remove entire re-term if lhs is '1':
    setNames(list(paste0("~.- (", form, ")")), paste0("(", form, ")"))
  } else {
    # Drop terms from lhs of random term:
    ll <- lapply(scope, function(scp) { # scp <- scope
      # If there are no other terms in lhs than scp set new_lhs to just ~1:
      new_lhs <- if(setequal(attr(terms(lhs), "term.labels"), scp)) "1" else {
        tmp <- deparse2(update.formula(lhs, paste("~.-", scp)))
        gsub("~", "", tmp, fixed=TRUE)
      }
      new_form <- paste0("(", new_lhs, " | ", rhs, ")")
      paste0("~.- (", form, ")", " + ", new_form)
    })
    names(ll) <- paste(scope, paste0("in (", form, ")"))
    ll
  }
  # Update original formula 'full_formula' with update_forms and return:
  lapply(update_forms, function(upd) {
    form <- update.formula(full_formula, upd)
    environment(form) <- environment(full_formula)
    form
  })
}

#' @importFrom stats logLik
get_logLik <- function(object, ...) {
  # Extract data.frame with "df" and "logLik" values from object.
  ll <- logLik(object, ...)
  data.frame("Df"=attr(ll, "df"), "logLik"=c(ll))
}

#' @importFrom stats pchisq
mk_LRtab <- function(x) {
  # Compute drop1-table with LR-tests
  # x: a 2-col data.frame with "Df" and "logLik"; 1st row is the full model
  chisq_pval <- function(q, df, ...) pchisq(q=q, ifelse(df > 0, df, NA_real_), ...)
  stopifnot(is.data.frame(x), colnames(x) == c("Df", "logLik"))
  res <- data.frame("npar" = x[, "Df"],
                    "logLik" = x[, "logLik"],
                    "AIC" = -2*x[, "logLik"] + 2*x[, "Df"],
                    "LRT" = NA_real_,
                    "Df"  = NA_real_,
                    "Pr(>Chisq)" = NA_real_, check.names = FALSE)
  if(nrow(x) >= 2) {
    res[-1, "LRT"] <- 2*(x[1, "logLik"] - x[-1, "logLik"])
    res[-1, "Df"] <- x[1, "Df"] - x[-1, "Df"]
    res[-1, "Pr(>Chisq)"] <-
      chisq_pval(res[-1, "LRT"], res[-1, "Df"], lower.tail=FALSE)
  }
  rownames(res) <- rownames(x)
  res
}

has_ranef <- function(form) {
  # Determine if formula 'form' contain random effect terms.
  if(is.character(form)) form <- deparse2(form)
  length(grep("|", form, fixed=TRUE)) > 0
}

has_terms <- function(form) {
  # Determine if formula 'form' contain any terms beyond intercept.
  length(attr(terms(form), "term.labels")) > 0
}

get_lhs <- function(ranef_term) {
  # Extract lhs in (lhs | rhs)
  if(!is.character(ranef_term)) ranef_term <- deparse2(ranef_term)
  lhs <- trimws(gsub("\\|.*$", "", ranef_term))
  form <- as.formula(paste0("~", lhs))
  form
  ## Add "1" for intercept if is suppressed:
  # FIXME: Only if there no other terms in lhs?
  # if(attr(terms(form), "intercept") == 1) form else
  #   as.formula(paste0("~1 + ", lhs))
}

get_rhs <- function(ranef_term) {
  # Extract rhs in (lhs | rhs)
  if(!is.character(ranef_term)) ranef_term <- deparse2(ranef_term)
  trimws(gsub("^.*\\|", "", ranef_term))
}
