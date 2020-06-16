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
# step.R - implementation of backward elimination for lmerModLmerTest objects

# ------- Contents: --------
#
# --- Generics: ---
#
# step
# get_model
#
# --- methods: ---
#
# step.lmerModLmerTest
# step.default
# get_model.step_list
# print.step_list
# plot.step_list
#
# --- other exported function: ---
#
# --- utility functions: ---
#
# ran_redTable
# fix_redTable
# reduce_random
# ranova_lm
# reduce_fixed
#


##############################################
######## step()
##############################################
#' Generic Step Function
#'
#' Generic step function with default method \code{stats::step}. This
#' construction ensures that \code{stats::step} still works on \code{lm}
#' objects etc. after loading the \pkg{lmerTest} package.
#'
#' @param object a model object.
#' @param ... currently not used.
#'
#' @author Rune Haubo B. Christensen
#' @seealso \code{\link[=step.lmerModLmerTest]{step}}
#' @export
#' @keywords internal
step <- function(object, ...) UseMethod("step")


##############################################
######## step.default()
##############################################
#' @rdname step
#' @export
#' @keywords internal
step.default <- function(object, ...) stats::step(object, ...)


##############################################
######## step.lmerModLmerTest()
##############################################
#' Backward Elimination for Linear Mixed Models
#'
#' Backward elimination of random-effect terms followed by backward elimination
#' of fixed-effect terms in linear mixed models.
#'
#' Tests of random-effects are performed using \code{\link{ranova}} (using
#' \code{reduce.terms = TRUE}) and tests of fixed-effects are performed using
#' \code{\link[=drop1.lmerModLmerTest]{drop1}}.
#'
#' The step method for \code{\link{lmer}} fits has a print method.
#'
#' @param object a fitted model object. For the \code{lmerModLmerTest} method
#' an \code{\link{lmer}} model fit (of class \code{"lmerModLmerTest"}.)
#' @param ddf the method for computing the denominator degrees of freedom and
#' F-statistics. \code{ddf="Satterthwaite"} (default) uses Satterthwaite's method;
#' \code{ddf="Kenward-Roger"} uses Kenward-Roger's method.
#' @param alpha.random alpha for random effects elimination
#' @param alpha.fixed alpha for fixed effects elimination
#' @param reduce.fixed reduce fixed effect structure? \code{TRUE} by default.
#' @param reduce.random reduce random effect structure? \code{TRUE} by default.
#' @param keep an optional character vector of fixed effect terms which should
#' not be considered for eliminated. Valid terms are given by
#' \code{attr(terms(object), "term.labels")}. Terms that are marginal to terms
#' in keep will also not be considered for eliminations.
#' @param ... currently not used.
#'
#' @return \code{step} returns a list with elements \code{"random"} and
#' \code{"fixed"} each
#' containing anova-like elimination tables. The \code{"fixed"} table is
#' based on \code{drop1} and the \code{"random"} table is
#' based on \code{ranova} (a \code{drop1}-like table for random effects). Both
#' tables have a column \code{"Eliminated"} indicating the order in which terms
#' are eliminated from the model with zero (\code{0}) indicating that the term
#' is not eliminated from the model.
#'
#' The \code{step} object also contains the final model as an attribute which
#' is extractable with \code{get_model(<step_object>)}.
#' @seealso \code{\link[=drop1.lmerModLmerTest]{drop1}} for tests of marginal
#' fixed-effect terms and \code{\link{ranova}} for a
#' \code{\link[=drop1.lmerModLmerTest]{drop1}}-like table of reduction of
#' random-effect terms.
#' @author Rune Haubo B. Christensen and Alexandra Kuznetsova
#' @export
#' @examples
#'
#' # Fit a model to the ham dataset:
#' fm <- lmer(Informed.liking ~ Product*Information+
#'              (1|Consumer) + (1|Product:Consumer)
#'            + (1|Information:Consumer), data=ham)
#'
#' # Backward elimination using terms with default alpha-levels:
#' (step_res <- step(fm))
#' final <- get_model(step_res)
#' anova(final)
#'
#' \dontrun{
#' # Fit 'big' model:
#' fm <- lmer(Informed.liking ~ Product*Information*Gender*Age +
#'              + (1|Consumer) + (1|Consumer:Product) +
#'              (1|Consumer:Information), data=ham)
#' step_fm <- step(fm)
#' step_fm # Display elimination results
#' final_fm <- get_model(step_fm)
#' }
#'
step.lmerModLmerTest <- function(object, ddf=c("Satterthwaite", "Kenward-Roger"),
                                 alpha.random=0.1, alpha.fixed=0.05,
                                 reduce.fixed=TRUE, reduce.random=TRUE,
                                 keep, ...) {
  # Check for and warn about deprecated arguments:
  ignored <- c("type", "fixed.calc", "lsmeans.calc", "difflsmeans.calc",
               "test.effs")
  dots <- list(...)
  for(nm in ignored) if(any(pmatch(names(dots), nm, nomatch = 0)))
    warning(paste0("Argument '", nm, "' is deprecated and ignored."))
  if(any(pmatch(names(dots), "keep.effs", nomatch = 0)))
    warning("Argument 'keep.effs' is deprecated: use 'keep' instead")

  # reduce random and fixed parts?
  if(!reduce.random) alpha.random <- 1
  if(!reduce.fixed) alpha.fixed <- 1
  if(missing(keep)) keep <- character(0L)
  # Reduce random and fixed parts:
  red_random <- eval.parent(reduce_random(object, alpha=alpha.random))
  model <- attr(red_random, "model")
  # 'model' may be 'lmerMod' rather than 'lmerModLmerTest', so we coerce to
  # 'lmerModLmerTest' if required:
  if(!inherits(model, "lmerModLmerTest"))
    model <- as_lmerModLmerTest(model)
  red_fixed <- eval.parent(reduce_fixed(model, ddf=ddf,
                                        alpha=alpha.fixed, keep=keep))
  # get 'reduction' tables:
  step_random <- ran_redTable(red_random)
  step_fixed <- fix_redTable(red_fixed)
  # organize results and return:
  step_list <- list(random=step_random, fixed=step_fixed)
  class(step_list) <- "step_list"
  attr(step_list, "model") <- attr(red_fixed, "model")
  attr(step_list, "drop1") <- attr(red_fixed, "drop1")
  step_list
}


##############################################
######## get_model()
##############################################
#' Extract Model from an Object
#'
#' @param x an object.
#' @param ... currently not used.
#'
#' @seealso \code{\link{get_model.step_list}}
#' @export
#' @keywords internal
get_model <- function(x, ...) UseMethod("get_model")


##############################################
######## get_model.step_list()
##############################################
#' @rdname step.lmerModLmerTest
#' @param x a step object.
#' @export
get_model.step_list <- function(x, ...) {
  attr(x, "model")
}


##############################################
######## print.step_list()
##############################################
#' @importFrom stats formula
#' @export
#' @keywords internal
print.step_list <- function(x, digits = max(getOption("digits") - 2L, 3L),
                            signif.stars = getOption("show.signif.stars"),
                            ...) {
  print(x[["random"]])
  cat("\n")
  print(x[["fixed"]])
  cat("\nModel found:", deparse(formula(attr(x, "model"))), sep="\n")
  invisible(x)
}


##############################################
######## plot.step_list()
##############################################
#' Plot LS-means for Backward Reduced Model
#'
#' Computes the LS-means for the final backward reduced model and passes these
#' to \code{\link{plot.ls_means}}.
#'
#' Error bars are confidence intervals - the default is 95% CI but the confidence
#' level can be changed.
#'
#' @param x a \code{step_list} object; the result of running
#' \code{\link[=step.lmerModLmerTest]{step}}.
#' @param y not used and ignored with a warning.
#' @param which optional character vector naming factors for which LS-means should
#' be plotted. If \code{NULL} (default) plots for all LS-means are generated.
#' @param mult if \code{TRUE} and there is more than one term for which to plot
#' LS-means the plots are organized in panels with \code{facet_wrap}.
#' @param pairwise pairwise differences of LS-means?
#' @param level confidence level.
#' @param ddf denominator degree of freedom method.
#' @param ... currently not used.
#'
#' @export
#' @author Rune Haubo B. Christensen and Alexandra Kuznetsova
#' @seealso \code{\link[=ls_means.lmerModLmerTest]{ls_means}} and
#' \code{\link{plot.ls_means}}
#' @keywords internal
#' @examples
#'
#' \dontrun{
#' # Fit example model:
#' tv <- lmer(Sharpnessofmovement ~ TVset * Picture +
#'              (1 | Assessor:TVset) + (1 | Assessor:Picture) +
#'              (1 | Assessor:Picture:TVset) + (1 | Repeat) + (1 | Repeat:Picture) +
#'              (1 | Repeat:TVset) + (1 | Repeat:TVset:Picture) + (1 | Assessor),
#'            data = TVbo)
#'
#' # Backward reduce the model:
#' (st <- step(tv)) # takes ~10 sec to run
#'
#' # Pairwise comparisons of LS-means for Picture and TVset:
#'   plot(st, which=c("Picture", "TVset"), pairwise = TRUE)
#' }
#'
plot.step_list <- function(x, y=NULL, which=NULL, pairwise=FALSE, mult=TRUE,
                           level=0.95, ddf=c("Satterthwaite", "Kenward-Roger"),
                           ...) {
  plot(ls_means(get_model(x), pairwise=pairwise, level=level, ddf=ddf),
       y=y, which=which, mult=mult)
}


##############################################
######## step utility functions below
##############################################

ran_redTable <- function(table) {
  aov <- attr(table, "ranova")[-1, , drop=FALSE]
  stopifnot(nrow(table) >= 1)
  tab <- rbind(cbind("Eliminated"=c(NA_real_, seq_len(nrow(table)-1)), table),
               cbind("Eliminated"=rep(0, nrow(aov)), aov))
  class(tab) <- c("anova", "data.frame")
  attr(tab, "heading") <- "Backward reduced random-effect table:\n"
  tab
}

fix_redTable <- function(table) {
  aov <- attr(table, "drop1")
  tab <- rbind(cbind("Eliminated"=seq_len(nrow(table)), table),
               cbind("Eliminated"=rep(0, nrow(aov)), aov))
  class(tab) <- c("anova", "data.frame")
  attr(tab, "heading") <- "Backward reduced fixed-effect table:"
  if(!is.null(ddf <- attr(table, "ddf"))) {
    ddf <- switch(ddf, "Satterthwaite" = "Satterthwaite",
                  "Kenward-Roger" = "Kenward-Roger")
    attr(tab, "heading") <-
      c(attr(tab, "heading"), paste("Degrees of freedom method:", ddf, "\n"))
  }
  tab
}


#' @importFrom stats formula update
#' @importFrom lme4 getME
reduce_random <- function(model, alpha=0.1) {
  ran <- ranova(model)
  reduced <- ran[1L, ]
  newfit <- model
  newform <- formula(model)
  forms <- attr(ran, "formulae")
  pvals <- ran[-1, "Pr(>Chisq)"]
  above <- (!is.na(pvals) & pvals > alpha)
  while(any(above)) {
    remove <- which.max(pvals)
    newform <- forms[[remove]]
    reduced <- rbind(reduced, ran[1 + remove, ])
    if(!has_ranef(newform)) { # If no random effects: fit with lm
      reml <- getME(newfit, "is_REML")
      lm_call <- get_lm_call(newfit, formula=newform)
      newfit <- eval.parent(as.call(lm_call))
      ran <- ranova_lm(newfit, REML=reml)
      break
    }
    newfit <- eval.parent(update(newfit, formula. = newform))
    # newfit <- update(newfit, formula = newform)
    ran <- ranova(newfit)
    forms <- attr(ran, "formulae")
    pvals <- ran[-1, "Pr(>Chisq)"]
    above <- (!is.na(pvals) & pvals > alpha)
  }
  attr(reduced, "model") <- newfit
  attr(reduced, "formula") <- newform
  attr(reduced, "ranova") <- ran
  reduced
}

ranova_lm <- function(model, REML=TRUE) {
  # Compute a ranova table for an lm-object only containing a '<none>' row
  # and the right header.
  aov <- mk_LRtab(get_logLik(model, REML=REML))
  rownames(aov) <- "<none>"
  head <- c("ANOVA-like table for random-effects: Single term deletions",
            "\nModel:", deparse(formula(model)))
  # attr(aov, "formulae") <- new_forms
  structure(aov, heading = head, class = c("anova", "data.frame"))
}

#' @importFrom stats nobs formula
reduce_fixed <- function(model, ddf=c("Satterthwaite", "Kenward-Roger"), alpha=0.05,
                         keep) {
  if(missing(keep)) keep <- character(0L)
  stopifnot(is.character(keep))
  term_names <- attr(terms(model), "term.labels")
  # Test validity of
  if(!all(keep %in% term_names)) {
    offending <- paste(setdiff(keep, term_names), collapse = " ")
    txt1 <- sprintf("Invalid 'keep' ignored: %s.", offending)
    txt2 <- sprintf("Valid terms are: %s.", paste(term_names, collapse = " "))
    warning(paste(txt1, txt2, sep="\n"), call. = FALSE)
  }
  ddf <- match.arg(ddf)
  aov <- if(inherits(model, "lmerMod")) drop1.lmerModLmerTest(model, ddf=ddf) else
    drop1(model, test="F")[-1L, , drop=FALSE]
  reduced <- aov[0L, ]
  newfit <- model
  newform <- orig_form <- formula(model)
  nobs_model <- nobs(model)
  terms <- rownames(aov)
  consider <- setdiff(terms, keep)
  pvals <- aov[consider, "Pr(>F)"]
  above <- (!is.na(pvals) & pvals > alpha)
  if(any(above)) while(any(above)) {
    remove <- consider[which.max(pvals)]
    newform <- rm_complete_terms(remove, orig_form, random = FALSE)[[1L]]
    reduced <- rbind(reduced, aov[remove, ])
    newfit <- eval.parent(update(newfit, formula = newform))
    # newfit <- update(newfit, formula = newform)
    nobs_newfit <- nobs(newfit)
    if(all(is.finite(c(nobs_model, nobs_newfit))) && nobs_newfit != nobs_model)
      stop("number of rows in use has changed: remove missing values?",
           call.=FALSE)
    aov <- if(inherits(newfit, "lmerMod")) drop1.lmerModLmerTest(newfit, ddf=ddf) else
      drop1(newfit, test="F")[-1L, , drop=FALSE]
    # aov <- drop1(newfit)
    orig_form <- formula(newfit)
    terms <- rownames(aov)
    consider <- setdiff(terms, keep)
    pvals <- aov[consider, "Pr(>F)"]
    above <- (!is.na(pvals) & pvals > alpha)
  }
  attr(reduced, "model") <- newfit
  attr(reduced, "formula") <- newform
  attr(reduced, "drop1") <- aov
  attr(reduced, "ddf") <- if(inherits(model, "lmerMod")) ddf else NULL
  reduced
}
