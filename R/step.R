
#' @rdname step.lmerModLmerTest
#' @export
step <- function(object, ...) UseMethod("step")

#' @rdname step.lmerModLmerTest
#' @export
step.default <- function(object, ...) stats::step(object, ...)

#' Backward Elimination for Linear Mixed Models
#'
#' @param object the model
#' @param ddf degrees of freedom
#' @param alpha.random alpha for random effects
#' @param alpha.fixed alpha for fixed effects
#' @param ... currently not used
#'
#' @return a list ...
#' @export
#'
step.lmerModLmerTest <- function(object, ddf=c("Satterthwaite", "KR"),
                                 alpha.random=0.1, alpha.fixed=0.05, ...) {
  red_random <- eval.parent(reduce_random(object, alpha=alpha.random))
  red_fixed <- eval.parent(reduce_fixed(attr(red_random, "model"), ddf=ddf,
                                        alpha=alpha.fixed))
  step_random <- ran_redTable(red_random)
  step_fixed <- fix_redTable(red_fixed)

  step_list <- list(random=step_random, fixed=step_fixed)
  class(step_list) <- "step_list"
  attr(step_list, "model") <- attr(red_fixed, "model")
  attr(step_list, "drop1") <- attr(red_fixed, "drop1")
  step_list
}

#' @rdname step.lmerModLmerTest
#' @export
get_model <- function(x, ...) UseMethod("get_model")

#' @rdname step.lmerModLmerTest
#' @param x a step object
#' @export
get_model.step_list <- function(x, ...) {
  attr(x, "model")
}


#' @importFrom stats formula
print.step_list <- function(x, digits = max(getOption("digits") - 2L, 3L),
                            signif.stars = getOption("show.signif.stars"),
                            ...) {
  print(x[["random"]])
  cat("\n")
  print(x[["fixed"]])
  cat("\nFound model:", deparse(formula(attr(x, "model"))), sep="\n")
  invisible(x)
}


ran_redTable <- function(table) {
  aov <- attr(table, "ranova")[-1, , drop=FALSE]
  stopifnot(nrow(table) >= 1)
  tab <- rbind(cbind("Step"=c(NA_real_, seq_len(nrow(table)-1)), table),
               cbind("Step"=rep(0, nrow(aov)), aov))
  class(tab) <- c("anova", "data.frame")
  attr(tab, "heading") <- "Backward reduced random-effect table:\n"
  tab
}


fix_redTable <- function(table) {
  aov <- attr(table, "drop1")
  tab <- rbind(cbind("Step"=seq_len(nrow(table)), table),
               cbind("Step"=rep(0, nrow(aov)), aov))
  class(tab) <- c("anova", "data.frame")
  attr(tab, "heading") <- "Backward reduced fixed-effect table:"
  if(!is.null(ddf <- attr(table, "ddf"))) {
    ddf <- switch(ddf, "Satterthwaite" = "Satterthwaite",
                  "KR" = "Kenward-Roger")
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
    newfit <- eval.parent(update(newfit, formula = newform))
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
  aov <- mk_LRtab(get_logLik(model, REML=REML))
  rownames(aov) <- "<none>"
  head <- c("ANOVA-like table for random-effects: Single term deletions",
            "\nModel:", deparse(formula(model)))
  # attr(aov, "formulae") <- new_forms
  structure(aov, heading = head, class = c("anova", "data.frame"))
}

#' @importFrom stats nobs formula
reduce_fixed <- function(model, ddf=c("Satterthwaite", "KR"), alpha=0.05) {
  ddf <- match.arg(ddf)
  aov <- if(inherits(model, "lmerMod")) drop1(model, ddf=ddf) else
    drop1(model, test="F")[-1L, , drop=FALSE]
  reduced <- aov[0L, ]
  newfit <- model
  newform <- orig_form <- formula(model)
  nobs_model <- nobs(model)
  terms <- rownames(aov)
  pvals <- aov[, "Pr(>F)"]
  above <- (!is.na(pvals) & pvals > alpha)
  if(any(above)) while(any(above)) {
    remove <- terms[which.max(pvals)]
    newform <- rm_complete_terms(remove, orig_form, random = FALSE)[[1L]]
    reduced <- rbind(reduced, aov[remove, ])
    newfit <- eval.parent(update(newfit, formula = newform))
    # newfit <- update(newfit, formula = newform)
    nobs_newfit <- nobs(newfit)
    if(all(is.finite(c(nobs_model, nobs_newfit))) && nobs_newfit != nobs_model)
      stop("number of rows in use has changed: remove missing values?",
           call.=FALSE)
    aov <- if(inherits(newfit, "lmerMod")) drop1(newfit, ddf=ddf) else
      drop1(newfit, test="F")[-1L, , drop=FALSE]
    # aov <- drop1(newfit)
    orig_form <- formula(newfit)
    terms <- rownames(aov)
    pvals <- aov[, "Pr(>F)"]
    above <- (!is.na(pvals) & pvals > alpha)
  }
  attr(reduced, "model") <- newfit
  attr(reduced, "formula") <- newform
  attr(reduced, "drop1") <- aov
  attr(reduced, "ddf") <- if(inherits(model, "lmerMod")) ddf else NULL
  reduced
}
