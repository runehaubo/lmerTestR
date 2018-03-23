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
# drop1.R - drop1 method for lmerModLmerTest objects

# ------- Contents: --------
#
# drop1.lmerModLmerTest
#
# --- Utility functions: ---
#
# get_Ldiffmat
# get_Ldiffmat2
#

##############################################
######## drop1.lmerModLmerTest
##############################################
#' Drop Marginal Terms from Model
#'
#' Computes the F-test for all marginal terms, i.e. terms that can be dropped
#' from the model while respecting the hierarchy of terms in the model.
#'
#' Simple marginal contrasts are used for all marginal terms unless the design
#' matrix is rank deficient. In that case (and if \code{force_get_contrasts} is
#' \code{TRUE}) the contrasts (i.e. restriction matrices on the design matrix
#' of the full model) are computed by comparison of the design matrices
#' for full and restricted models. The set of marginal terms considered for
#' dropping are computed using \code{drop.scope(terms(object))}.
#'
#' Since all tests are based on tests of contrasts in the full model, no
#' models are being (re)fitted.
#'
#' @param object an \code{\link{lmer}} model fit (of class
#' \code{"lmerModLmerTest"}.)
#' @param scope optional character vector naming terms to be dropped from the
#' model. Note that only marginal terms can be dropped. To see which terms are
#' marginal, use \code{drop.scope(terms(object))}.
#' @param ddf the method for computing the denominator degrees of freedom and
#' F-statistics. \code{ddf="Satterthwaite"} (default) uses Satterthwaite's method;
#' \code{ddf="Kenward-Roger"} uses Kenward-Roger's method.
#' \code{ddf = "lme4"} returns the \code{drop1} table for \code{merMod} objects
#' as defined in package \pkg{lme4}.
#' @param force_get_contrasts enforce computation of contrast matrices by a
#' method in which the design matrices for full and restricted models are
#' compared.
#' @param ... currently not used.
#'
#' @author Rune Haubo B. Christensen
#' @seealso \code{\link{ranova}} for tests of marginal random terms.
#' @return An anova-like table with F-tests of marginal terms.
#' @export
#'
#' @importFrom stats drop1 drop.scope terms formula
#' @examples
#'
#' # Basic usage:
#' fm <- lmer(angle ~ recipe + temp + (1|recipe:replicate), cake)
#' drop1(fm)
#' drop1(fm, ddf="Kenward-Roger") # Alternative DenDF and F-test method
#' drop1(fm, ddf="lme4", test="Chi") # Asymptotic Likelihood ratio tests
#'
#' # Consider a rank-deficient design matrix:
#' fm <- lmer(angle ~ recipe + temp + temperature + (1|recipe:replicate), cake)
#' # Here temp accounts for the linear effect of temperature, and
#' # temperature is an (ordered) factor that accounts for the remaining
#' # variation between temperatures (4 df).
#' drop1(fm)
#' # While temperature is in the model, we cannot test the effect of dropping
#' # temp. After removing temperature we can test the effect of dropping temp:
#' drop1(lmer(angle ~ recipe + temp + (1|recipe:replicate), cake))
#'
#' # Polynomials:
#' # Note that linear terms should usually not be dropped before squared terms.
#' # Therefore 'Days' should not be dropped before 'I(Days^2)' despite it being
#' # tested here:
#' fm <- lmer(Reaction ~ Days + I(Days^2) + (Days|Subject), sleepstudy)
#' drop1(fm)
#' # Using poly() provides a test of the whole polynomial structure - not a
#' # separate test for the highest order (squared) term:
#' fm <- lmer(Reaction ~ poly(Days, 2) + (Days|Subject), sleepstudy)
#' drop1(fm)
#'
drop1.lmerModLmerTest <- function(object, scope, ddf=c("Satterthwaite", "Kenward-Roger", "lme4"),
                                  force_get_contrasts=FALSE, ...) {
  ddf <- match.arg(ddf)
  if(ddf == "lme4") return(NextMethod())
  marg_terms <- drop.scope(terms(object))
  if(missing(scope)) scope <- marg_terms else {
    if(length(scope) == 0 || !is.character(scope))
      stop("'scope' should be a character vector naming terms to be dropped")
    if(!all(scope %in% marg_terms))
      stop("Only marginal terms can be dropped from the model")
  }
  # Get contrasts for marginal terms:
  X <- model.matrix(object)
  Llist <- get_contrasts_marginal(object)
  if(length(scope)) {
    Llist <- Llist[scope] # retain contrasts for terms in scope
    if(!is.null(attr(X, "col.dropped")) || force_get_contrasts) {
      # Compute L directly if model is rank deficient or force_get_contrasts is TRUE:
      orig_form <- formula(object)
      new_forms <- lapply(rm_complete_terms(scope, orig_form, random=FALSE), nobars)
      # Compute list of contrast matrices as 'diffs' to orig. X:
      Llist <- if(!length(new_forms)) list() else
        lapply(new_forms, function(form) {
          suppressWarnings(x <- model.matrix(form[-2], data=model.frame(object),
                                             contrasts.arg = attr(X, "contrasts")))
          L <- get_Ldiffmat2(x, X) # L may be length 0 if x == X (rank-deficint fits.)
          if(!length(L)) rep(NA_real_, ncol(X)) else L
        })
    }
  }
  # Compute anova-like table:
  aov <- rbindall(lapply(Llist, function(L) contestMD(object, L, ddf = ddf)))
  # Format results:
  method <- switch(ddf, "Satterthwaite" = "Satterthwaite's",
                   "Kenward-Roger" = "Kenward-Roger's")
  attr(aov, "heading") <-
    c(paste("Single term deletions using", method, "method:"),
      "\nModel:", deparse(formula(object)))
  attr(aov, "hypotheses") <- Llist
  attr(aov, "ddf") <- ddf
  class(aov) <- c("anova", "data.frame")
  aov
}


get_Ldiffmat <- function(A0, A) {
  Rank <- function(X) qr(X)$rank
  Q <- qr.Q(qr(cbind(A0, A)))
  rA0 <- Rank(A0)
  rA <- Rank(A)
  set <- if(rA0 < rA) (rA0+1):rA else numeric(0L)
  Q2 <- Q[, set, drop=FALSE]
  L <- t(Q2) %*% A
  L <- t(qr.Q(qr(t(L)))) # Orthonormalize contrast
  L
}

#' @importFrom stats .lm.fit resid
get_Ldiffmat2 <- function(X0, X) {
  # X : design matrix for the full model
  # X0: design matrix for the restricted model
  # R is the residual of the orthogonal projection of X on X0, thus
  # R is orthogonal to X0 and a subspace of X, and
  # Lt is a restriction matrix on X.
  R <- resid(.lm.fit(x=X0, y=X))
  R <- R[, colSums(abs(R)) > 1e-8]
  Lt <- crossprod(X, R)
  Lt[] <- zapsmall(qr.Q(qr(Lt))) # orthonormalize contrasts
  t(Lt)
}
