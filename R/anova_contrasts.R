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
# anova_contrasts.R - functions of the form get_contrasts_xxxx() used by anova
#  to get contrasts for model terms.

# Functions in this file:

# Standard contrast functions:
# get_contrast_type3              # type = 3
# get_contrast_type2_unfolded     # type = 2
# get_contrast_type1              # type = 1
# get_contrast_marginal           # type = marginal
#
# get_contrast_yates              # type = yates
# get_contrast_type2              # type = 2b


##############################################
######## get_contrasts_type3
##############################################
#' Contrasts for Type III Tests
#'
#' @param model model object.
#' @param which optional character vector naming terms for which to compute the
#' the contrasts.
#'
#' @return list of contrast matrices.
#' @importFrom stats terms
#' @keywords internal
get_contrasts_type3 <- function(model, which=NULL) {
  term_names <- attr(terms(model), "term.labels")
  # Extract original design matrix:
  Xorig <- model.matrix(model)
  # Assumes Xorig is full (column) rank
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(Xorig) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))

  # Extract contrast coding in Xorig:
  codings <- unlist(attr(Xorig, "contrast"))
  # If only treatment contrasts are used we can just return the type 3
  # contrasts for contr.treatment coding:
  if(length(codings) > 0 &&
     all(is.character(codings)) && all(codings %in% c("contr.treatment")))
    return(extract_contrasts_type3(model, X=Xorig))
  # otherwise we need to map the type III contrasts to whatever contrast
  # coding was used:
  X <- get_model_matrix(model, type="remake", contrasts="contr.treatment")
  # Ensure that X is full (column) rank:
  X <- ensure_full_rank(X, silent=TRUE, test.ans=FALSE)
  # Extract contrasts assuming contr.treatment coding:
  type3ctr <- extract_contrasts_type3(model, X=X)
  map <- zapsmall(ginv(X) %*% Xorig) # Maps between contrast codings
  rownames(map) <- colnames(X)
  lapply(type3ctr[which], function(L) L %*% map)
}


##############################################
######## get_contrasts_type2_unfolded
##############################################
#' @importFrom stats model.matrix terms
get_contrasts_type2_unfolded <- function(model, which=NULL) {
  # Computes the 'genuine type II contrast' for all terms that are
  # contained in other terms. For all terms which are not contained in other
  # terms, the simple marginal contrast is computed.
  X <- model.matrix(model)
  Terms <- terms(model)
  term_names <- attr(Terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))

  is_contained <- containment(model)
  do_marginal <- names(is_contained)[sapply(is_contained, length) == 0L]
  do_type2 <- setdiff(term_names, do_marginal)

  if(!length(do_marginal)) list() else
    Llist <- get_contrasts_marginal(model, which=do_marginal)
  if(length(do_type2))
    Llist <- c(Llist, get_contrasts_type2(model, which=do_type2))
  Llist[term_names]
}


##############################################
######## get_contrasts_type1
##############################################
#' Type I ANOVA table contrasts
#'
#' @param model a model object with \code{terms} and \code{model.matrix} methods.
#'
#' @return List of contrast matrices - one contrast matrix for each model term.
#' @importFrom stats setNames
#' @author Rune Haubo B. Christensen
#'
#' @keywords internal
get_contrasts_type1 <- function(model) {
  terms <- terms(model)
  X <- model.matrix(model)
  p <- ncol(X)
  if(p == 0L) return(list(matrix(numeric(0L), nrow=0L))) # no fixef
  if(p == 1L && attr(terms, "intercept")) # intercept-only model
    return(list(matrix(numeric(0L), ncol=1L)))
  # Compute 'normalized' doolittle factorization of XtX:
  L <- if(p == 1L) matrix(1L) else t(doolittle(crossprod(X))$L)
  dimnames(L) <- list(colnames(X), colnames(X))
  # Determine which rows of L belong to which term:
  ind.list <- term2colX(terms, X)[attr(terms, "term.labels")]
  lapply(ind.list, function(rows) L[rows, , drop=FALSE])
}


##############################################
######## get_contrasts_marginal
##############################################
#' @importFrom stats model.matrix terms
get_contrasts_marginal <- function(model, which=NULL) {
  # Computes marginal contrasts.
  #
  # No tests of conformity with coefficients are implemented
  #
  # returns a list
  X <- model.matrix(model)
  terms <- terms(model)
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))
  ## FIXME: test use of 'which' arg.

  # Compute map from terms to columns in X and contrasts matrix
  term2colX <- term2colX(terms, X)
  L <- structure(diag(ncol(X)), dimnames = list(colnames(X), colnames(X)))

  # Extract contrast for each term - return as named list:
  which <- setNames(as.list(which), which)
  lapply(which, function(term) {
    L[term2colX[[term]], , drop=FALSE]
  })
}


##############################################
######## get_contrasts_yates
##############################################
get_contrasts_yates <- function(model) {
  # Is this really type 4?
  X <- model.matrix(model)
  Terms <- terms(model)
  term_names <- attr(Terms, "term.labels")

  is_contained <- containment(model)
  do_marginal <- names(is_contained)[sapply(is_contained, length) == 0L]
  not_marginal <- setdiff(term_names, do_marginal)
  # Split not_marginal in do_yates and do_type2:
  do_yates <- need_yates(model)
  do_type2 <- setdiff(not_marginal, do_yates)

  if(!length(do_marginal)) list() else
    Llist <- get_contrasts_marginal(model, which=do_marginal)
  if(length(do_yates))
    Llist <- c(Llist, get_yates_contrast(model, which=do_yates))
  if(length(do_type2)) {
    data_classes <- attr(terms(model, fixed.only=FALSE), "dataClasses")
    Llist <- c(Llist, get_contrasts_type2(model, which=do_type2))
  }
  Llist[term_names]
}


##############################################
######## get_contrasts_type2
##############################################
get_contrasts_type2 <- function(model, which=NULL) {
  # Computes the type 2 contrasts - either for all terms or for those
  # included in 'which' (a chr vector naming model terms).
  # returns a list
  X <- model.matrix(model)
  terms <- terms(model)
  data_classes <- attr(terms(model, fixed.only=FALSE), "dataClasses")
  if(is.null(asgn <- attr(X, "assign")))
    stop("design matrix 'X' should have a non-null 'assign' attribute")
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)

  # Compute containment:
  is_contained <- containment(model)

  # Compute term asignment list: map from terms to columns in X
  has_intercept <- attr(terms, "intercept") > 0
  col_terms <- if(has_intercept) c("(Intercept)", term_names)[asgn + 1] else
    term_names[asgn[asgn > 0]]
  if(!length(col_terms) == ncol(X)) # should never happen.
    stop("An error happended when computing Type II contrasts")
  term2colX <- split(seq_along(col_terms), col_terms)[unique(col_terms)]

  # Compute contrast for each term - return as named list:
  lapply(which, function(term) {
    # Reorder the cols in X to [, unrelated_to_term, term, contained_in_term]
    cols_term <- unlist(term2colX[c(term, is_contained[[term]])])
    Xnew <- cbind(X[, -cols_term, drop=FALSE], X[, cols_term, drop=FALSE])
    # Compute order of terms in Xnew:
    newXcol_terms <- c(col_terms[-cols_term], col_terms[cols_term])
    # Compute Type I contrasts for the reordered X:
    Lc <- t(doolittle(crossprod(Xnew))$L)
    dimnames(Lc) <- list(colnames(Xnew), colnames(Xnew))
    # Extract rows for term and get original order of columns:
    Lc[newXcol_terms == term, colnames(X), drop=FALSE]
  })
}
