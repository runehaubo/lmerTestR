# anova_contrasts.R - functions of the form get_contrasts_xxxx() used by anova
#  to get contrasts for model terms.

# Functions in this file:

# Standard contrast functions:
# get_contrast_type3              # type = 3
# get_contrast_type2_unfolded     # type = 2
# get_contrast_type1              # type = 1
# get_contrast_marginal           # type = marginal
#
# get_contrast_type3old           # type = 3c
# get_contrast_type3b             # type = 3b / yates
# get_contrast_type2              # type = 2b
# get_contrast_yates # == type 4?

##############################################
######## get_contrasts_type3
##############################################
#' @importFrom MASS ginv
#' @importFrom stats terms resid lm.fit
get_contrasts_type3 <- function(model, which=NULL) {
  Terms <- terms(model)
  X <- model.matrix(model)
  term_names <- attr(Terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))

  # Check that appropriate contrast coding is used:
  codings <- unlist(attr(X, "contrast"))
  if(length(codings) > 0 &&
     !all(is.character(codings) && codings %in% c("contr.treatment", "contr.SAS")))
    warning("Type 3c only applies if all factors are coded with 'contr.treatment' ",
            "or 'contr.SAS'")
  # Get 'complete' design matrix:
  rdX <- get_rdX(model, do.warn = TRUE, use_contrasts = TRUE)
  # cols for aliased coefs should be removed in X; not in rdX.
  # This makes ginv(X) unique!
  L <- zapsmall(t(MASS::ginv(X) %*% rdX))
  dimnames(L) <- list(names(attr(rdX, "param")), colnames(X))

  # Orthogonalize contrasts for terms which are contained in other terms:
  map <- term2colX(Terms, X)
  is_contained <- containment(model)
  # Orthogonalize higher order terms before lower order terms:
  terms_order <- attr(Terms, "order")
  orthog_order <- term_names[order(terms_order, decreasing = TRUE)]
  for(term in orthog_order) {
    if(length(contains <- is_contained[[term]]) > 0) {
      L[, map[[term]]] <-
        zapsmall(resid(lm.fit(x=L[, map[[contains]], drop=FALSE],
                              y=L[, map[[term]], drop=FALSE])))
    }
  }
  # # Get indicator for which rows in L to keep:
  # keep <- if(all(colnames(X) %in% rownames(L))) colnames(X) else {
  #   param_ind <- term2colX(Terms, rdX)
  #   which <- if(attr(Terms, "intercept") == 1) c("(Intercept)", term_names) else
  #     term_names
  #   unlist(lapply(which, function(term) {
  #     param_ind[[term]][seq_along(map[[term]])]
  #   }))
  # }
  # Return list of contrast matrices - one for each term:
  lapply(map[which], function(term) t(L[colnames(X), term, drop=FALSE]))
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
#' @param keep_intercept defaults to \code{FALSE}. If \code{TRUE} a contrast
#' for the intercept is included in the returned list of contrast matrices.
#'
#' @return List of contrast matrices - one contrast matrix for each model term.
#' @importFrom stats setNames
#' @author Rune Haubo B. Christensen
#'
#' @keywords internal
get_contrasts_type1 <- function(model, keep_intercept = FALSE) {
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
  asgn <- attr(X, "assign")
  stopifnot(!is.null(asgn))
  term_labels <- attr(terms, "term.labels")
  term_names <- c("(Intercept)", term_labels)
  term_names <- term_names[1 + unique(asgn)] # order appropriately
  # Compute list of row indicators for L matrix:
  ind.list <- setNames(split(1L:p, asgn), nm=term_names)
  ind.list <- ind.list[term_labels] # rm intercept if present
  lapply(ind.list, function(rows) L[rows, , drop=FALSE])
  # FIXME: implement keep_intercept = TRUE. Do we need it, want it?
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
######## get_contrasts_type3old
##############################################
#' Compute Type III Contrast Matrices
#'
#' Experimental - not extensively tested.
#'
#' @param model a model object; lmerMod or lmerModLmerTest
#'
#' @return a list of contrast matrices for each model term.
#' @author Rune Haubo B. Christensen
#'
#' @keywords internal
get_contrasts_type3old <- function(model) {
  # Compute list of contrast matrices for type 3 tests.
  #
  # Get original design matrix:
  X <- model.matrix(model)
  # Catch boundary cases:
  p <- ncol(X)
  if(p == 0L) return(list(matrix(numeric(0L), nrow=0L))) # no fixef
  if(p == 1L && attr(terms(model), "intercept")) # intercept-only model
    return(list(matrix(numeric(0L), ncol=1L)))
  # Check that valid contrasts are used for factors:
  if(length(contr <- attr(X, "contrasts")) > 0) {
    is_valid_contrast <-
      sapply(contr, '%in%', c("contr.treatment", "contr.SAS"))
    if(!all(is_valid_contrast))
      stop("invalid contrasts found for type 3: Use 'contr.treatment' or 'contr.SAS' for all factors")
  }
  # get rank-deficient design-matrix X:
  rdX <- get_rdX(model)
  is_coef <- attr(rdX, "is_coef")
  # Compute terms to be tested (term_names):
  Terms <- terms(model, fixed.only=TRUE)
  term_names <- attr(Terms, "term.labels")
  ord <- unique(attr(X, "assign"))
  term_names <- term_names[ord] # order terms (when are they not?)
  # Compute 'general set' L:
  L <- general_L(rdX, is_coef)
  # Compute Type 3 Lc for each term:
  data_classes <- attr(terms(model, fixed.only=FALSE), "dataClasses")
  term_names <- setNames(as.list(term_names), term_names)
  Lc_list <- lapply(term_names, function(nm) {
    Lc <- contrast_type3SAS(nm, Terms, data_classes, L)
    if(!is.matrix(Lc)) Lc <- matrix(Lc, ncol=length(Lc),
                                    dimnames = list(NULL, names(Lc)))
    Lc <- Lc[, is_coef, drop=FALSE]
    colnames(Lc) <- colnames(X)
    Lc
  })
}


##############################################
######## get_contrasts_type3b
##############################################
get_contrasts_type3b <- function(model) {
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
  factor_mat <- attr(terms, "factors")
  # list: for each term 'T' a vector of terms in which 'T' is contained:
  contained_in <- setNames(lapply(term_names, function(term) {
    term_names[relatives(data_classes, term, term_names, factor_mat)]
  }), term_names)

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
    cols_term <- unlist(term2colX[c(term, contained_in[[term]])])
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
