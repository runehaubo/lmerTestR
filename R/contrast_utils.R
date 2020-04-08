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
# contrast-utils.R - utility functions for contrasts, terms and anova

# -------- Contents: --------
#
# containment
# term_contain
# relatives
# doolittle
# ensure_full_rank
# get_rdX
# extract_contrasts_type3
# get_yates_contrast


##############################################
######## containment()
##############################################
#' Determine the Containment Structure for All Terms in a Model
#'
#' See \code{\link{term_contain}} for details about containment.
#'
#' @param object a model object, e.g. of class \code{lm} or \code{merMod}.
#'
#' @return a list with one element for each term in the model. Each element/term
#' is a character vector of terms that the term is contained in.
#' @importFrom stats terms setNames
#' @keywords internal
containment <- function(object) { # lm or merMod
  # For all terms 'T' in object compute the terms
  # Return a list:
  #   for each term 'T' a vector of terms that contain 'T'.
  terms <- terms(object)
  data_classes <- attr(terms(object, fixed.only=FALSE), "dataClasses")
  # Note: need fixed.only for merMod objects to get dataClasses
  term_names <- attr(terms, "term.labels")
  factor_mat <- attr(terms, "factors")
  lapply(setNames(term_names, term_names), function(term) {
    term_names[term_contain(term, factor_mat, data_classes, term_names)]
  })
}

##############################################
######## term_contain()
##############################################
#' Determine which Terms Contain a Term
#'
#' The definition of \emph{containment} follows from the SAS documentation on
#' "The Four Types of Estimable Functions".
#'
#' Containment is defined for two model terms, say, F1 and F2 as:
#' F1 is contained in F2 (F2 contains F1) if
#' \enumerate{
#' \item F1 and F2 involve the same continuous variables (if any)
#' \item F2 involve more factors than F1
#' \item All factors in F1 (if any) are part of F2
#' }
#' The intercept, though not really a model term, is defined by SAS to be
#' contained in all factor terms, but it is not contained in any
#' effect involving a continuous variable.
#'
#' @param term character; name of a model term and one of \code{term_names}.
#' @param factors the result of \code{attr(terms_object, "factors")}.
#' @param dataClasses the result of
#' \code{attr(terms(model, fixed.only=FALSE), "dataClasses")}. Note that
#' \code{fixed.only=FALSE} is only needed for \code{merMod} objects, but does
#' no harm for \code{lm} objects.
#' @param term_names the result of \code{attr(terms_object, "term.labels")}.
#'
#' @return a logical vector indicating for each term in \code{term_names} if
#' it contains \code{term}.
#' @importFrom stats setNames
#' @keywords internal
term_contain <- function(term, factors, dataClasses, term_names) {
  get_vars <- function(term)
    # Extract vector of names of all variables in a term
    rownames(factors)[factors[, term] == 1]
  contain <- function(F1, F2) {
    # Returns TRUE if F1 is contained in F2 (i.e. if F2 contains F1)
    # F1, F2: Names of terms, i.e. attr(terms_object, "term.labels")
    all(vars[[F1]] %in% vars[[F2]]) && # all variables in F1 are also in F2
      length(setdiff(vars[[F2]], vars[[F1]])) > 0L && # F2 involve more variables than F1
      setequal(numerics[[F1]], numerics[[F2]]) # F1 and F2 involve the same covariates (if any)
  }
  # Get (named) list of all variables in terms:
  vars <- lapply(setNames(term_names, term_names), get_vars)
  # Get (named) list of all _numeric_ variables in all terms:
  numerics <- lapply(vars, function(varnms)
    varnms[which(dataClasses[varnms] == "numeric")])
  # Check if 'term' is contained in each model term:
  sapply(term_names, function(term_nm) contain(term, term_nm))
}


##############################################
######## relatives()
##############################################
# relatives <- function(classes.term, term, term_names, factors) {
#   ## checks if the terms have the same number of covariates (if any)
#   checkCovContain <- function(term1, term2) {
#     num.numeric <- which(classes.term=="numeric")
#     num.numeric.term1 <- which((num.numeric %in% which(factors[,term1]!=0))==TRUE)
#     num.numeric.term2 <- which((num.numeric %in% which(factors[,term2]!=0))==TRUE)
#     if((length(num.numeric.term1)>0 && length(num.numeric.term2)>0)||
#        (length(num.numeric.term1)==0 && length(num.numeric.term2)==0))
#       return(all(num.numeric.term2 == num.numeric.term1))
#     else
#       return(FALSE)
#   }
#   is.relative <- function(term1, term2) {
#     all(!(factors[, term1] & (!factors[, term2]))) && checkCovContain(term1, term2)
#   }
#   if(length(term_names) == 1) return(NULL)
#   which.term <- which(term == term_names)
#   (1:length(term_names))[-which.term][sapply(term_names[-which.term],
#                                              function(term2) is.relative(term, term2))]
# }


##############################################
######## doolittle()
##############################################
#' Doolittle Decomposition
#'
#' @param x a numeric square matrix with at least 2 columns/rows.
#' @param eps numerical tolerance on the whether to normalize with components
#' in \code{L} with the diagonal elements of \code{U}.
#'
#' @return a list with two matrices of the same dimension as \code{x}:
#' \item{L}{lower-left unit-triangular matrix}
#' \item{U}{upper-right triangular matrix (\emph{not} unit-triangular)}
#'
#' @keywords internal
doolittle <- function(x, eps = 1e-6) {
  if(!is.matrix(x) || ncol(x) != nrow(x) || !is.numeric(x))
    stop("argument 'x' should be a numeric square matrix")
  stopifnot(ncol(x) > 1L)
  n <- nrow(x)
  L <- U <- matrix(0, nrow=n, ncol=n)
  diag(L) <- rep(1, n)
  for(i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for(j in 1:n) {
      U[i,j] <- x[i,j]
      if (im1 > 0) {
        for(k in 1:im1) {
          U[i,j] <- U[i,j] - L[i,k] * U[k,j]
        }
      }
    }
    if ( ip1 <= n ) {
      for ( j in ip1:n ) {
        L[j,i] <- x[j,i]
        if ( im1 > 0 ) {
          for ( k in 1:im1 ) {
            L[j,i] <- L[j,i] - L[j,k] * U[k,i]
          }
        }
        L[j, i] <- if(abs(U[i, i]) < eps) 0 else L[j,i] / U[i,i]
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list( L=L, U=U )
}


##############################################
######## ensure_full_rank()
##############################################
#' Ensure a Design Matrix has Full (Column) Rank
#'
#' Determine and drop redundant columns using the \code{\link{qr}}
#' decomposition.
#'
#' @param X a design matrix as produced by \code{model.matrix}.
#' @param tol \code{qr} tolerance.
#' @param silent throw message if columns are dropped from \code{X}? Default
#' is \code{FALSE}.
#' @param test.ans Test if the resulting/returned matrix has full rank? Default
#' is \code{FALSE}.
#'
#' @return A design matrix in which redundant columns are dropped
#' @keywords internal
ensure_full_rank <- function(X, tol = 1e-7, silent = FALSE, test.ans = FALSE) {
  ### works if ncol(X) >= 0 and nrow(X) >= 0
  ## test and match arguments:
  stopifnot(is.matrix(X))
  silent <- as.logical(silent)[1]
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(X, tol = tol, LAPACK = FALSE)
  if(qr.X$rank == ncol(X)) {
    ## return X if X has full column rank
    return(X)
  }
  if(!silent) ## message the no. dropped columns:
    message(gettextf("Design is column rank deficient so dropping %d coef",
                     ncol(X) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of X:
  keep <- with(qr.X, pivot[seq_len(rank)])
  newX <- X[, keep, drop = FALSE]
  sel <- with(qr.X, pivot[-seq_len(rank)])
  ## Copy old attributes:
  if(!is.null(contr <- attr(X, "contrasts"))) attr(newX, "contrasts") <- contr
  if(!is.null(asgn <- attr(X, "assign")))  attr(newX, "assign") <- asgn[-sel]
  ## did we succeed? stop-if-not:
  if(test.ans && qr.X$rank != qr(newX)$rank)
    stop(gettextf("Determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX)
}


##############################################
######## get_rdX()
##############################################
#' Compute the 'Full' Rank-Deficient Design Matrix
#'
#'
#' @param model a model object; lmerMod or lmerModLmerTest.
#' @param do.warn throw a message if there is no data for some factor
#' combinations.
#'
#' @return the rank-deficien design matrix
#' @author Rune Haubo B. Christensen
#' @keywords internal
#'
#' @importFrom stats as.formula model.frame terms model.matrix
get_rdX <- function(model, do.warn=TRUE) {
  # Compute rank-deficient design-matrix X usign contr.treatment coding.
  #
  # model: terms(model), model.frame(model), fixef(model)
  Terms <- terms(model, fixed.only=TRUE)
  term_names <- attr(Terms, "term.labels")
  df <- model.frame(model)
  # Compute rank-deficient (full) design-matrix, X:
  rdXi <- if(length(term_names)) lapply(term_names, function(trm) {
    form <- as.formula(paste0("~ 0 + ", trm))
    model.matrix(form, data=df) # no contrast arg
  }) else list(model.matrix(~ 1, data=df)[, -1, drop=FALSE])
  rdX <- do.call(cbind, rdXi)
  param_names <- unlist(lapply(rdXi, colnames))
  # Potentially add intercept:
  has_intercept <- attr(Terms, "intercept") != 0
  if(has_intercept) {
    rdX <- cbind('(Intercept)'=rep(1, nrow(rdX)), rdX)
    param_names <- c("(Intercept)", param_names)
  }
  colnames(rdX) <- param_names
  # Warn/message if there are cells without data:
  is_zero <- which(colSums(rdX) == 0)
  if(do.warn && length(is_zero)) {
    txt <- sprintf("Missing cells for: %s. ",
                   paste(param_names[is_zero], collapse = ", "))
    # warning(paste(txt, "\nInterpret type III hypotheses with care."), call.=FALSE)
    message(paste(txt, "\nInterpret type III hypotheses with care."))
  }
  rdX
}


##############################################
######## extract_contrasts_type3
##############################################
#' @importFrom MASS ginv
#' @importFrom stats terms resid lm.fit
extract_contrasts_type3 <- function(model, X=NULL) {
  # Computes contrasts for type III tests with reference to treatment contrast coding
  # X: Optional full rank design matrix in contr.treatment coding
  Terms <- terms(model)
  term_names <- attr(Terms, "term.labels")
  if(is.null(X)) {
    X <- get_model_matrix(model, type="remake", contrasts="contr.treatment")
    X <- ensure_full_rank(X)
  }
  # Get 'complete' design matrix:
  rdX <- get_rdX(model, do.warn = TRUE) # treatment contrasts
  # cols for aliased coefs should be removed in X; not in rdX.
  # This makes ginv(X) unique!
  L <- zapsmall(t(MASS::ginv(X) %*% rdX)) # basic contrast matrix
  dimnames(L) <- list(colnames(rdX), colnames(X))

  # Orthogonalize contrasts for terms which are contained in other terms:
  map <- term2colX(Terms, X)
  is_contained <- containment(model)
  # Orthogonalize higher order terms before lower order terms:
  terms_order <- attr(Terms, "order")
  orthog_order <- term_names[order(terms_order, decreasing = TRUE)]
  for(term in orthog_order) {
    # if term is contained in other terms:
    if(length(contains <- is_contained[[term]]) > 0) {
      # orthogonalize cols in L for 'term' wrt. cols that contain 'term':
      L[, map[[term]]] <-
        zapsmall(resid(lm.fit(x=L[, unlist(map[contains]), drop=FALSE],
                              y=L[, map[[term]], drop=FALSE])))
    }
  }
  # Keep rows in L corresponding to model coefficients:
  L <- L[colnames(X), , drop=FALSE]
  # Extract list of contrast matrices from L - one for each term:
  Llist <- lapply(map[term_names], function(term) t(L[, term, drop=FALSE]))
  # Keep all non-zero rows:
  lapply(Llist, function(L) L[rowSums(abs(L)) > 1e-8, , drop=FALSE])
}


##############################################
######## get_yates_contrast()
##############################################
get_yates_contrast <- function(model, which=NULL) {
  term_names <- attr(terms(model), "term.labels")
  if(is.null(which)) which <- term_names
  stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)

  var_list <- get_var_list(model)
  grid <- get_min_data(model)
  form <- formula(model)[-2]
  if(inherits(model, "lmerMod")) form <- nobars(form)
  coef_nm <- if(inherits(model, "lmerMod")) colnames(model.matrix(model)) else
    names(coef(model))[!is.na(coef(model))]
  uX <- model.matrix(form, data=grid)

  # Compute LS-means contrast:
  Llist <- lapply(which, function(term) {
    Lt <- model.matrix(formula(paste0("~ 0 + ", term)), data=grid)
    wts <- 1/colSums(Lt)
    # Lt * c(Lt %*% wts)
    # L <- diag(wts) %*% t(Lt)
    L <- t(sweep(Lt, 2, wts, "*"))
    L %*% uX
  })
  # Check estimability:
  XX <- model.matrix(terms(model), data=model.frame(model))
  # Restore contrast coding here.
  nullspaceX <- nullspace(XX)
  not_estim <- sapply(Llist, function(L)
    any(!is_estimable(L, nullspace = nullspaceX)))
  if(any(not_estim))
    warning(sprintf("Yates contrast is not uniquely defined for: %s",
                    paste(names(Llist[not_estim]), collapse = ", ")),
            call. = FALSE)

  # Make contrast for joint test of contrast among LS-means:
  lapply(Llist, function(L) {
    (t(get_trts(rownames(L))) %*% L)[, coef_nm, drop=FALSE]
  })
}


