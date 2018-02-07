# type12.R

###############################################################################
######## Type I anova functions
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
  L <- if(p == 1L) matrix(1L) else t(doolittle(crossprod(X))$L)
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
  # FIXME: implement keep_intercept = TRUE
}

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

###############################################################################
######## Type II anova functions
###############################################################################

##############################################
######## get_contrasts_type2
##############################################
get_contrasts_type2 <- function(X, terms, data_classes) {
  # If model has at most one term return Type I contrasts:
  term_names <- attr(terms, "term.labels")
  if(ncol(X) <= 1L || length(term_names) <= 1L)
    return(get_contrasts_type1(X, terms))

  # Compute containment:
  factor_mat <- attr(terms, "factors")
  # list: for each term 'T' a vector of terms in which 'T' is contained:
  contained_in <- setNames(lapply(term_names, function(term) {
    term_names[relatives(data_classes, term, term_names, factor_mat)]
  }), term_names)

  # Compute term asignment list: map from terms to columns in X
  has_intercept <- attr(terms, "intercept") > 0
  if(is.null(asgn <- attr(X, "assign")))
    stop("design matrix 'X' should have a non-null 'assign' attribute")
  col_terms <- if(has_intercept) c("(Intercept)", term_names)[asgn + 1] else
    term_names[asgn[asgn > 0]]
  if(!length(col_terms) == ncol(X)) # should never happen.
    stop("An error happended when computing Type II contrasts")
  term2colX <- split(seq_along(col_terms), col_terms)[unique(col_terms)]

  # Compute contrast for each term - return as named list:
  lapply(setNames(as.list(term_names), term_names), function(term) {
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

