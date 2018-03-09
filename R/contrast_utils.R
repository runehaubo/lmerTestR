# contrast-utils.R - utility functions for contrasts, terms and anova

# -------- Contents: --------
#
# containment
# relatives
# doolittle
# ensure_full_rank
# get_rdX
# extract_contrasts_type3
# get_yates_contrast


##############################################
######## containment()
##############################################
containment <- function(object) { # lm or merMod
  # For all terms 'T' in object compute the terms
  # Return a list:
  #   for each term 'T' a vector of terms that contain 'T'.
  terms <- terms(object)
  data_classes <- attr(terms(object, fixed.only=FALSE), "dataClasses")
  # Note: need fixed.only for merMod objects to get dataClasses
  term_names <- attr(terms, "term.labels")
  factor_mat <- attr(terms, "factors")
  setNames(lapply(term_names, function(term) {
    term_names[relatives(data_classes, term, term_names, factor_mat)]
  }), term_names)
}


##############################################
######## relatives()
##############################################
#' Find which model terms 'Contain' a term
#'
#' Experimental - not extensively tested in this context.
#'
#' @param classes.term ?
#' @param term ?
#' @param term_names ?
#' @param factors ?
#'
#' @return ?
#' @author Lifted from old lmerTest package
#' @keywords internal
relatives <- function(classes.term, term, term_names, factors) {
  ## checks if the terms have the same number of covariates (if any)
  checkCovContain <- function(term1, term2) {
    num.numeric <- which(classes.term=="numeric")
    num.numeric.term1 <- which((num.numeric %in% which(factors[,term1]!=0))==TRUE)
    num.numeric.term2 <- which((num.numeric %in% which(factors[,term2]!=0))==TRUE)
    if((length(num.numeric.term1)>0 && length(num.numeric.term2)>0)||
       (length(num.numeric.term1)==0 && length(num.numeric.term2)==0))
      return(all(num.numeric.term2 == num.numeric.term1))
    else
      return(FALSE)
  }
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2]))) && checkCovContain(term1, term2)
  }
  if(length(term_names) == 1) return(NULL)
  which.term <- which(term == term_names)
  (1:length(term_names))[-which.term][sapply(term_names[-which.term],
                                             function(term2) is.relative(term, term2))]
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
        zapsmall(resid(lm.fit(x=L[, map[[contains]], drop=FALSE],
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


