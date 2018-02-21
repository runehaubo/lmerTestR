# contrast-utils.R - utility functions for contrasts, terms and anova

# Functions:

# containment
# relatives
# doolittle
# get_rdX
# general_L
# term2colX
# need_yates
# no_yates
# numeric_terms
# get_yates_contrast
# contrast_type3SAS


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
#' @author Lifted from \code{lmerTest:::makeContrastType3SAS}
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
######## get_rdX()
##############################################
#' Compute the 'Full' Rank-Deficient Design Matrix
#'
#' Experimental - not extensively tested.
#'
#' @param model a model object; lmerMod or lmerModLmerTest.
#' @param do.warn throw a warning if there is no data for some factor
#' combinations.
#' @param use_contrasts if \code{TRUE} contrasts coding from the original design
#' matrix will be used. If \code{FALSE} the settings in \code{options()$contrasts}
#' will be used.
#'
#' @return the rank-deficien design matrix, \code{rdX} with attributes
#' \item{param}{a vector of _parameters_ in the over-parameterized model
#' notation with length equal to the number of columns in \code{rdX}.}
#' \item{is_coef}{numeric vector with indices on the columns in \code{rdX}
#' which corresponds to the columns in the original full-rank \code{X} and
#' therefore also to the coefficient vector.}
#' @author Rune Haubo B. Christensen
#' @keywords internal
#'
#' @importFrom stats as.formula model.frame terms model.matrix setNames
#' @importFrom lme4 fixef
get_rdX <- function(model, do.warn=TRUE, use_contrasts=FALSE) {
  # Compute rank-deficient design-matrix X:
  #
  # model: terms(model), model.frame(model), fixef(model)
  Terms <- terms(model, fixed.only=TRUE)
  term_names <- attr(Terms, "term.labels")
  # X <- model.matrix(model)
  df <- model.frame(model)
  # Compute rank-deficient (full) design-matrix, X:
  # FIXME: Consider models without terms.
  if(!use_contrasts) {
    rdXi <- lapply(term_names, function(trm)
      model.matrix(as.formula(paste0("~ 0 + ", trm)), data=df))
  } else {
    contrasts <- attr(model.matrix(model), "contrasts")
    rdXi <- lapply(term_names, function(trm) {
      form <- as.formula(paste0("~ 0 + ", trm))
      Contrast <- if(!is.null(contrasts) && trm %in% names(contrasts))
        contrasts[trm] else NULL
      model.matrix(form, contrasts.arg = Contrast, data=df)
    })
  }
  # FIXME: Appears not to do the right thing for interactions with ordered factors.
  rdX <- do.call(cbind, rdXi)
  param_names <- unlist(lapply(rdXi, colnames))
  colnames(rdX) <- rep(term_names, vapply(rdXi, ncol, integer(1L)))
  # Potentially add intercept:
  has_intercept <- attr(Terms, "intercept") != 0
  if(has_intercept) {
    rdX <- cbind('(Intercept)'=rep(1, nrow(rdX)), rdX)
    param_names <- c("(Intercept)", param_names)
  }
  # Warn if there are cells without data:
  if(do.warn && any(colSums(rdX) == 0))
    warning("Missing cells for some factor levels: ",
            "Interpret type III hypotheses with care", call.=FALSE)
  # FIXME: mention name of offending factor-level in warning.
  # Add param and is_coef attributes to rdX and return:
  param <- setNames(numeric(ncol(rdX)), param_names)
  coef <- fixef(model)
  is_coef <- which(names(param) %in% names(coef))
  param[is_coef] <- coef
  attr(rdX, "param") <- param
  attr(rdX, "is_coef") <- is_coef
  rdX
}

##############################################
######## general_L()
##############################################
#' Compute a 'General Set of Estimable Functions'
#'
#' Computes the so-called L-matrix (in SAS terminology).
#' Experimental - not extensively tested.
#'
#' @param rdX a 'full' rank-deficient design matrix.
#' @param is_coef an indicator vector of the same length as \code{fixef(model)}
#' which indicates which columns in \code{rdX} correspond \code{fixef(model)}.
#'
#' @return The L-matrix
#'
#' @author Rune Haubo B. Christensen
#' @keywords internal
general_L <- function(rdX, is_coef) {
  # Compute a/the 'general set of estimable functions' (the L-matrix)
  #
  # rdX: rank-deficient design-matrix X with pp columns, where pp is the
  #   number of parameters in the over-parameterized model.
  # is_coef: indicator vector: which of pp cols(rdX) correspond to the model
  #   coefficients.
  xtx <- crossprod(rdX)
  g2 <- array(0, dim=dim(xtx))
  g2[is_coef, is_coef] <- solve(xtx[is_coef, is_coef])
  g2[abs(g2) < 1e-10] <- 0
  #general set of estimable function
  L <- g2 %*% xtx
  L[abs(L) < 1e-6] <- 0
  L
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

  # Make contrast for joint test of Compute contrast among LS-means:
  lapply(Llist, function(L) {
    (t(get_trts(rownames(L))) %*% L)[, coef_nm, drop=FALSE]
  })
}


##############################################
######## contrast_type3SAS
##############################################
#' Compute Type III Contrast Matrix for a Model Term
#'
#' Experimental - not extensively tested.
#'
#' Compute the type III qi x pp contrast matrix for 'term' where qi is the
#' numerater df for 'term' and pp is the number of parameters in the
#' over-parameterized model notation (and also equals the number of columns in
#' the 'full' rank-deficient design matrix). if qi=1 the result is a qi-vector.
#'
#' @param term the name of a model term; one of
#' \code{attr(terms(model), "term.labels")}.
#' @param terms a terms object (using \code{fixed.only=TRUE}).
#' @param data_classes \code{attr(terms(model, fixed.only=FALSE), "dataClasses")}.
#' @param L the pp x pp matrix with pp-p zero-rows indicating the relation between the
#' pp parameters and p coefficients.
#' @param eps tolererance for non-zero pivots
#'
#' @return a pp-vector or qi x pp matrix with the contrast for the term.
#' @author Rune Haubo B. Christensen based on
#' \code{lmerTest:::makeContrastType3SAS}
#' @keywords internal
contrast_type3SAS <- function(term, terms, data_classes, L, eps=1e-8) {
  # Apply rule 1 (Goodnight 1976)
  #
  # Find all terms that contain term_name:
  factor_mat <- attr(terms,"factors")
  term_names <- attr(terms, "term.labels")
  term_cols <- which(colnames(L) == term) # cols in L for the term
  # Terms which are relatives of 'name'/'term'
  num.relate <- relatives(data_classes, term, term_names, factor_mat)
  # related_terms <- term_names[num.relate]
  # related_cols <- ...
  # unrelated_cols <- ...

  # FIXME: Describe this:
  if(length(num.relate) == 0)
    colnums <- setdiff(1:ncol(L), term_cols)
  if(length(num.relate) > 0) {
    cols.contain <- NULL
    for(i in 1:length(num.relate))
      # cols.contain are the cols in rdX for terms which contain 'name'
      cols.contain <- c(cols.contain,
                        which(colnames(L) == term_names[num.relate[i]]))
    # colnums are the cols in rdX which are unrelated to 'name'
    colnums <- setdiff(1:ncol(L), c(term_cols, cols.contain))
  }

  # For each column among those which are unrelated 'name':
  # Seems to orthogonalize all terms that are not T and does not contain T
  #  relative to T:
  for(colnum in colnums) {
    # pivots are the non-zero entries in L[, colnum]
    pivots <- which(abs(L[, colnum]) > eps)
    if(length(pivots) > 0) {
      # Normalize the L[pivots[1L], ] with L[pivots[1], colnum]
      L[pivots[1], ] <- L[pivots[1], ] / L[pivots[1], colnum]
      nonzeros <- setdiff(pivots, pivots[1])
      if(length(nonzeros) > 0) {
        for(nonzero in nonzeros) {
          L[nonzero, ] <- L[nonzero, ]-L[nonzero, colnum]*L[pivots[1], ]
        }
      }
      L[pivots[1], ] <- rep(0, ncol(L))
    }
  }

  # Select the non-zero rows in L:
  nums <- which(apply(L, 1, function(y) sum(abs(y))) != 0)
  L <- L[nums, , drop=FALSE]
  # If L is a vector then return, otherwise proceed to orthogonalization:
  if(nrow(L) <= 1L) return(L)

  # Orthogonalization:
  if(length(term_cols) > 1)
    # Determine the rows of L which are _indirectly_ related to 'name', i.e.
    # those for which only 'related' terms contribute:
    zero.rows <- which(apply(L[, term_cols], 1, function(y) sum(abs(y))) == 0)
  else
    zero.rows <- which(L[, term_cols] == 0)

  for(zero.row in zero.rows) { # for each 'related' row in L:
    w <- L[zero.row, ]
    for(i in setdiff(1:nrow(L), zero.row)) { # for all other rows:
      if(sum(abs(L[i, ])) != 0) # if the row is non-zero:
        # orthogonalize the i'th row relative to w:
        L[i, ] <- L[i, ] - ((w %*% L[i, ]) / (w %*% w)) %*% w
    }
    L[zero.row, ] <- rep(0, ncol(L)) # set 'w' to zero
  }
  # zero small entries:
  L[abs(L) < 1e-6] <- 0
  # Keep non-zero rows:
  nonzero <- which(apply(L, 1, function(y) sum(abs(y))) != 0)
  L <- L[nonzero, , drop=FALSE]
  L
}


