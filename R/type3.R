# type3.R - utilities for type III anova tables

##############################################
######## get_contrasts_type3
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
get_contrasts_type3 <- function(model) {
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
  Lc_list <- lapply(term_names, function(nm) {
    Lc <- contrast_type3SAS(nm, Terms, data_classes, L)
    if(!is.matrix(Lc)) Lc <- matrix(Lc, ncol=length(Lc))
    Lc[, is_coef, drop=FALSE]
  })
  names(Lc_list) <- term_names
  Lc_list
}

##############################################
######## get_rdX
##############################################
#' Compute the 'Full' Rank-Deficient Design Matrix
#'
#' Experimental - not extensively tested.
#'
#' @param model a model object; lmerMod or lmerModLmerTest
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
get_rdX <- function(model) {
  # Compute rank-deficient design-matrix X:
  #
  # model: terms(model), model.frame(model), fixef(model)
  Terms <- terms(model, fixed.only=TRUE)
  term_names <- attr(Terms, "term.labels")
  # X <- model.matrix(model)
  df <- model.frame(model)
  # Compute rank-deficient (full) design-matrix, X:
  # FIXME: Consider models without terms.
  rdXi <- lapply(term_names, function(trm)
    model.matrix(as.formula(paste0("~ 0 + ", trm)), data=df))
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
  if(any(colSums(rdX) == 0))
    warning("Missing cells for some factor levels:\n,
            Interpret type III hypotheses with care")
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
######## general_L
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
  for(colnum in colnums) {
    # pivots are the non-zero entries in L[, colnum]
    pivots <- which(abs(L[, colnum]) > eps)
    if(length(pivots) > 0) {
      # Normalize the L[pivots[1L], ] with L[pivots[1], colnum]
      L[pivots[1], ] <- L[pivots[1], ] / L[pivots[1], colnum]
      nonzeros <- setdiff(pivots, pivots[1])
      if(length(nonzeros) != 0) {
        for(nonzero in nonzeros) {
          L[nonzero, ] <- L[nonzero, ]-L[nonzero, colnum]*L[pivots[1], ]
        }
      }
      L[pivots[1], ] <- rep(0, ncol(L))
    }
  }

  # Select the non-zero rows in L:
  nums <- which(apply(L, 1, function(y) sum(abs(y))) != 0)
  L <- L[nums,]
  # If L is a vector then return, otherwise proceed to orthogonalization:
  if(is.vector(L)) return(L)

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
  L <- L[nonzero, ]
  L
}

##############################################
######## relatives
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
