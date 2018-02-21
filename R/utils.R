# utils.R - Utility functions


##############################################
######## qform
##############################################

#' Compute Quadratic Form
#'
#' Efficiently computes \eqn{x' A x} - or in R-notation:
#'
#' Length of \code{x} should equal the number of rows and columns of \code{A}.
#'
#' @param x a numeric vector
#' @param A a symmetric numeric matrix
#'
#' @return a numerical scalar
#' @keywords internal
qform <- function(x, A) {
  sum(x * (A %*% x)) # quadratic form: x'Ax
}

##############################################
######## rbindall
##############################################

#' \code{rbind} Multiple Objects
#'
#' @param ... objects to be \code{rbind}'ed - typically matrices or vectors
#'
#' @keywords internal
rbindall <- function(...) do.call(rbind, ...)

cbindall <- function(...) do.call(cbind, ...)

##############################################
######## cond
##############################################
cond <- function(X) with(eigen(X, only.values=TRUE), max(values) / min(values))

##############################################
######## safeDeparse
##############################################
safeDeparse <- function(expr, width.cutoff=500L, backtick = mode(expr) %in%
                          c("call", "expression", "(", "function"),
                        control = c("keepInteger","showAttributes", "keepNA"),
                        nlines = -1L) {
  deparse(expr=expr, width.cutoff=width.cutoff, backtick=backtick,
          control=control, nlines=nlines)
}

##############################################
######## containment
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

