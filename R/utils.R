# utils.R - Utility functions

# ------- Contents: --------
#
# --- utility functions: ---
#
# qform
# rbindall
# cond
# safeDeparse
# waldCI
#

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
######## waldCI
##############################################
#' @importFrom stats qt
waldCI <- function(estimate, se, df=Inf, level=0.95) {
  stopifnot(length(level) == 1,
            is.numeric(level),
            level > 0, level < 1)
            # all(se > 0))
  alpha <- (1 - level)/2
  fac <- qt(alpha, df=df, lower.tail = FALSE)
  res <- cbind(lower = estimate - se * fac,
               upper = estimate + se * fac)
  if(!is.null(names(estimate))) rownames(res) <- names(estimate)
  res
}

# waldCI(setNames(1, "est"), .2)

