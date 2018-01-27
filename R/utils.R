# utils.R - Utility functions

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

##############################################
######## cond
##############################################
cond <- function(X) with(eigen(X, only.values=TRUE), max(values) / min(values))

