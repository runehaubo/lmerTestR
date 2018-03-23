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
# estimability.R - functions for assessing model estimability

# ------- Contents: --------
#
# is_estimable
# nullspace
#


##############################################
######## is_estimable
##############################################
#' Estimability of Contrasts
#'
#' Computes the estimability of a vector or matrix of contrasts (i.e. linear
#' functions of the coefficients) from the nullspace of a design matrix or
#' potentially directly from the design matrix.
#'
#' @param contrast a numeric matrix where each row is a contrast vector for
#' which estimability is computed. The matrix should have as many columns as
#' there are columns in the design matrix (which equals the number of
#' coefficients). If \code{contrast} is a vector it is coerced to a matrix.
#' @param nullspace the nullspace of the design matrix.
#' @param X design matrix.
#' @param tol tolerance for determining if a contrast is orthogonal to the
#   nullspace.
#'
#' @return a logical vector of length \code{nrow(contrast)} determining if each
#' contrast is estimable
#' @importFrom stats setNames
#' @keywords internal
#' @seealso \code{\link{nullspace}}
#'
#' @author Rune Haubo B. Christensen
#' @keywords internal
#' @examples
#'
#' # FIXME: We need some examples here
#'
is_estimable <- function(contrast, nullspace=NULL, X=NULL,
                         tol=sqrt(.Machine$double.eps)) {
  if(!is.matrix(contrast)) contrast <- matrix(contrast, ncol=length(contrast))
  N <- if(!is.null(nullspace)) { # get nullspace
    nullspace
  } else if(!is.null(X)) {
    nullspace(X)
  } else {
    stop("Need non-null 'nullspace' or 'X' to compute estimability")
  }
  if(ncol(contrast) != nrow(N))
    stop(sprintf("'contrast' has %i columns: expecting %i columns",
                 ncol(contrast), nrow(N)))
  # Determine estimability:
  res <- if(length(N) == 0) rep(TRUE, nrow(contrast)) else
    c(abs(rowSums(contrast %*% N)) < tol)
  setNames(res, rownames(contrast))
}
#
# XX <- model.matrix(terms(model), data=model.frame(model))
# nullspaceX <- nullspace(XX)
# is_estimable(Llist$DAY, nullspaceX)
# is_estimable(c(Llist$DAY[1, ]), nullspaceX)
# is_estimable(Llist$DAY, X=XX)
# NCOL(0:1)
#
# X <- model.matrix(model)
# str(Llist$DAY[, -9] %*% nullspace(X))
# is_estimable(Llist$DAY[, -9], X=X)
# is_estimable(0:1, X=X)
# contrast <- 0:1
# nrow(matrix(0:1, ncol=2))
# rep(TRUE, 1)
#
# length(Llist$DAY[, -9] %*% nullspace(X))
# apply(Llist$DAY[, -9] %*% nullspace(X), 1, length)
# length(nullspace(X))

##############################################
######## nullspace
##############################################
#' Nullspace
#'
#' Compute the (right or left) nullspace of matrix using a (semi-complete)
#' Singular Value Decomposition.
#'
#' This implementation is fastest on matrices with more rows
#' than columns such as a typical design matrix for a linear model.
#'
#' @param A a numeric matrix.
#' @param type \code{"right"} (default) gives is the standard nullspace,
#' \code{"left"} gives left nullspace of \code{A}.
#' @param tol tolerance multiple of the first singular value to determine if
#' subsequent singular values are (sufficiently) positive to be determined
#' greater than zero.
#'
#' @return a matrix with as many rows as there are columns in \code{A}. The
#' number of columns (which may be zero) determine the dimensionality of the
#' nullspace of \code{A}.
#' @author Rune Haubo B. Christensen
#'
#' @keywords internal
#' @examples
#'
#' # FIXME: We need some examples here
#'
nullspace <- function(A, type = c("right", "left"),
                      tol=sqrt(.Machine$double.eps)) {
  # Compute the right (standard and default) or left null space of a matrix A.
  # using SVD.
  type <- match.arg(type)
  if(type == "left") return(nullspace(t(A), type="right", tol=tol))
  if(length(A) == 0L) return(matrix(numeric(0L))) # length(A) == 0 if any(dim(A) == 0)
  svdA <- svd(A, nv = ncol(A))
  tol <- 1e-8
  positive <- svdA$d > max(tol * svdA$d[1L], 0)
  rank <- sum(positive)
  set <- if(rank == 0) 1:ncol(A) else -(1:rank)
  svdA$v[, set, drop=FALSE]
}


