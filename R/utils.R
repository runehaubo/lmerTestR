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
# utils.R - Utility functions

# ------- Contents: --------
#
# --- utility functions: ---
#
# qform
# rbindall
# cond
# safeDeparse, deparse2
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

deparse2 <- function(x) paste(safeDeparse(x), collapse = " ")

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

