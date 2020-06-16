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
# contest.R - contrast tests using Satterthwaites or KR ddf

# ------- Contents: --------
#
# --- Generics: ---
#
# contest
# contest1D
# contestMD
#
# --- methods: ---
#
# contest.lmerModLmerTest
# contest1D.lmerModLmerTest
# contestMD.lmerModLmerTest
# contest.lmerMod
# contest1D.lmerMod
# contestMD.lmerMod
#
# --- other exported function: ---
#
# calcSatterth
#
# --- utility functions: ---
#
# get_KR1D
# get_Fstat_ddf


##############################################
######## Generics for contest, contest1D and contestMD
##############################################
#' Generic Contrast Test Functions
#'
#' Generic functions for tests contrasts.
#'
#' @param L a contrast vector or matrix.
#' @param model a model object.
#' @param ... additional arguments passed to methods.
#'
#' @export
#' @author Rune Haubo B. Christensen
#' @seealso contest methods for \code{\link{lmer}} objects:
#' \code{\link[=contest.lmerModLmerTest]{contest}},
#' \code{\link[=contest1D.lmerModLmerTest]{contest1D}}, and
#' \code{\link[=contestMD.lmerModLmerTest]{contestMD}}.
#' @keywords internal
contest <- function(model, L, ...) UseMethod("contest")

#' @rdname contest
#' @export
contest1D <- function(model, L, ...) UseMethod("contest1D")

#' @rdname contest
#' @export
contestMD <- function(model, L, ...) UseMethod("contestMD")


##############################################
######## contest()
##############################################
#' Test of Contrasts
#'
#' Tests of vector or matrix contrasts for \code{\link{lmer}} model fits.
#'
#' If the design matrix is rank deficient, \code{lmer} drops columns for the
#' aliased coefficients from the design matrix and excludes the corresponding
#' aliased coefficients from \code{fixef(model)}. When estimability is checked
#' the original rank-deficient design matrix is recontructed and therefore
#' \code{L} contrast vectors need to include elements for the aliased
#' coefficients. Similarly when \code{L} is a matrix, its number of columns
#' needs to match that of the reconstructed rank-deficient design matrix.
#'
#' @param L a contrast vector or matrix or a list of these.
#' The \code{length}/\code{ncol} of each contrasts should equal
#' \code{length(fixef(model))}.
#' @param model a model object fitted with \code{lmer} from package
#' \pkg{lmerTest}, i.e., an object of class \code{\link{lmerModLmerTest}}.
#' @param rhs right-hand-side of the statistical test, i.e. the hypothesized
#' value (a numeric scalar).
#' @param ddf the method for computing the denominator degrees of freedom.
#' \code{ddf="Kenward-Roger"} uses Kenward-Roger's method.
#' @param confint include columns for lower and upper confidence limits? Applies
#' when \code{joint} is \code{FALSE}.
#' @param level confidence level.
#' @param joint make an F-test of potentially several contrast vectors? If
#' \code{FALSE} single DF t-tests are applied to each vector or each row of
#' contrasts matrices.
#' @param collect collect list of tests in a matrix?
#' @param check_estimability check estimability of contrasts? Only single DF
#' contrasts are checked for estimability thus requiring \code{joint = FALSE} to
#' take effect. See details section for necessary adjustments to \code{L} when
#' estimability is checked with rank deficient design matrices.
#' @param ... passed to \code{\link{contestMD}}.
#'
#' @return a \code{data.frame} or a list of \code{data.frame}s.
#' @export
#' @seealso \code{\link[=contestMD.lmerModLmerTest]{contestMD}} for multi
#' degree-of-freedom contrast tests,
#' and \code{\link[=contest1D.lmerModLmerTest]{contest1D}} for tests of
#' 1-dimensional contrasts.
#' @author Rune Haubo B. Christensen
#' @importFrom stats coef model.matrix setNames
#'
#' @examples
#'
#' data("sleepstudy", package="lme4")
#' fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
#'            sleepstudy)
#' # F-test of third coeffcients - I(Days^2):
#' contest(fm, c(0, 0, 1))
#' # Equivalent t-test:
#' contest(fm, L=c(0, 0, 1), joint=FALSE)
#' # Test of 'Days + I(Days^2)':
#' contest(fm, L=diag(3)[2:3, ])
#' # Other options:
#' contest(fm, L=diag(3)[2:3, ], joint=FALSE)
#' contest(fm, L=diag(3)[2:3, ], joint=FALSE, collect=FALSE)
#'
#' # Illustrate a list argument:
#' L <- list("First"=diag(3)[3, ], "Second"=diag(3)[-1, ])
#' contest(fm, L)
#' contest(fm, L, collect = FALSE)
#' contest(fm, L, joint=FALSE, confint = FALSE)
#' contest(fm, L, joint=FALSE, collect = FALSE, level=0.99)
#'
#' # Illustrate testing of estimability:
#' # Consider the 'cake' dataset with a missing cell:
#' data("cake", package="lme4")
#' cake$temperature <- factor(cake$temperature, ordered=FALSE)
#' cake <- droplevels(subset(cake, temperature %in% levels(cake$temperature)[1:2] &
#'                             !(recipe == "C" & temperature == "185")))
#' with(cake, table(recipe, temperature))
#' fm <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
#' fixef(fm)
#' # The coefficient for recipeC:temperature185 is dropped:
#' attr(model.matrix(fm), "col.dropped")
#' # so any contrast involving this coefficient is not estimable:
#' Lmat <- diag(6)
#' contest(fm, Lmat, joint=FALSE, check_estimability = TRUE)
#'
contest.lmerModLmerTest <- function(model, L, rhs=0, joint=TRUE, collect=TRUE, confint=TRUE,
                                    level=0.95, check_estimability=FALSE,
                                    ddf=c("Satterthwaite", "Kenward-Roger", "lme4"), ...) {
  ddf <- match.arg(ddf)
  if(!(is_list <- is.list(L))) L <- list(L)
  if(joint) {
    res <- lapply(L, function(l) contestMD(model, l, ddf=ddf, rhs=rhs, ...))
  } else { # joint is FALSE:
    if(check_estimability) {
      coef_nm <- if(inherits(model, "lmerMod")) colnames(model.matrix(model)) else
        names(coef(model))[!is.na(coef(model))]
      XX <- get_model_matrix(model, type="remake", contrasts="restore")
      keep_coef <- match(coef_nm, colnames(XX), 0L)
      nullspaceX <- nullspace(XX)
    }
    res <- lapply(L, function(l) {
      if(!is.matrix(l)) l <- matrix(l, ncol=length(l))
      if(check_estimability) {
        if(ncol(l) != ncol(XX))
          stop(sprintf("Contrast has length/ncol %i, expecting length/ncol %i when checking estimability.",
                       ncol(l), ncol(XX)))
        estim <- is_estimable(l, nullspace = nullspaceX)
        l[!estim, ] <- NA_real_ # set unestimable contrasts to NA
        l <- l[, keep_coef, drop=FALSE] # drop aliased coefs
      }
      l <- lapply(setNames(1:nrow(l), rownames(l)), function(i) l[i, ])
      rbindall(lapply(l, function(ll)
        contest1D(model, ll, rhs=rhs, ddf=ddf, confint=confint,
                  level=level)))
    })
  }
  if(collect) rbindall(res) else res
}


##############################################
######## contest1D()
##############################################
#' Contrast Tests in 1D
#'
#' Compute the test of a one-dimensional (vector) contrast in a
#' linear mixed model fitted with lmer from package \pkg{lmerTest}.
#' The contrast should specify a linear function of the
#' mean-value parameters, beta. The Satterthwaite or Kenward-Roger method is
#' used to compute the (denominator) df for the t-test.
#'
#' The t-value and associated p-value is for the hypothesis
#' \eqn{L' \beta = \mathrm{rhs}}{L' \beta = rhs} in which rhs may be non-zero
#' and \eqn{\beta} is \code{fixef(model)}.
#' The estimated value (\code{"Estimate"}) is \eqn{L' \beta} with associated
#' standard error and (optionally) confidence interval.
#'
#' @param L a numeric (contrast) vector of the same length as
#' \code{fixef(model)}.
#' @param model a model object fitted with \code{lmer} from package
#' \pkg{lmerTest}, i.e., an object of class \code{\link{lmerModLmerTest}}.
#' @param rhs right-hand-side of the statistical test, i.e. the hypothesized
#' value (a numeric scalar).
#' @param ddf the method for computing the denominator degrees of freedom.
#' \code{ddf="Kenward-Roger"} uses Kenward-Roger's method.
#' @param confint include columns for lower and upper confidence limits?
#' @param level confidence level.
#' @param ... currently not used.
#'
#' @return A \code{data.frame} with one row and columns with \code{"Estimate"},
#' \code{"Std. Error"}, \code{"t value"}, \code{"df"}, and \code{"Pr(>|t|)"}
#' (p-value). If \code{confint = TRUE} \code{"lower"} and \code{"upper"} columns
#' are included before the p-value column.
#' @export
#' @seealso \code{\link[=contest.lmerModLmerTest]{contest}} for a flexible
#' and general interface to tests of contrasts among fixed-effect parameters.
#' \code{\link[=contestMD.lmerModLmerTest]{contestMD}} is also available as a
#' direct interface for tests of multi degree-of-freedom contrast.
#' @author Rune Haubo B. Christensen
#' @importFrom stats pt
#'
#' @examples
#'
#' # Fit model using lmer with data from the lme4-package:
#' data("sleepstudy", package="lme4")
#' fm <- lmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy)
#'
#' # Tests and CI of model coefficients are obtained with:
#' contest1D(fm, c(1, 0), confint=TRUE) # Test for Intercept
#' contest1D(fm, c(0, 1), confint=TRUE) # Test for Days
#'
#' # Tests of coefficients are also part of:
#' summary(fm)
#'
#' # Illustrate use of rhs argument:
#' contest1D(fm, c(0, 1), confint=TRUE, rhs=10) # Test for Days-coef == 10
#'
#'
contest1D.lmerModLmerTest <- function(model, L, rhs=0,
                                      ddf=c("Satterthwaite", "Kenward-Roger"),
                                      confint=FALSE, level = 0.95, ...) {
  mk_ttable <- function(estimate, se, ddf) {
    tstat <- (estimate - rhs)/se
    pvalue <- 2 * pt(abs(tstat), df = ddf, lower.tail = FALSE)
    if(confint) {
      ci <- waldCI(estimate, se, ddf, level=level)
      data.frame("Estimate"=estimate, "Std. Error"=se, "df"=ddf,
                 "t value"=tstat,
                 lower=unname(ci[, "lower"]), upper=unname(ci[, "upper"]),
                 "Pr(>|t|)"=pvalue, check.names=FALSE)
    } else
      data.frame("Estimate"=estimate, "Std. Error"=se, "df"=ddf,
                 "t value"=tstat, "Pr(>|t|)"=pvalue, check.names=FALSE)
  }
  method <- match.arg(ddf)
  if(is.matrix(L)) L <- drop(L)
  stopifnot(is.numeric(L), length(L) == length(model@beta),
            is.numeric(rhs), length(rhs) == 1L)
  if(length(L) == 0L) {
    o <- numeric(0L)
    return(mk_ttable(o, o, o))
  }
  if(any(is.na(L))) return(mk_ttable(NA_real_, NA_real_, NA_real_))
  estimate <- sum(L * model@beta) # contrast estimate
  if(method == "Kenward-Roger") { # Handle KR method:
    ans <- get_KR1D(model, L) # get var(contrast) and ddf
    if(!ans$error) {
      return(mk_ttable(estimate=estimate, se=sqrt(ans$var_con), ddf=ans$ddf))
    } else {
      warning("Unable to compute Kenward-Roger t-test: using Satterthwaite instead",
              call.=FALSE)
      if(!inherits(model, "lmerModLmerTest")) model <- as_lmerModLmerTest(model)
    }
  } # method == "Satterthwaite" proceeds:
  var_con <- qform(L, model@vcov_beta) # variance of contrast
  # Compute denominator DF:
  grad_var_con <-
    vapply(model@Jac_list, function(x) qform(L, x), numeric(1L)) # = {L' Jac L}_i
  satt_denom <- qform(grad_var_con, model@vcov_varpar) # g'Ag
  ddf <- drop(2 * var_con^2 / satt_denom) # denominator DF
  # return t-table:
  mk_ttable(estimate, sqrt(var_con), ddf)
}

get_KR1D <- function(model, L) {
  # Compute var(contrast) and ddf using KR-method via the pbkrtest package
  if(!getME(model, "is_REML"))
    stop("Kenward-Roger's method is only available for REML model fits",
         call.=FALSE)
  if(!requireNamespace("pbkrtest", quietly = TRUE))
    stop("pbkrtest package required for Kenward-Roger's method",
         call.=FALSE)
  ## Add warning as faulty results have been seen with R version 3.3.2 cf https://github.com/hojsgaard/pbkrtest/issues/1
  ## It may also be related to the Matrix version: an unstated dependency in pbkrtest.
  if(getRversion() < "3.3.2")
    warning("Kenward-Roger may give faulty results with R <= 3.3.2")
  vcov_beta_adj <- try(pbkrtest::vcovAdj(model), silent=TRUE) # Adjusted vcov(beta)
  if(inherits(vcov_beta_adj, "try-error")) return(list(error=TRUE))
  var_con_adj <- qform(L, as.matrix(vcov_beta_adj)) # variance of contrast
  ddf <- try(pbkrtest::Lb_ddf(L=L, V0=vcov(model),
                              Vadj=vcov_beta_adj), silent=TRUE) # vcov_beta_adj need to be dgeMatrix!
  if(inherits(ddf, "try-error")) return(list(error=TRUE))
  list(var_con=var_con_adj, ddf=ddf, error=FALSE)
}


##############################################
######## contestMD()
##############################################
#' Multiple Degrees-of-Freedom Contrast Tests
#'
#' Compute the multi degrees-of-freedom test in a linear mixed model fitted
#' by \code{\link{lmer}}. The contrast (L) specifies a linear function of the
#' mean-value parameters, beta. Satterthwaite's method is used to compute the
#' denominator df for the F-test.
#'
#' The F-value and associated p-value is for the hypothesis
#' \eqn{L \beta = \mathrm{rhs}}{L \beta = rhs} in which rhs may be non-zero
#' and \eqn{\beta} is \code{fixef(model)}.
#'
#' Note: NumDF = row-rank(L) is determined automatically so row rank-deficient L
#' are allowed. One-dimensional contrasts are also allowed (L has 1 row).
#'
#' @param L a contrast matrix with nrow >= 1 and ncol ==
#' \code{length(fixef(model))}.
#' @param model a model object fitted with \code{lmer} from package
#' \pkg{lmerTest}, i.e., an object of class \code{\link{lmerModLmerTest}}.
#' @param rhs right-hand-side of the statistical test, i.e. the hypothesized
#' value. A numeric vector of length \code{nrow(L)} or a numeric scalar.
#' @param ddf the method for computing the denominator degrees of freedom and
#' F-statistics. \code{ddf="Kenward-Roger"} uses Kenward-Roger's method.
#' @param eps tolerance on eigenvalues to determine if an eigenvalue is
#' positive. The number of positive eigenvalues determine the rank of
#' L and the numerator df of the F-test.
#' @param ... currently not used.
#'
#' @return a \code{data.frame} with one row and columns with \code{"Sum Sq"},
#' \code{"Mean Sq"}, \code{"F value"}, \code{"NumDF"} (numerator df),
#' \code{"DenDF"} (denominator df) and \code{"Pr(>F)"} (p-value).
#' @export
#' @seealso \code{\link[=contest.lmerModLmerTest]{contest}} for a flexible and
#' general interface to tests of contrasts among fixed-effect parameters.
#' \code{\link[=contest1D.lmerModLmerTest]{contest1D}} is a direct interface for
#' tests of 1-dimensional contrasts.
#' @author Rune Haubo B. Christensen
#' @importFrom stats pf
#' @importFrom MASS ginv
#'
#' @examples
#'
#' data("sleepstudy", package="lme4")
#' fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
#'            sleepstudy)
#'
#' # Define 2-df contrast - since L has 2 (linearly independent) rows
#' # the F-test is on 2 (numerator) df:
#' L <- rbind(c(0, 1, 0), # Note: ncol(L) == length(fixef(fm))
#'            c(0, 0, 1))
#'
#' # Make the 2-df F-test of any effect of Days:
#' contestMD(fm, L)
#'
#' # Illustrate rhs argument:
#' contestMD(fm, L, rhs=c(5, .1))
#'
#' # Make the 1-df F-test of the effect of Days^2:
#' contestMD(fm, L[2, , drop=FALSE])
#' # Same test, but now as a t-test instead:
#' contest1D(fm, L[2, , drop=TRUE])
#'
contestMD.lmerModLmerTest <- function(model, L, rhs=0,
                                      ddf=c("Satterthwaite", "Kenward-Roger"),
                                      eps=sqrt(.Machine$double.eps), ...) {
  mk_Ftable <- function(Fvalue, ndf, ddf, sigma, Fscale=1) {
    MS <- Fvalue * sigma^2
    Fvalue <- Fvalue * Fscale
    pvalue <- pf(q=Fvalue, df1=ndf, df2=ddf, lower.tail=FALSE)
    data.frame("Sum Sq"=MS*ndf, "Mean Sq"=MS, "NumDF"=ndf, "DenDF"=ddf,
               "F value"=Fvalue, "Pr(>F)"=pvalue, check.names = FALSE)
  }
  if(!is.matrix(L)) L <- matrix(L, ncol=length(L))
  stopifnot(is.matrix(L), is.numeric(L),
            ncol(L) == length(model@beta))
  if(length(rhs) == 1L) rhs <- rep(rhs, nrow(L))
  stopifnot(is.numeric(rhs), length(rhs) == nrow(L))
  method <- match.arg(ddf)
  if(nrow(L) == 0L) { # May happen if there are no fixed effects
    x <- numeric(0L)
    return(mk_Ftable(x, x, x, x))
  }
  if(any(is.na(L))) return(mk_Ftable(NA_real_, NA_real_, NA_real_, NA_real_))
  if(method == "Kenward-Roger") {
    if(!getME(model, "is_REML"))
      stop("Kenward-Roger's method is only available for REML model fits",
           call.=FALSE)
    if(!requireNamespace("pbkrtest", quietly = TRUE))
      stop("pbkrtest package required for Kenward-Roger's method",
           call.=FALSE)
    if(getRversion() < "3.3.2") # See comments above.
      warning("Kenward-Roger may give faulty results with R <= 3.3.2")
    if(qr(L)$rank < nrow(L) && !all(rhs == 0))
      warning("Contrast is rank deficient and test may be affected")
    betaH <- if(all(rhs == 0)) 0 else drop(MASS::ginv(L) %*% rhs)
    x <- try(pbkrtest::KRmodcomp(model, L, betaH=betaH)$test, silent = TRUE)
    if(inherits(x, "try-error")) { # Handle try-error
      warning("Unable to compute Kenward-Roger F-test: using Satterthwaite instead",
              call.=FALSE)
      if(!inherits(model, "lmerModLmerTest")) model <- as_lmerModLmerTest(model)
    } else { # return F-table if we can compute the KR F-test:
      return(mk_Ftable(Fvalue=x["FtestU", "stat"], ndf=x[1L, "ndf"],
                       ddf=x[1L, "ddf"], sigma=sigma(model),
                       Fscale=x["Ftest", "F.scaling"]))
    }
    # NOTE on the KR method:
    # It seems there is no easy way to calculate the scaling of the F-value,
    # so we will have to resort to "KRmodcomp(model, L)" for each of the k terms in
    # the anova table. This is highly ineffective since the same vcovAdj(model)
    # has to be compute k times.
  } # method == "Satterthwaite" proceeds:
  if(nrow(L) == 1L) { # 1D case:
    res <- contest1D(model, drop(L), rhs=rhs, confint=FALSE)
    return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=1L, ddf=res$df,
                     sigma=model@sigma))
  } # multi-D case proceeds:
  beta <- model@beta
  # Adjust beta for rhs:
  if(!all(rhs == 0)) beta <- beta - drop(MASS::ginv(L) %*% rhs)
  # Compute Var(L beta) and eigen-decompose:
  VLbeta <- L %*% model@vcov_beta %*% t(L) # Var(contrast) = Var(Lbeta)
  eig_VLbeta <- eigen(VLbeta)
  P <- eig_VLbeta$vectors
  d <- eig_VLbeta$values
  tol <- max(eps * d[1], 0)
  pos <- d > tol
  q <- sum(pos) # rank(VLbeta)
  if(q < nrow(L) && !all(rhs == 0))
    warning("Contrast is rank deficient and test may be affected")
  if(q <= 0) { # shouldn't happen if L is a proper contrast
    x <- numeric(0L)
    return(mk_Ftable(x, x, x, x))
  }
  PtL <- crossprod(P, L)[1:q, ]
  if(q == 1) { # 1D case:
    res <- contest1D(model, PtL, rhs=rhs[1L], confint=FALSE)
    return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=q, ddf=res$df,
                     sigma=model@sigma))
  } # multi-D case proceeds:
  # Compute t-squared values and F-value:
  t2 <- drop(PtL %*% beta)^2 / d[1:q]
  Fvalue <- sum(t2) / q
  # Compute q-list of gradients of (PtL)' cov(beta) (PtL) wrt. varpar vector:
  grad_PLcov <- lapply(1:q, function(m) {
    vapply(model@Jac_list, function(J) qform(PtL[m, ], J), numeric(1L))
  })
  # Compute degrees of freedom for the q t-statistics:
  nu_m <- vapply(1:q, function(m) {
    2*(d[m])^2 / qform(grad_PLcov[[m]], model@vcov_varpar) }, numeric(1L)) # 2D_m^2 / g'Ag
  # Compute ddf for the F-value:
  ddf <- get_Fstat_ddf(nu_m, tol=1e-8)
  mk_Ftable(Fvalue, ndf=q, ddf=ddf, sigma=model@sigma)
}


##############################################
######## get_Fstat_ddf()
##############################################
#' Compute denominator df for F-test
#'
#' From a vector of denominator df from independent t-statistics (\code{nu}),
#' the denominator df for the corresponding F-test is computed.
#'
#' Note that if any \code{nu <= 2} then \code{2} is returned. Also, if all nu
#' are within tol of each other the simple average of the nu-vector is returned.
#' This is to avoid downward bias.
#'
#' @param nu vector of denominator df for the t-statistics
#' @param tol tolerance on the consequtive differences between elements of nu to
#  determine if mean(nu) should be returned
#'
#' @author Rune Haubo B. Christensen
#'
#' @return the denominator df; a numerical scalar
#' @keywords internal
get_Fstat_ddf <- function(nu, tol=1e-8) {
  # Computes denominator df for an F-statistic that is derived from a sum of
  # squared t-statistics each with nu_m degrees of freedom.
  #
  # nu : vector of denominator df for the t-statistics
  # tol: tolerance on the consequtive differences between elements of nu to
  #      determine if mean(nu) should be returned.
  #
  # Result: a numeric scalar
  #
  # Returns nu if length(nu) == 1. Returns mean(nu) if all(abs(diff(nu)) < tol;
  # otherwise ddf appears to be downward biased.
  fun <- function(nu) {
    if(any(nu <= 2)) 2 else {
      E <- sum(nu / (nu - 2))
      2 * E / (E - (length(nu))) # q = length(nu) : number of t-statistics
    }
  }
  stopifnot(length(nu) >= 1,
            # all(nu > 0), # returns 2 if any(nu < 2)
            all(sapply(nu, is.numeric)))
  if(length(nu) == 1L) return(nu)
  if(all(abs(diff(nu)) < tol)) return(mean(nu))
  if(!is.list(nu)) fun(nu) else vapply(nu, fun, numeric(1L))
}


##############################################
######## calcSatterth()
##############################################
#' @rdname contestMD.lmerModLmerTest
#' @export
calcSatterth <- function(model, L) {
  stopifnot(inherits(model, "lmerMod"))
  if(!inherits(model, "lmerModLmerTest")) {
    message("Coercing model to class 'lmerModLmerTest'")
    model <- as_lmerModLmerTest(model)
    if(!inherits(model, "lmerModLmerTest"))
      stop("Failed to coerce model to class 'lmerModLmerTest'")
  }
  x <- contestMD(model, L)
  list("denom"=x[["DenDF"]], "Fstat"=as.matrix(x[["F value"]]),
       "pvalue"=as.matrix(x[["Pr(>F)"]]), "ndf"=x[["NumDF"]])
}
# m <- lmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy)
# L <- cbind(0,1) ## specify contrast vector
# contestMD(m, L)
# calcSatterth(m, L)


##############################################
######## lmerMod methods for contest, contest1D and contestMD
##############################################
#' @rdname contest.lmerModLmerTest
#' @export
contest.lmerMod <- function(model, L, rhs=0, joint=TRUE, collect=TRUE, confint=TRUE,
                            level=0.95, check_estimability=FALSE,
                            ddf=c("Satterthwaite", "Kenward-Roger", "lme4"), ...) {
  ddf <- match.arg(ddf)
  # For Satterthwaite we need to compute stuff - not for K-R:
  if(ddf == "Satterthwaite") model <- as_lmerModLmerTest(model)
  # Use lmerModLmerTest method:
  eval.parent(contest.lmerModLmerTest(model, L=L, joint=joint, collect=collect,
                                      confint=confint, level=level,
                                      check_estimability=check_estimability,
                                      ddf=ddf, rhs=rhs, ...))
}

#' @rdname contest1D.lmerModLmerTest
#' @export
contest1D.lmerMod <- function(model, L, rhs=0,
                              ddf=c("Satterthwaite", "Kenward-Roger"),
                              confint=FALSE, level = 0.95, ...) {
  ddf <- match.arg(ddf)
  # For Satterthwaite we need to compute stuff - not for K-R:
  if(ddf == "Satterthwaite") model <- as_lmerModLmerTest(model)
  # Use lmerModLmerTest method:
  eval.parent(contest1D.lmerModLmerTest(model, L=L, rhs=rhs, ddf=ddf,
                                        confint=confint, level=level))
}

#' @rdname contestMD.lmerModLmerTest
#' @export
contestMD.lmerMod <- function(model, L, rhs=0,
                              ddf=c("Satterthwaite", "Kenward-Roger"),
                              eps=sqrt(.Machine$double.eps), ...) {
  ddf <- match.arg(ddf)
  # For Satterthwaite we need to compute stuff - not for K-R:
  if(ddf == "Satterthwaite") model <- as_lmerModLmerTest(model)
  # Use lmerModLmerTest method:
  eval.parent(contestMD.lmerModLmerTest(model, L=L, rhs=rhs, ddf=ddf, eps=eps))
}
