# contest.R - contrast tests using Satterthwaites df

# contest <- function(L, model, joint=FALSE) {
#   # Test a contrast L vector/matrix in model
#   #
#   # L - vector
#   # L - matrix (joint = c(TRUE, FALSE))
#   # L - list of vectors
#   # L - list of matrices (joint = c(TRUE, FALSE))
#
# }

##############################################
######## contest1D()
##############################################
#' Contrast Tests in 1D
#'
#' Compute the test of a one-dimensional (vector) contrast in a
#' linear mixed model fitted with lmer from package \pkg{lmerTestR}.
#' The contrast should specify a linear function of the
#' mean-value parameters, beta. Satterthwaite's method is used to compute the
#' denominator df for the t-test.
#'
#' @param L a numeric (contrast) vector of the same length as
#' \code{fixef(model)}.
#' @param model a model object fitted with \code{lmer} from package
#' \pkg{lmerTestR}, i.e., an object of class \code{\link{lmerModLmerTest}}.
#'
#' @return A \code{data.frame} with one row and columns with \code{"Estimate"},
#' \code{"Std. Error"}, \code{"t value"}, \code{"df"}, and \code{"Pr(>|t|)"}
#' (p-value)
#' @export
#' @seealso \code{\link{contestMD}} for multi degree-of-freedom contrast tests.
#' @author Rune Haubo B. Christensen
#' @importFrom stats pt
#'
#' @examples
#'
#' # Fit model using lmer with data from the lme4-package:
#' fm <- lmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy)
#' # Note that summary do not contain any tests/p-values:
#' coef(summary(fm))
#' # We can get these tests with:
#' contest1D(c(1, 0), fm) # Test for Intercept
#' contest1D(c(0, 1), fm) # Test for Days
#'
contest1D <- function(L, model) {
  if(is.matrix(L)) L <- drop(L)
  stopifnot(is.numeric(L),
            length(L) == length(model@parlist$beta))
  estimate <- sum(L * model@parlist$beta)
  grad_var_con <-
    vapply(model@Jac_list, function(x) qform(L, x), numeric(1L)) # = {L' Jac L}_i
  var_con <- qform(L, model@parlist$vcov_beta)
  se.estimate <- sqrt(var_con)
  satt_denom <- qform(grad_var_con, model@A) # g'Ag
  ddf <- drop(2 * var_con^2 / satt_denom) # denominator DF
  tstat <- estimate/se.estimate
  pvalue <- 2 * pt(abs(tstat), df = ddf, lower.tail = FALSE)
  data.frame("Estimate"=estimate, "Std. Error"=se.estimate, "df"=ddf,
             "t value"=tstat, "Pr(>|t|)"=pvalue, check.names=FALSE)
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
#' Note: NumDF = row-rank(L) is determined automatically so row rank-deficient L
#' are allowed. One-dimensional contrasts are also allowed (L has 1 row).
#'
#' @param L a contrast matrix with nrow >= 1 and ncol ==
#' \code{length(fixef(model))}.
#' @param model a model object fitted with \code{lmer} from package
#' \pkg{lmerTestR}, i.e., an object of class \code{\link{lmerModLmerTest}}.
#' @param eps tolerance on eigenvalues to determine if an eigenvalue is
#' positive. The number of positive eigenvalues determine the rank of
#' L and the numerator df of the F-test.
#'
#' @return a \code{data.frame} with one row and columns with \code{"Sum Sq"},
#' \code{"Mean Sq"}, \code{"F value"}, \code{"NumDF"} (numerator df),
#' \code{"DenDF"} (denominator df) and \code{"Pr(>F)"} (p-value).
#' @export
#' @seealso \code{\link{contest1D}} for tests of 1-dimensional contrasts.
#' @author Rune Haubo B. Christensen
#' @importFrom stats pf
#'
#' @examples
#'
#' fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
#'            sleepstudy)
#' # Define 2-df contrast - since L has 2 (linearly independent) rows
#' # the F-test is on 2 (numerator) df:
#' L <- rbind(c(0, 1, 0), # Note: ncol(L) == length(fixef(fm))
#'            c(0, 0, 1))
#' # Make the 2-df F-test of any effect of Days:
#' contestMD(L, fm)
#' # Make the 1-df F-test of the effect of Days^2:
#' contestMD(L[2, , drop=FALSE], fm)
#' # Same test, but now as a t-test instead:
#' contest1D(L[2, , drop=TRUE], fm)
#'
contestMD <- function(L, model, eps=sqrt(.Machine$double.eps)) {
  mk_Ftable <- function(Fvalue, ndf, ddf, sigma) {
    MS <- Fvalue * sigma^2
    pvalue <- pf(q=Fvalue, df1=ndf, df2=ddf, lower.tail=FALSE)
    data.frame("Sum Sq"=MS*ndf, "Mean Sq"=MS, "NumDF"=ndf, "DenDF"=ddf,
               "F value"=Fvalue, "Pr(>F)"=pvalue, check.names = FALSE)
  }
  stopifnot(is.matrix(L), is.numeric(L),
            ncol(L) == length(model@parlist$beta),
            nrow(L) >= 1)
  if(nrow(L) == 1) { # 1D case:
    res <- contest1D(drop(L), model)
    return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=1L, ddf=res$df,
                     sigma=model@parlist$sigma))
  } # multi-D case proceeds:
  # Compute Var(L beta) and eigen-decompose:
  VLbeta <- L %*% model@parlist$vcov_beta %*% t(L) # Var(contrast) = Var(Lbeta)
  eig_VLbeta <- eigen(VLbeta)
  tol <- max(eps * eig_VLbeta$values[1], 0)
  pos <- eig_VLbeta$values > tol
  q <- sum(pos) # rank(VLbeta)
  P <- eig_VLbeta$vectors
  d <- eig_VLbeta$values
  PtL <- crossprod(P, L)[1:q, ]
  if(q <= 0) { # shouldn't happen if L is a proper contrast
    return(mk_Ftable(NA, NA, NA, NA))
  }
  if(q == 1) { # 1D case:
    res <- contest1D(PtL, model)
    return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=q, ddf=res$df,
                     sigma=model@parlist$sigma))
  } # multi-D case proceeds:
  # Compute t-squared values and F-value:
  t2 <- drop(PtL %*% model@parlist$beta)^2 / d[1:q]
  Fvalue <- sum(t2) / q
  # Compute q-list of gradients of (PtL)' cov(beta) (PtL) wrt. varpar vector:
  grad_PLcov <- lapply(1:q, function(m) {
    vapply(model@Jac_list, function(J) qform(PtL[m, ], J), numeric(1L))
  })
  # Compute degrees of freedom for the q t-statistics:
  nu_m <- vapply(1:q, function(m) {
    2*(d[m])^2 / qform(grad_PLcov[[m]], model@A) }, numeric(1L)) # 2D_m^2 / g'Ag
  # Compute ddf for the F-value:
  ddf <- get_Fstat_ddf(nu_m, tol=1e-8)
  mk_Ftable(Fvalue, ndf=q, ddf=ddf, sigma=model@parlist$sigma)
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
