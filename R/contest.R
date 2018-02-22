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
#' @param ddf the method for computing the denominator degrees of freedom.
#' \code{ddf="KR"} uses Kenward-Roger's method.
#' @param confint include columns for lower and upper confidence limits?
#' @param level confidence level.
#'
#' @return A \code{data.frame} with one row and columns with \code{"Estimate"},
#' \code{"Std. Error"}, \code{"t value"}, \code{"df"}, and \code{"Pr(>|t|)"}
#' (p-value). If \code{confint = TRUE} \code{"lower"} and \code{"upper"} columns
#' are included before the p-value column.
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
contest1D <- function(L, model, ddf=c("Satterthwaite", "KR"), confint=FALSE,
                      level = 0.95) {
  mk_ttable <- function(estimate, se, ddf, confint) {
    tstat <- estimate/se
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
  if(!inherits(model, "lmerModLmerTest"))
    stop("'model' has to be of class lmerModLmerTest")
  method <- match.arg(ddf)
  if(is.matrix(L)) L <- drop(L)
  stopifnot(is.numeric(L),
            length(L) == length(model@beta))
  if(length(L) == 0L) {
    o <- numeric(0L)
    return(mk_ttable(o, o, o, confint))
  }
  if(any(is.na(L))) return(mk_ttable(NA_real_, NA_real_, NA_real_, confint))
  estimate <- sum(L * model@beta) # contrast estimate
  if(method == "KR") { # Handle KR method:
    ans <- get_KR1D(L, model) # get var(contrast) and ddf
    if(!ans$error) {
      return(mk_ttable(estimate=estimate, se=sqrt(ans$var_con), ddf=ans$ddf,
                       confint))
    } else {
      warning("Unable to compute Kenward-Roger t-test: using Satterthwaite instead",
              call.=FALSE)
    }
  } # method == "Satterthwaite" proceeds:
  var_con <- qform(L, model@vcov_beta) # variance of contrast
  # Compute denominator DF:
  grad_var_con <-
    vapply(model@Jac_list, function(x) qform(L, x), numeric(1L)) # = {L' Jac L}_i
  satt_denom <- qform(grad_var_con, model@vcov_varpar) # g'Ag
  ddf <- drop(2 * var_con^2 / satt_denom) # denominator DF
  # return t-table:
  mk_ttable(estimate, sqrt(var_con), ddf, confint)
}

get_KR1D <- function(L, model) {
  # Compute var(contrast) and ddf using KR-method via the pbkrtest package
  if(!getME(model, "is_REML"))
    stop("Kenward-Roger's method is only available for REML model fits",
         call.=FALSE)
  if(!requireNamespace("pbkrtest", quietly = TRUE))
    stop("pbkrtest package required for Kenward-Roger's method",
         call.=FALSE)
  vcov_beta_adj <- try(pbkrtest::vcovAdj(model), silent=TRUE) # Adjusted vcov(beta)
  if(inherits(vcov_beta_adj, "try-error")) return(list(error=TRUE))
  var_con_adj <- qform(L, as.matrix(vcov_beta_adj)) # variance of contrast
  ddf <- try(pbkrtest::Lb_ddf(L=L, V0=model@vcov_beta,
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
#' Note: NumDF = row-rank(L) is determined automatically so row rank-deficient L
#' are allowed. One-dimensional contrasts are also allowed (L has 1 row).
#'
#' @param L a contrast matrix with nrow >= 1 and ncol ==
#' \code{length(fixef(model))}.
#' @param model a model object fitted with \code{lmer} from package
#' \pkg{lmerTestR}, i.e., an object of class \code{\link{lmerModLmerTest}}.
#' @param ddf the method for computing the denominator degrees of freedom and
#' F-statistics. \code{ddf="KR"} uses Kenward-Roger's method.
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
contestMD <- function(L, model, ddf=c("Satterthwaite", "KR"),
                      eps=sqrt(.Machine$double.eps)) {
  mk_Ftable <- function(Fvalue, ndf, ddf, sigma, Fscale=1) {
    MS <- Fvalue * sigma^2
    Fvalue <- Fvalue * Fscale
    pvalue <- pf(q=Fvalue, df1=ndf, df2=ddf, lower.tail=FALSE)
    data.frame("Sum Sq"=MS*ndf, "Mean Sq"=MS, "NumDF"=ndf, "DenDF"=ddf,
               "F value"=Fvalue, "Pr(>F)"=pvalue, check.names = FALSE)
  }
  if(!inherits(model, "lmerModLmerTest"))
    stop("'model' has to be of class lmerModLmerTest")
  if(!is.matrix(L)) L <- matrix(L, ncol=length(L))
  stopifnot(is.matrix(L), is.numeric(L),
            ncol(L) == length(model@beta))
  method <- match.arg(ddf)
  if(nrow(L) == 0L) { # May happen if there are no fixed effects
    x <- numeric(0L)
    return(mk_Ftable(x, x, x, x))
  }
  if(method == "KR") {
    if(!getME(model, "is_REML"))
      stop("Kenward-Roger's method is only available for REML model fits",
           call.=FALSE)
    if(!requireNamespace("pbkrtest", quietly = TRUE))
      stop("pbkrtest package required for Kenward-Roger's method",
           call.=FALSE)
    x <- try(pbkrtest::KRmodcomp(model, L)$test, silent = TRUE)
    if(inherits(x, "try-error")) { # Handle try-error
      warning("Unable to compute Kenward-Roger F-test: using Satterthwaite instead",
              call.=FALSE)
    } else { # return F-table if we can compute the KR F-test:
      return(mk_Ftable(Fvalue=x["FtestU", "stat"], ndf=x[1L, "ndf"],
                       ddf=x[1L, "ddf"], sigma=model@sigma,
                       Fscale=x["Ftest", "F.scaling"]))
    }
    # NOTE on the KR method:
    # It seems there is no easy way to calculate the scaling of the F-value,
    # so we will have to resort to "KRmodcomp(model, L)" for each of the k terms in
    # the anova table. This is highly ineffective since the same vcovAdj(model)
    # has to be compute k times.
  } # method == "Satterthwaite" proceeds:
  if(nrow(L) == 1L) { # 1D case:
    res <- contest1D(drop(L), model)
    return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=1L, ddf=res$df,
                     sigma=model@sigma))
  } # multi-D case proceeds:
  # Compute Var(L beta) and eigen-decompose:
  VLbeta <- L %*% model@vcov_beta %*% t(L) # Var(contrast) = Var(Lbeta)
  eig_VLbeta <- eigen(VLbeta)
  tol <- max(eps * eig_VLbeta$values[1], 0)
  pos <- eig_VLbeta$values > tol
  q <- sum(pos) # rank(VLbeta)
  P <- eig_VLbeta$vectors
  d <- eig_VLbeta$values
  if(q <= 0) { # shouldn't happen if L is a proper contrast
    x <- numeric(0L)
    return(mk_Ftable(x, x, x, x))
  }
  PtL <- crossprod(P, L)[1:q, ]
  if(q == 1) { # 1D case:
    res <- contest1D(PtL, model)
    return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=q, ddf=res$df,
                     sigma=model@sigma))
  } # multi-D case proceeds:
  # Compute t-squared values and F-value:
  t2 <- drop(PtL %*% model@beta)^2 / d[1:q]
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
