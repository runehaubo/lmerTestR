##############################################
######## lmerModLmerTest class
##############################################

#' Represent Linear Mixed-Effects Models
#'
#' The \code{lmerModLmerTest} class extends \code{lmerMod} (which extends
#' \code{merMod}) from the \pkg{lme4}-package.
#'
#' @slot vcov_varpar a numeric matrix holding the asymptotic variance-covariance
#' matrix of the variance parameters (including sigma).
#' @slot Jac_list a list of gradient matrices (Jacobians) for the gradient of
#' the variance-covariance of beta with respect to the variance parameters,
#' where beta are the mean-value parameters available in \code{fixef(object)}.
#' @slot vcov_beta a numeric matrix holding the asymptotic variance-covariance
#' matrix of the fixed-effect regression parameters (beta).
#' @slot sigma the residual standard deviation.
#'
#' @seealso \code{\link[lme4]{lmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @author Rune Haubo B. Christensen
#' @importClassesFrom lme4 lmerMod
#'
#' @return An object of class \code{lmerModLmerTest} with slots as in
#' \code{lmerMod} objects (see \code{\link[lme4]{merMod}}) and a few
#' additional slots as described in the slots section.
lmerModLmerTest <-
  setClass("lmerModLmerTest",
           contains = c("lmerMod"),
           representation = representation(vcov_varpar = "matrix",
                                           Jac_list = "list",
                                           vcov_beta = "matrix",
                                           sigma = "numeric"))

##############################################
######## lmer()
##############################################
#' Fit Linear Mixed-Effects Models
#'
#' This function overloads \code{\link[lme4]{lmer}} from the \pkg{lme4}-package
#' (\code{lme4::lmer}) and adds a couple of slots needed for the computation of
#' Satterthwaite denominator degrees of freedom. All arguments are the same as
#' for \code{lme4::lmer} and all the usual \code{lmer}-methods work.
#'
#' For details about \code{lmer} see \code{\link[lme4]{lmer}}
#' (\code{help(lme4::lmer)}). The description of all arguments is taken
#' unedited from the \pkg{lme4}-package.
#'
#' In cases when a valid \code{lmer}-object
#' (\code{lmerMod}) is produced, but when the computations needed for
#' Satterthwaite df fails, the \code{lmerMod} object is returned - not an
#' \code{lmerModLmerTest} object.
#'
#' @inheritParams lme4::lmer
#'
#' @return an S4 object of class \code{"lmerModLmerTest"}
#' @export
#' @importFrom lme4 lmerControl
#' @importFrom methods as new
#' @seealso \code{\link[lme4]{lmer}} and \code{\link{lmerModLmerTest}}
#' @author Rune Haubo B. Christensen and Alexandra Kuznetsova for the overload
#' in \pkg{lmerTest} -- \pkg{lme4}-authors for the underlying implementation
#' in \pkg{lme4}.
#'
#' @examples
#'
#' data("sleepstudy", package="lme4")
#' m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' class(m) # lmerModLmerTest
#'
lmer <- function(formula, data = NULL, REML = TRUE,
                 control = lmerControl(), start = NULL, verbose = 0L,
                 subset, weights, na.action, offset, contrasts = NULL,
                 devFunOnly = FALSE, ...) {
  mc <- match.call()
  mc[[1L]] <- quote(lme4::lmer)
  model <- eval.parent(mc)
  # Make an lmerModLmerTest object:
  model <- as_lmerModLmerTest(model)
  # Restore the right 'call' in model:
  model@call[[1L]] <- quote(lmer)
  return(model)
}

##############################################
######## as_lmerModLmerTest()
##############################################
#' Coerce lmerMod Objects to lmerModLmerTest
#'
#' Coercing an lme4::lmer model-object (of class 'lmerMod') to a model-object
#' of class 'lmerModLmerTest' involves computing the covariance
#' matrix of the variance parameters and the gradient (Jacobian) of cov(beta)
#' with respect to the variance parameters.
#'
#' @param model and lmer model-object (of class 'lmerMod') -- the result of a
#' call to \code{lme4::lmer()}
#' @param tol tolerance for determining of eigenvalues are negative, zero or
#' positive
#'
#' @return an object of class \code{'lmerModLmerTest'} which sets the following
#' slots:
#' \item{vcov_varpar}{the asymptotic covariance matrix of the variance parameters
#' (theta, sigma).}
#' \item{Jac_list}{list of Jacobian matrices; gradients of vcov(beta) with
#' respect to the variance parameters.}
#' \item{vcov_beta}{the asymptotic covariance matrix of the fixed-effect
#' regression parameters (beta; vcov(beta)).}
#' \item{sigma}{the residual standard deviation.}
#'
#' @seealso the class definition in \code{\link{lmerModLmerTest}}) and
#' \code{\link{lmer}}
#'
#' @importFrom numDeriv hessian jacobian
#' @importFrom stats vcov update sigma
#' @importFrom lme4 getME
#'
#' @author Rune Haubo B. Christensen
#' @export
#'
#' @examples
#' m <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' bm <- lmerTest:::as_lmerModLmerTest(m)
#' slotNames(bm)
#'
#' @keywords internal
as_lmerModLmerTest <- function(model, tol=1e-8) {
  if(!inherits(model, "lmerMod"))
    stop("model not of class 'lmerMod': cannot coerce to class 'lmerModLmerTest")
  # Extract deviance function and REML indicator
  # if 'control' is not set we suppress potential message about rank deficient X:
  devfun <- if("control" %in% names(as.list(getCall(model))))
    update(model, devFunOnly=TRUE) else
      update(model, devFunOnly=TRUE,
             control=lmerControl(check.rankX = "silent.drop.cols"))
  is_reml <- getME(model, "is_REML")
  # Coerce 'lme4-model' to 'lmerModLmerTest':
  res <- as(model, "lmerModLmerTest")
  # Set relevant slots of the new model object:
  res@sigma <- sigma(model)
  res@vcov_beta <- as.matrix(vcov(model))
  varpar_opt <- unname(c(res@theta, res@sigma))
  # Compute Hessian:
  h <- numDeriv::hessian(func=devfun_vp, x=varpar_opt, devfun=devfun,
                         reml=is_reml)
  # Eigen decompose the Hessian:
  eig_h <- eigen(h, symmetric=TRUE)
  evals <- eig_h$values
  neg <- evals < -tol
  pos <- evals > tol
  zero <- evals > -tol & evals < tol
  if(sum(neg) > 0) { # negative eigenvalues
    eval_chr <- if(sum(neg) > 1) "eigenvalues" else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
    warning(sprintf("Model failed to converge with %d negative %s: %s",
                    sum(neg), eval_chr, evals_num), call.=FALSE)
  }
  # Note: we warn about negative AND zero eigenvalues:
  if(sum(zero) > 0) { # some eigenvalues are zero
    eval_chr <- if(sum(zero) > 1) "eigenvalues" else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
    warning(sprintf("Model may not have converged with %d %s close to zero: %s",
                    sum(zero), eval_chr, evals_num))

  }
  # Compute vcov(varpar):
  pos <- eig_h$values > tol
  q <- sum(pos)
  # Using the Moore-Penrose generalized inverse for h:
  h_inv <- with(eig_h, {
    vectors[, pos, drop=FALSE] %*% diag(1/values[pos], nrow=q) %*%
      t(vectors[, pos, drop=FALSE]) })
  res@vcov_varpar <- 2 * h_inv # vcov(varpar)
  # Compute Jacobian of cov(beta) for each varpar and save in list:
  Jac <- numDeriv::jacobian(func=get_covbeta, x=varpar_opt, devfun=devfun)
  res@Jac_list <- lapply(1:ncol(Jac), function(i)
    array(Jac[, i], dim=rep(length(res@beta), 2))) # k-list of jacobian matrices
  return(res)
}


##############################################
######## update.lmerModLmerTest()
##############################################
# #' @importFrom stats update
# #' @importFrom methods as
# #' @method update lmerModLmerTest
# #' @export
# #' @keywords internal
# update.lmerModLmerTest <- function(object, formula., ..., evaluate=TRUE) {
#   if(!evaluate) return(update(as(object, "lmerMod"), formula.=formula., ...,
#                               evaluate = evaluate))
#   model <- eval.parent(update(as(object, "lmerMod"), formula.=formula., ...,
#                               evaluate = evaluate))
#   return(as_lmerModLmerTest(model))
# }


##############################################
######## devfun_vp()
##############################################
#' Compute Deviance of an LMM as a Function of Variance Parameters
#'
#' This function is used for extracting the asymptotic variance-covariance matrix
#'   of the variance parameters.
#'
#' @param varpar variance parameters; \code{varpar = c(theta, sigma)}.
#' @param devfun deviance function as a function of theta only.
#' @param reml if \code{TRUE} the REML deviance is computed;
#'   if \code{FALSE}, the ML deviance is computed.
#'
#' @return the REML or ML deviance.
#' @author Rune Haubo B. Christensen
#' @keywords internal
devfun_vp <- function(varpar, devfun, reml) {
  nvarpar <- length(varpar)
  sigma2 <- varpar[nvarpar]^2
  theta <- varpar[-nvarpar]
  df_envir <- environment(devfun)
  devfun(theta) # Evaluate deviance function at varpar
  n <- nrow(df_envir$pp$V)
  # Compute deviance for ML:
  dev <- df_envir$pp$ldL2() + (df_envir$resp$wrss() + df_envir$pp$sqrL(1))/sigma2 +
    n * log(2 * pi * sigma2)
  if(!reml) return(dev)
  # Adjust if REML is used:
  RX <- df_envir$pp$RX() # X'V^{-1}X ~ crossprod(RX^{-1}) = cov(beta)^{-1} / sigma^2
  dev + 2*c(determinant(RX)$modulus) - ncol(RX) * log(2 * pi * sigma2)
}

##############################################
######## get_covbeta()
##############################################
#' Compute cov(beta) as a Function of varpar of an LMM
#'
#' At the optimum cov(beta) is available as vcov(lmer-model). This function
#' computes cov(beta) at non (RE)ML estimates of \code{varpar}.
#'
#' @inheritParams devfun_vp
#'
#' @return cov(beta) at supplied varpar values.
#' @author Rune Haubo B. Christensen
#' @keywords internal
get_covbeta <- function(varpar, devfun) {
  nvarpar <- length(varpar)
  sigma <- varpar[nvarpar] # residual std.dev.
  theta <- varpar[-nvarpar] # ranef var-par
  devfun(theta) # evaluate REML or ML deviance 'criterion'
  df_envir <- environment(devfun) # extract model environment
  sigma^2 * tcrossprod(df_envir$pp$RXi()) # vcov(beta)
}
