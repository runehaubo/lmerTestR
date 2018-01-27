##############################################
######## lmerModLmerTest class
##############################################

#' Represent Linear Mixed-Effects Models
#'
#' The \code{lmerModLmerTest} class extends \code{lmerMod} (which extends
#' \code{merMod}) from the \pkg{lme4}-package.
#'
#' @slot A a numeric matrix holding the asymptotic variance-covariance matrix
#' of the variance paramters (including sigma).
#' @slot Jac_list a list of gradient matrices (Jacobians) for the gradient of
#' the variance-covariance of beta with respect to the variance parameters,
#' where beta are the mean-value parameters available in \code{fixef(object)}.
#' @slot parlist a list of parameters used internally in \pkg{lmerTestR} for
#' computing Satterthwaite's denominator degrees-of-freedom.
#'
#' @seealso \code{\link[lme4]{lmer}} and \code{\link[lme4]{merMod}}
#' @author Rune Haubo B. Christensen
#' @importClassesFrom lme4 lmerMod
#'
#' @return An object of class \code{lmerModLmerTest} with slots as in
#' \code{lmerMod} objects (see \code{\link[lme4]{merMod}}) and a few
#' additional slots as described in the slots section.
lmerModLmerTest <-
  setClass("lmerModLmerTest",
           contains = c("lmerMod"),
           representation = representation(A = "matrix",
                                           Jac_list = "list",
                                           parlist = "list"))

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
#' (\code{help(lme4::lmer)}).
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
#' @author Rune Haubo B. Christensen
#'
#' @examples
#'
#' data("sleepstudy", package="lme4")
#' m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' class(m)
#'
lmer <- function(formula, data = NULL, REML = TRUE,
                 control = lmerControl(), start = NULL, verbose = 0L,
                 subset, weights, na.action, offset, contrasts = NULL,
                 devFunOnly = FALSE, ...) {
  mc <- match.call()
  mc[[1]] <- quote(lme4::lmer)
  model <- eval.parent(mc)
  if(!inherits(model, "lmerMod")) stop("A problem occured")
  # Make an lmerModLmerTest object:
  mm <- as(model, "lmerModLmerTest")
  # Assign relevant objects to slots:
  bm <- boost_lmer(model)
  mm@parlist <- bm$parlist
  mm@A <- bm$A
  mm@Jac_list <- bm$Jac_list
  return(mm)
}

##############################################
######## boost_lmer()
##############################################
#' Boost an lmer Model-Object
#'
#' Boosting an lme4::lmer model-object involves computing the covariance
#' matrix of the variance parameters and the gradient (Jacobian) of cov(beta)
#' with respect to the variance parameters.
#'
#' @param model and lmer model-object -- the result of a call to \code{lme4::lmer()}
#' @param tol tolerance for determining of eigenvalues are negative, zero or
#' positive
#'
#' @return a list with components
#' \item{parlist}{list of parameter estimates including beta,
#' theta, sigma (vectors) and vcov(beta) (matrix)}
#' \item{A}{the asymptotic covariance matrix of the variance parameters
#' (theta, sigma)}
#' \item{Jac_list}{list of Jacobian matrices; gradients of vcov(beta) with
#' respect to the variance parameters }
#'
#' @importFrom numDeriv hessian jacobian
#' @importFrom stats vcov update sigma
#' @importFrom lme4 getME fixef
#'
#' @author Rune Haubo B. Christensen
#'
#' @examples
#' m <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' bm <- lmerTestR:::boost_lmer(m)
#' names(bm)
#'
#' @keywords internal
boost_lmer <- function(model, tol=1e-8) {
  if(!inherits(model, "lmerMod"))
    stop("'model' should be the result of lme4::lmer(), ie., an 'lmerMod'-object")
  res <- list()
  res$parlist <- list(beta=fixef(model),
                      theta=getME(model, "theta"),
                      sigma=sigma(model),
                      vcov_beta = as.matrix(vcov(model)))
  devfun <- update(model, devFunOnly=TRUE)
  is_reml <- getME(model, "is_REML")
  varpar_opt <- unname(c(res$parlist$theta, res$parlist$sigma))
  # Compute Hessian:
  h <- numDeriv::hessian(func=devfun_vp, x=varpar_opt, devfun=devfun, reml=is_reml)
  # Eigen decompose the Hessian:
  eig_h <- eigen(h, symmetric=TRUE)
  if(any(eig_h$values < -tol))
    stop("Model did not converge: Hessian has negative eigenvalues")
  if(any(abs(eig_h$values) < tol))
    warning("Model may not have converged: Hessian is singular")
  # Compute vcov(varpar):
  pos <- eig_h$values > tol
  q <- sum(pos)
  # Using the Moore-Penrose generalized inverse for h:
  h_inv <- with(eig_h, {
    vectors[, pos, drop=FALSE] %*% diag(1/values[pos], nrow=q) %*%
      t(vectors[, pos, drop=FALSE]) })
  res$A <- 2 * h_inv # vcov(varpar)
  # Compute Jacobian of cov(beta) for each varpar and save in list:
  Jac <- numDeriv::jacobian(func=get_covbeta, x=varpar_opt, devfun=devfun)
  res$Jac_list <- lapply(1:ncol(Jac), function(i)
    array(Jac[, i], dim=rep(length(res$parlist$beta), 2))) # k-list of jacobian matrices
  res
}

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
