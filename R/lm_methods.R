# lm_methods.R - experimental methods for linear models
#
#

#' @importFrom utils as.roman
#' @importFrom stats model.matrix terms setNames formula
anova2 <- function(object) {
  Terms <- terms(object)
  data_classes <- attr(Terms, "dataClasses")
  Llist <- get_contrasts_type2(object)
  table <- do.call("rbind", lapply(Llist, contest_lm, model=object))
  ## Add row for residuals:
  df <- object$df.residual
  RSS <- sum(object$residuals^2)
  resid <- setNames(data.frame(df, RSS, RSS/df, NA, NA,
                               row.names = "Residuals"), names(table))
  table <- rbind(table, resid)
  attr(table, "heading") <-
    paste("Type", as.roman(2L), "Analysis of Variance Table\n",
          paste("\nResponse:", deparse(formula(object)[[2L]])))
  class(table) <- c("anova", "data.frame")
  table
}

#' @importFrom stats coef vcov
contest_lm <- function(L, model, eps=1e-8) {
  mk_Ftable <- function(Fvalue, ndf, ddf, RSS) {
    sigma2 <- RSS / ddf
    MS <- Fvalue * sigma2
    pvalue <- pf(q=Fvalue, df1=ndf, df2=ddf, lower.tail=FALSE)
    ans <- data.frame("Df"=ndf, "Sum Sq"=MS*ndf, "Mean Sq"=MS,
                      "F value"=Fvalue, "Pr(>F)"=pvalue, check.names = FALSE)
  }
  # if(!inherits(model, "lmerModLmerTest"))
  #   stop("'model' has to be of class lmerModLmerTest")
  if(!is.matrix(L)) L <- matrix(L, ncol=length(L))
  stopifnot(is.matrix(L), is.numeric(L),
            ncol(L) == length(beta <- coef(model)))
  if(nrow(L) == 0L) { # May happen if there are no fixed effects
    x <- numeric(0L)
    return(mk_Ftable(x, x, x, x))
  }
  # if(nrow(L) == 1L) { # 1D case:
  #   res <- contest1D(drop(L), model)
  #   return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=1L, ddf=res$df,
  #                    sigma=model@sigma))
  # } # multi-D case proceeds:
  # Compute Var(L beta) and eigen-decompose:
  vcov_beta <- vcov(model)
  VLbeta <- L %*% vcov_beta %*% t(L) # Var(contrast) = Var(Lbeta)
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
  # if(q == 1) { # 1D case:
  #   res <- contest1D(PtL, model)
  #   return(mk_Ftable(Fvalue=res[["t value"]]^2, ndf=q, ddf=res$df,
  #                    sigma=model@sigma))
  # } # multi-D case proceeds:
  # Compute t-squared values and F-value:
  t2 <- drop(PtL %*% beta)^2 / d[1:q]
  Fvalue <- sum(t2) / q
  # Compute q-list of gradients of (PtL)' cov(beta) (PtL) wrt. varpar vector:
  mk_Ftable(Fvalue, ndf=q, ddf=model$df.residual, RSS=sum(model$residuals^2))
}


get_D <- function(rdX) {
  # Compute the general Dependency matrix - aka. 'general form of estimable
  # functions' in SAS terminology.
  # This implementation uses the Moore-Penrose generalized inverse.
  D <- array(0, dim=c(ncol(rdX), ncol(rdX)))
  dimnames(D) <- list(names(attr(rdX, "param")), names(attr(rdX, "param")))
  D[1, 1] <- 1
  set <- 1
  for(i in 2:ncol(D)) {
    y <- rdX[, i]
    xx <- rdX[, set]
    b <- MASS::ginv(xx) %*% y
    resid <- y - xx %*% b
    RSS <- sum(resid^2)
    if(RSS < 1e-8) {
      b[abs(b) < 1e-8] <- 0
      D[i, set] <- b
    } else {
      D[i, i] <- 1
      set <- c(set, i)
    }
  }
  zapsmall(D)
}

#' @importFrom stats lm residuals
get_D2 <- function(rdX) {
  # Compute the general Dependency matrix - aka. 'general form of estimable
  # functions' in SAS terminology.
  # This implementation uses QR via lm().
  D <- array(0, dim=c(ncol(rdX), ncol(rdX)))
  dimnames(D) <- list(names(attr(rdX, "param")), names(attr(rdX, "param")))
  D[1, 1] <- 1
  set <- 1
  for(i in 2:ncol(D)) {
    y <- rdX[, i]
    xx <- rdX[, set]
    m <- lm(y ~ 0 + xx)
    b <- coef(m)
    RSS <- sum(residuals(m)^2)
    if(RSS < 1e-8) {
      b[abs(b) < 1e-8] <- 0
      D[i, set] <- b
    } else {
      D[i, i] <- 1
      set <- c(set, i)
    }
  }
  zapsmall(D)
}

