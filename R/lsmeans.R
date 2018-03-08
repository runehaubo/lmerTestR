

##############################################
######## ls_means()
##############################################
#' LS-means for lmerTest Model Fits
#'
#' Computes LS-means or pairwise differences of LS-mean for all factors in a
#' linear mixed model. \code{lsmeansLT} is provided as an alias for
#' \code{ls_means} for backward compatibility.
#'
#' Confidence intervals and p-values are based on the t-distribution using
#' degrees of freedom based on Satterthwaites or Kenward-Roger methods.
#'
#' LS-means is SAS terminology for predicted/estimated marginal means, i.e. means
#' for levels of factors which are averaged over the levels of other factors in
#' the model. A flat (i.e. unweighted) average is taken which gives equal weight
#' to all levels of each of the other factors. Numeric/continuous variables are
#' set at their mean values. See \pkg{emmeans} package
#' for more options and greater flexibility.
#'
#' LS-means contrasts are checked for estimability and unestimable contrasts appear
#' as \code{NA}s in the resulting table.
#'
#' LS-means objects (of class \code{"ls_means"} have a print method).
#'
#' @param model a model object fitted with \code{\link{lmer}} (of class
#' \code{"lmerModLmerTest"}).
#' @param which optional character vector naming factors for which LS-means should
#' be computed. If \code{NULL} (default) LS-means for all factors are computed.
#' @param level confidence level.
#' @param ddf method for computation of denominator degrees of freedom.
#' @param pairwise compute pairwise differences of LS-means instead?
#' @param ... currently not used.
#'
#' @return An LS-means table in the form of a \code{data.frame}. Formally an object
#' of class \code{c("ls_means", "data.frame")} with a number of attributes set.
#' @author Rune Haubo B. Christensen and Alexandra Kuznetsova
#' @seealso \code{\link{show_contrasts}} for display of the underlying LS-means
#' contrasts.
#' @export
#'
#' @examples
#'
#' # Get data and fit model:
#' data("cake", package="lme4")
#' model <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
#'
#' # Compute LS-means:
#' ls_means(model)
#'
#' # Compute pairwise differences of LS-means for each factor:
#' ls_means(model, pairwise=TRUE)
#' difflsmeans(model) # Equivalent.
#'
ls_means.lmerModLmerTest <- function(model, which=NULL, level=0.95,
                                    ddf=c("Satterthwaite", "Kenward-Roger"),
                                    pairwise=FALSE, ...) {
  ddf <- match.arg(ddf)
  Llist <- lsmeans_contrasts(model, which=which)
  coef_nm <- if(inherits(model, "lmerMod")) colnames(model.matrix(model)) else
    names(coef(model))[!is.na(coef(model))]
  # Need nullspace of _remade_ model matrix to check estimability:
  XX <- get_model_matrix(model, type="remake", contrasts="restore")
  nullspaceX <- nullspace(XX)

  # Pairwise differences:
  if(pairwise == TRUE) # Adjust contrasts to compute pairwise diffs:
    Llist <- lapply(Llist, function(L)
      crossprod(as.matrix(get_pairs(rownames(L))), L))

  # Compute LS-means:
  if(length(Llist) == 0) {
    means <- contest1D(rep(NA_real_, length(coef_nm)), ddf=ddf, model=model,
                       confint=TRUE, level=level)[0L, , drop=FALSE]
  } else
    means <- rbindall(lapply(names(Llist), function(var) {
      L <- Llist[[var]]
      # Check estimability before computing the contrast:
      estim <- is_estimable(L, nullspace = nullspaceX)
      L[!estim, ] <- NA_real_ # set unestimable contrasts to NA
      L <- L[, coef_nm, drop=FALSE] # drop aliased coefs
      # Evaluate contrasts:
      tab <- rbindall(lapply(1:nrow(L), function(i)
        contest1D(L[i, ], model=model, ddf=ddf, confint=TRUE, level=level)))
      rownames(tab) <- rownames(L)
      tab
    }))
  attr(means, "confidence_level") <- level
  attr(means, "ddf") <- ddf
  attr(means, "contrasts") <- Llist
  attr(means, "heading") <- "Least Squares Means table:\n"
  class(means) <- c("ls_means", "data.frame")
  means
}

##############################################
######## ls_means()
##############################################
#' LS-means Generic Function
#'
#' @param model a model object.
#' @param ... parsed on to methods.
#'
#' @export
#' @author Rune Haubo B. Christensen
#' @seealso \code{\link{ls_means.lmerModLmerTest}}
#' @keywords internal
ls_means <- function(model, ...) UseMethod("ls_means")

##############################################
######## difflsmeans()
##############################################
#' @rdname ls_means
#' @export
#' @seealso \code{\link{difflsmeans.lmerModLmerTest}}
#' @keywords internal
difflsmeans <- function(model, ...) UseMethod("difflsmeans")

##############################################
######## lsmeansLT()
##############################################
#' @rdname ls_means
#' @export
#' @seealso \code{\link{lsmeansLT.lmerModLmerTest}}
#' @keywords internal
lsmeansLT <- function(model, ...) UseMethod("lsmeansLT")

##############################################
######## lsmeansLT.lmerModLmerTest()
##############################################
#' @rdname ls_means.lmerModLmerTest
#' @export
lsmeansLT.lmerModLmerTest <- ls_means.lmerModLmerTest

##############################################
######## difflsmeans.lmerModLmerTest()
##############################################
#' @rdname ls_means.lmerModLmerTest
#' @export
difflsmeans.lmerModLmerTest <- function(model, which=NULL, level=0.95,
                        ddf=c("Satterthwaite", "Kenward-Roger"), ...) {
  ls_means(model, which=which, level=level, ddf=ddf, pairwise = TRUE)
}


##############################################
######## lsmeans_contrasts()
##############################################
lsmeans_contrasts <- function(model, which=NULL) {
  stopifnot(inherits(model, "lmerModLmerTest"))
  factor_terms <- attr(terms(model), "term.labels")[!numeric_terms(model)]
  if(is.null(which)) which <- factor_terms
  stopifnot(is.character(which), all(which %in% factor_terms))
  which <- setNames(as.list(which), which)

  # Get minimal 'unique rows' design matrix:
  grid <- get_min_data(model)
  form <- formula(model)[-2]
  if(inherits(model, "lmerMod")) form <- nobars(form)
  Contr <- attr(model.matrix(model), "contrasts")
  uX <- model.matrix(form, data=grid, contrasts.arg=Contr)
  # Get utilities needed to compute the LS-means contrasts:
  var_names <- names(get_var_list(model))
  factor_mat <- attr(terms(model), "factors")
  Contrasts <- .getXlevels(terms(model), grid)
  Contrasts[] <- "contr.treatment"

  # Compute LS-means contrast:
  Llist <- lapply(which, function(term) {
    vars_in_term <- factor_mat[var_names, term] == 1
    Lt <- model.matrix(formula(paste0("~ 0 + ", term)), data=grid,
                       contrasts.arg=Contrasts[vars_in_term])
    wts <- 1/colSums(Lt)
    # Lt * c(Lt %*% wts)
    # L <- diag(wts) %*% t(Lt)
    L <- t(sweep(Lt, 2, wts, "*"))
    L %*% uX
  })
  Llist
}


##############################################
######## print.ls_means
##############################################
#' @importFrom stats printCoefmat
print.ls_means <- function(x, digits = max(getOption("digits") - 2L, 3L),
                          signif.stars = getOption("show.signif.stars"),
                          ...) {
  if(!is.null(heading <- attr(x, "heading")))
    cat(heading, sep = "\n")
  if(nrow(x) > 0) {
    dig.df <- 1
    x[, "df"] <- round(x[, "df"], dig.df)
  }
  printCoefmat(x, digits=digits, signif.stars = signif.stars,
               has.Pvalue = TRUE, cs.ind=c(1:2, 5:6), tst.ind=4)
  if(!is.null(ci_level <- attr(x, "confidence_level")))
    cat(paste0("\n  Confidence level: ", format(100*ci_level, digits=2), "%\n"))
  if(!is.null(ddf <- attr(x, "ddf")))
    cat("  Degrees of freedom method:", ddf, "\n")
  invisible(x)
}


##############################################
######## show_contrasts
##############################################
#' Show LS-means Contrasts
#'
#' Extracts the contrasts which defines the LS-mean contrasts.
#'
#' @param object an \code{ls_means} object.
#' @param fractions display contrasts as fractions rather than decimal numbers?
#' @param names include row and column names of the contrasts matrices?
#'
#' @return a list of contrast matrices; one matrix for each model term.
#' @export
#' @author Rune Haubo B. Christensen
#' @importFrom MASS fractions
#' @seealso \code{\link[=ls_means.lmerModLmerTest]{ls_means}} for computation of LS-means.
#'
#' @examples
#'
#' data("cake", package="lme4")
#' model <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
#' (lsm <- ls_means(model))
#' show_contrasts(lsm)
#'
show_contrasts <- function(object, fractions=FALSE, names=TRUE) {
  tests <- attr(object, "contrasts")
  # FIXME: Maybe this should be a generic with a method for anova objects?
  if(is.null(tests))
    stop("'object' does not have an 'contrasts' attribute")
  if(fractions) tests <- lapply(tests, MASS::fractions)
  if(names) tests else lapply(tests, unname)
}

