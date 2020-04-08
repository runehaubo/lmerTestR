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
# lsmeans.R - lsmeans methods for lmerTest::lmer model fits

# ------- Contents: --------
#
# --- Generics: ---
#
# ls_means
# difflsmeans
# lsmeansLT
#
# --- methods: ---
#
# ls_means.lmerModLmerTest
# difflsmeans.lmerModLmerTest
# lsmeansLT.lmerModLmerTest
# print.ls_means
# plot.ls_means
# as.data.frame.ls_means
#
# show_tests.ls_means
#
# --- other exported function: ---
#
# show_contrasts
#
# --- utility functions: ---
#
# lsmeans_contrasts
#

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
#' @seealso \code{\link[=show_tests.ls_means]{show_tests}} for display of the
#' underlying LS-means contrasts.
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
#' # Get LS-means contrasts:
#' show_tests(ls_means(model))
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
    means <- contest1D(model, rep(NA_real_, length(coef_nm)), ddf=ddf,
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
        contest1D(model, L[i, ], ddf=ddf, confint=TRUE, level=level)))
      rownames(tab) <- rownames(L)
      tab
    }))
  attr(means, "response") <- deparse(formula(model)[[2]])
  attr(means, "confidence_level") <- level
  attr(means, "ddf") <- ddf
  attr(means, "hypotheses") <- Llist
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
#' @export
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
######## show_tests.ls_means
##############################################
#' Show LS-means Hypothesis Tests and Contrasts
#'
#' Extracts the contrasts which defines the LS-mean hypothesis tests.
#'
#' @param object an \code{ls_means} object.
#' @param fractions display contrasts as fractions rather than decimal numbers?
#' @param names include row and column names of the contrasts matrices?
#' @param ... currently not used.
#'
#' @return a list of contrast matrices; one matrix for each model term.
#' @export
#' @author Rune Haubo B. Christensen
#' @importFrom MASS fractions
#' @seealso \code{\link[=ls_means.lmerModLmerTest]{ls_means}} for computation of
#' LS-means and \code{\link[=show_tests.anova]{show_tests}} for \code{anova}
#' objects.
#'
#' @examples
#'
#' data("cake", package="lme4")
#' model <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
#'
#' # LS-means:
#' (lsm <- ls_means(model))
#'
#' # Contrasts for LS-means estimates and hypothesis tests:
#' show_tests(lsm)
#'
show_tests.ls_means <- function(object, fractions=FALSE, names=TRUE, ...)
  NextMethod() # use default method

##############################################
######## plot.ls_means
##############################################
#' Bar Plots of LS-Means
#'
#' Bar plots of LS-means using the \pkg{ggplot2} package.
#'
#' @param x an \code{\link{ls_means}} object.
#' @param y not used and ignored with a warning.
#' @param which optional character vector naming factors for which LS-means should
#' be plotted. If \code{NULL} (default) plots for all LS-means are generated.
#' @param mult if \code{TRUE} and there is more than one term for which to plot
#' LS-means the plots are organized in panels with \code{facet_wrap}.
#' @param ... currently not used.
#'
#' @return generates the desired plots and invisibly returns the plot objects.
#' @author Rune Haubo B. Christensen
#' @seealso \code{\link{ls_means.lmerModLmerTest}}
#' @export
#' @importFrom graphics plot
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar theme element_text
#' @importFrom ggplot2 scale_fill_manual xlab ylab facet_wrap rel
#' @keywords internal
#'
#' @examples
#'
#' # Fit example model with 2 factors:
#' data("cake", package="lme4")
#' cake$Temp <- factor(cake$temperature, ordered = FALSE)
#' model <- lmer(angle ~ recipe * Temp + (1|recipe:replicate), cake)
#'
#' # Extract LS-means:
#' (lsm <- ls_means(model))
#'
#' # Multi-frame plot of the LS-means
#' plot(lsm)
#'
#' # Compute list of 'single frame' plots:
#' res <- plot(lsm, mult=FALSE)
#'
#' # Display each plot separately:
#' plot(res[[1]])
#' plot(res[[2]])
#'
#' # Example with pairwise differences of LS-means:
#' (lsm <- ls_means(model, pairwise = TRUE))
#' plot(lsm, which="Temp")
#'
plot.ls_means <- function(x, y=NULL, which=NULL, mult=TRUE, ...) {
  Estimate <- col.bars <- lower <- term <- upper <- NULL # so that r cmd check can see them
  get_plot <- function(d, response="") { # basic plot function
    ggplot(d, aes(x=levels, y = Estimate, fill = col.bars)) +
      geom_bar(position = "dodge", stat = "identity") +
      geom_errorbar(aes(ymin = lower, ymax = upper), colour="black", width=.1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
            axis.title.y = element_text(size = rel(1.4)),
            axis.text = element_text(size = rel(1)),
            legend.text = element_text(size = rel(1)),
            legend.title = element_text(size = rel(1)))  +
      scale_fill_manual(
        values=c("NS" = "grey", "p-value < 0.01" = "orange",
                 "p-value < 0.05" = "yellow", "p-value < 0.001" = "red"),
        name="Significance") + ylab(response)
  }
  get_color_values <- function(x) {
    if(x<0.001) return("p-value < 0.001")
    if(x<0.01) return("p-value < 0.01")
    if(x<0.05) return("p-value < 0.05")
    return("NS")
  }
  # Check for and warn about deprecated arguments:
  dots <- list(...)
  ignored <- c("main", "cex")
  for(nm in ignored) if(any(pmatch(names(dots), nm, nomatch = 0)))
    warning(paste0("Argument '", nm, "' is deprecated and ignored."))
  if(any(pmatch(names(dots), "effs", nomatch = 0)))
    warning("Argument 'effs' is deprecated: use 'which' instead.")
  if(!is.null(y)) warning("Argument 'y' is defunct and ignored.")

  # Get data for plotting:
  plotdata <- as.data.frame(x, add_levels = TRUE)
  plotdata <- # Add significance information for colors:
    cbind(plotdata, col.bars=sapply(plotdata[, "Pr(>|t|)"], get_color_values))

  # Subset plotdata for terms
  if(!is.null(which)) {
    stopifnot(is.character(which), length(which) >= 1L,
              all(sapply(which, length) > 0L))
    term_names <- unique(as.character(plotdata[["term"]]))
    valid <- which %in% term_names
    if(!all(valid)) {
      warning(sprintf("The following terms are invalid and ignored: %s.",
                      paste(which[!valid], collapse = ", ")))
    }
    plotdata <- subset(plotdata, term %in% which[valid])
  }
  if(nrow(plotdata) == 0L) stop("No LS-means to plot.")

  # Generate plots:
  if(mult && length(unique(as.character(plotdata[["term"]]))) > 1L) {
    res <- get_plot(plotdata, response=attr(x, "response")) + xlab("") +
      facet_wrap( ~ term, scales="free")
    print(res)
  } else {
    plotdata <- split(plotdata, plotdata$term)
    res <- lapply(1:length(plotdata), function(i)
      get_plot(plotdata[[i]], response=attr(x, "response")) +
        xlab(names(plotdata)[i])
    )
    names(res) <- names(plotdata)
    for(obj in res) print(obj)
  }
  invisible(res)
}

##############################################
######## as.data.frame.ls_means
##############################################
#' Coerce \code{ls_means} Objects to \code{data.frame}s
#'
#' @param x an \code{\link{ls_means}} object.
#' @param add_levels add \code{term} and \code{levels} columns to returned
#' \code{data.frame}?
#' @param ... currently not used.
#'
#' @export
#' @author Rune Haubo B. Christensen
#' @seealso \code{\link{ls_means.lmerModLmerTest}}
#' @keywords internal
#' @examples
#'
#' # Fit example model:
#' data("cake", package="lme4")
#' cake$Temp <- factor(cake$temperature, ordered = FALSE)
#' model <- lmer(angle ~ recipe + Temp + (1|recipe:replicate), cake)
#'
#' # Extract LS-means:
#' head(lsm <- ls_means(model))
#'
#' # Coerce LS-means objects to data.frames:
#' head(as.data.frame(lsm))
#' head(as.data.frame(lsm, add_levels=FALSE))
#'
as.data.frame.ls_means <- function(x, ..., add_levels=TRUE) {
  # Function to compute levels of terms including interaction:terms
  get_levels <- function(term, levels) {
    fun <- function(term, levels) { # workhorse
      strng <- paste(paste0("^", unlist(strsplit(term, ":"))), collapse = "|")
      sapply(strsplit(levels, ":"), function(txt)
        paste(gsub(strng, "", txt), collapse = ":"))
    }
    if(all(grepl(" - ", levels))) # pairwise contrasts
      sapply(strsplit(levels, " - "), function(lev)
        paste(fun(term, lev), collapse = " - ")) else fun(term, levels)
  }
  if(!add_levels) return(structure(x, class="data.frame"))
  contrasts <- attr(x, "hypotheses")
  term_names <- names(contrasts)
  lsm_levels <- lapply(1:length(term_names), function(i)
    get_levels(term_names[i], rownames(contrasts[[i]]))
  )
  class(x) <- "data.frame"
  cbind(term = rep(term_names, sapply(lsm_levels, length)),
        levels=unlist(lsm_levels), x, stringsAsFactors=FALSE)
}

