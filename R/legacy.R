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
# legacy.R - support for lecacy 'merModLmerTest' objects.

# ------- Contents: --------
#
# --- Classes: ---
#
# merModLmerTest
#
# --- methods: ---
#
# anova.merModLmerTest
# summary.merModLmerTest
# ls_means.merModLmerTest
# lsmeansLT.merModLmerTest
# difflsmeans.merModLmerTest
# drop1.merModLmerTest
#

##############################################
######## merModLmerTest class
##############################################
#' Legacy lmerTest representation of Linear Mixed-Effects Models
#'
#' The \code{merModLmerTest} class extends \code{lmerMod} (which extends
#' \code{merMod}) from the \pkg{lme4}-package.
#'
#' @export
#' @keywords internal
#' @author Rune Haubo B. Christensen
#' @importClassesFrom lme4 lmerMod
merModLmerTest <- setClass("merModLmerTest", contains = c("lmerMod"))

##############################################
######## anova method for merModLmerTest
##############################################
#' Methods for Legacy lmerTest Objects
#'
#' Methods are defined for legacy lmerTest objects of class
#' \code{merModLmerTest} generated with \pkg{lmerTest} version \code{< 3.0-0}.
#' These methods are defined by interfacing code for \code{lmerModLmerTest}
#' methods and therefore behaves like these methods do (which may differ from
#' the behavior of \pkg{lmerTest} version \code{< 3.0-0}.)
#'
#' @inheritParams anova.lmerModLmerTest
#' @param ... for the anova method optionally additional models; for other
#' methods see the corresponding \code{lmerModLmerTest} methods for details.
#' @rdname legacy
#' @aliases legacy
#' @keywords internal
#' @author Rune Haubo B. Christensen
#' @export
#' @examples
#' # Load model fits fm1 and fm2 generated with lmerTest version 2.3-37:
#' load(system.file("testdata","legacy_fits.RData", package="lmerTest"))
#'
#' # Apply some methods defined by lmerTest:
#' anova(fm1)
#' summary(fm1)
#' contest(fm1, c(0, 1))
#' contest(fm1, c(0, 1), joint=FALSE)
#' drop1(fm1)
#' ranova(fm1)
#'
#' # lme4-methods also work:
#' fixef(fm1)
#'
#' # Ditto for second model fit:
#' anova(fm2)
#' summary(fm2)
#' ls_means(fm2)
#' difflsmeans(fm2)
anova.merModLmerTest <- function(object, ..., type = c("III", "II", "I", "3", "2", "1"),
                                 ddf = c("Satterthwaite", "Kenward-Roger", "lme4")) {
  class(object) <- "lmerMod"
  dots <- list(...)
  models <- if (length(dots))
    sapply(dots, is, "merModLmerTest") | sapply(dots, is, "lmerModLmerTest") |
    sapply(dots, is, "merMod") | sapply(dots, is, "lm")
  else logical(0)
  if(any(models)) return(NextMethod())
  df <- match.arg(ddf)
  if (df == "lme4")
    return(anova(object, ...))

  object <- as_lmerModLmerTest(object)
  anova(object, ..., type=type, ddf=ddf)
}

##############################################
######## summary method for merModLmerTest
##############################################
#' @rdname legacy
#' @export
summary.merModLmerTest <- function(object, ...,
                                   ddf=c("Satterthwaite", "Kenward-Roger", "lme4")) {
  class(object) <- "lmerMod"
  object <- as_lmerModLmerTest(object)
  summary.lmerModLmerTest(object=object, ..., ddf=ddf)
}

##############################################
######## ls_means method for merModLmerTest
##############################################
#' @rdname legacy
#' @inheritParams ls_means.lmerModLmerTest
#' @export
ls_means.merModLmerTest <- function(model, which=NULL, level=0.95,
                                    ddf=c("Satterthwaite", "Kenward-Roger"),
                                    pairwise=FALSE, ...) {
  class(model) <- "lmerMod"
  model <- as_lmerModLmerTest(model)
  ls_means(model=model, which=which, level=level, ddf=ddf, pairwise=pairwise)
}

##############################################
######## lsmeansLT method for merModLmerTest
##############################################
#' @rdname legacy
#' @export
lsmeansLT.merModLmerTest <- ls_means.merModLmerTest

##############################################
######## difflsmeans method for merModLmerTest
##############################################
#' @rdname legacy
#' @export
difflsmeans.merModLmerTest <- function(model, which=NULL, level=0.95,
                                       ddf=c("Satterthwaite", "Kenward-Roger"), ...) {
  ls_means(model, which=which, level=level, ddf=ddf, pairwise = TRUE)
}

##############################################
######## drop1 method for merModLmerTest
##############################################
#' @rdname legacy
#' @inheritParams drop1.lmerModLmerTest
#' @export
drop1.merModLmerTest <- function(object, scope, ddf=c("Satterthwaite", "Kenward-Roger", "lme4"),
                                  force_get_contrasts=FALSE, ...) {
  class(object) <- "lmerMod"
  object <- as_lmerModLmerTest(object)
  drop1(object=object, scope=scope, ddf=ddf, force_get_contrasts=FALSE, ...)
}

##############################################
######## step method for merModLmerTest
##############################################
#' @rdname legacy
#' @inheritParams step.lmerModLmerTest
#' @export
step.merModLmerTest <- function(object, ddf=c("Satterthwaite", "Kenward-Roger"),
                                alpha.random=0.1, alpha.fixed=0.05,
                                reduce.fixed=TRUE, reduce.random=TRUE,
                                keep, ...) {
  class(object) <- "lmerMod"
  object <- as_lmerModLmerTest(object)
  step(object, ddf=ddf, alpha.random=alpha.random, alpha.fixed=alpha.fixed,
       reduce.fixed=reduce.fixed, reduce.random=reduce.random,
       keep=keep, ...)
}
