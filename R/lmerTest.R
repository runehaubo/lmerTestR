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
# lmerTest.R - package documentation page

#' lmerTest: Tests in Linear Mixed Effects Models
#'
#' The \pkg{lmerTest} package provides p-values in type I, II or III
#' \code{anova} and \code{summary}
#' tables for linear mixed models (\code{\link{lmer}} model fits cf. \pkg{lme4})
#' via Satterthwaite's degrees of freedom method; a Kenward-Roger method is also
#' available via the \pkg{pbkrtest} package.
#' Model selection and assessment methods include \code{\link{step}},
#' \code{\link{drop1}}, anova-like tables for random effects (\code{\link{ranova}}),
#' least-square means (LS-means; \code{\link{ls_means}})
#' and tests of linear contrasts of fixed effects (\code{\link{contest}}).
#'
#'
#' @section Key Functions and Methods:
#'
#' \describe{
#' \item{lmer}{overloads \code{lme4::lmer} and produced an object of class
#' \code{lmerModLmerTest} which inherits from \code{lmerMod}. In addition to
#' computing the model (using \code{lme4::lmer}), \code{lmerTest::lmer}
#' computes a couple of components needed for the evaluation of Satterthwaite's
#' denominator degrees of freedom.}
#' \item{anova}{anova method for \code{\link{lmer}} model fits produces
#' type I, II, and III anova tables for fixed-effect terms with
#' Satterthwaite and Kenward-Roger methods for denominator degrees of freedom
#' for F-tests.}
#' \item{summary}{summary method for \code{\link{lmer}} model fits adds
#' denominator degrees of freedom and p-values to the coefficient table.}
#' \item{ranova}{anova-like table of random effects via likelihood ratio tests
#' with methods for both \code{lmerMod} and \code{lmerModLmerTest} objects.
#' \code{ranova} can either test reduction of random-effect terms to simpler
#' structures or it can test removal of entire random-effect terms.}
#' \item{drop1}{F-tests of fixed-effect terms using Satterthwaite or
#' Kenward-Roger methods for denominator degrees of freedom. These 'single term
#' deletion' tables are useful for model selection and tests of marginal terms.
#' Compared to the likelihood ratio tests of \code{lme4::drop1} the F-tests and
#' p-values of \code{lmerTest::drop1} are more accurate and considerably faster
#' since no additional model fitting is required.}
#' \item{contest}{tests of contrasts, i.e. tests of linear functions of the
#' fixed-effect coefficients. A user-friendly interface for tests of contrasts
#' with outputs either as a summary-like table of t-tests or an anova-like table
#' of F-tests (or a list of either). Contrasts can optionally be tested for
#' estimability. Contrasts are allowed to be rank-deficient as the rank is
#' automatically detected and appropriate adjustments made. Methods for
#' \code{lmerModLmerTest} as well as \code{lmerMod} objects -- the latter avoids
#' the Satterthwaite specific computations when the Kenward-Roger method is used.}
#' \item{show_test}{a function which operates on anova tables and LS-means tables
#' makes it possible to see exactly which
#' functions of the coefficients are being tested. This is helpful when
#' differences between type I, II and III anova tables are being considered and
#' discussed.}
#' \item{ls_means}{computes the so-called least-squares means (classical Yates
#' contrasts) as well as pairwise differences of these.}
#' \item{step}{performs automatic backward model selection of fixed and random
#' parts of the linear mixed model.}
#' \item{as_lmerModLmerTest}{an explicit coerce function from class
#' \code{lmerMod} to \code{lmerModLmerTest}.}
#' }
#'
#' @section Details:
#' The computational approach is to let \code{lmerTest::lmer} compute the
#' Hessian and derivatives needed for evaluation of degrees of freedom and
#' t- and F-tests and to store these in the model object. The
#' Hessian and derivaties are therefore computed only once per model fit
#' and reused with each call to \code{anova}, \code{summary}, etc. Evaluation of
#' t and F-tests does not involve model re-fitting.
#'
#' \code{lmerTest::lmer} roughly amounts to calling \code{lme4::lmer} followed by
#' \code{lmerTest::as_lmerModLmerTest}, so for computationally intensive model
#' fits it can make sense to use \code{lme4::lmer} rather than \code{lmerTest:lmer}
#' if computational time is an issue and summary tables and anova tables will
#' not be needed.
#'
#' @author Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
#'
#' @references
#'
#'   Alexandra Kuznetsova, Per B. Brockhoff and Rune H. B. Christensen (2017)
#'   lmerTest Package: Tests in Linear Mixed Effects Models.
#'   \emph{Journal of Statistical Software}, 82(13), 1--26. doi:10.18637/jss.v082.i13
#'
#'
#' @docType package
#' @name lmerTest-package
#' @aliases lmerTest
#'
#' @examples
#'
#' ## load lmerTest package
#' library(lmerTest)
#'
#' ## Fit linear mixed model to the ham data:
#' fm <- lmer(Informed.liking ~ Gender + Information * Product + (1 | Consumer) +
#'              (1 | Consumer:Product), data=ham)
#'
#' ## Summary including coefficient table with p-values for t-statistics using
#' ## Satterthwaite's method for denominator degrees of freedom:
#' summary(fm)
#'
#' ## Type III anova table with p-values for F-tests based on Satterthwaite's
#' ## method:
#' (aov <- anova(fm))
#'
#' ## Inspect the contrast matrix for the Type III test of Product:
#' show_tests(aov, fractions = TRUE)$Product
#'
#' ## Choose type II anova table with Kenward-Roger method for the F-test:
#' \dontrun{
#' if(requireNamespace("pbkrtest", quietly = TRUE))
#'   anova(fm, type=2, ddf="Kenward-Roger")
#' }
#'
#' ## Anova-like table of random-effect terms using likelihood ratio tests:
#' ranova(fm)
#'
#' ## F-tests of 'single term deletions' for all marginal terms:
#' drop1(fm)
#'
#' ## Least-Square means and pairwise differences:
#' (lsm <- ls_means(fm))
#' ls_means(fm, which = "Product", pairwise = TRUE)
#'
#' ## ls_means also have plot and as.data.frame methods:
#' \dontrun{
#' plot(lsm, which=c("Product", "Information"))
#' as.data.frame(lsm)
#' ## Inspect the LS-means contrasts:
#' show_tests(lsm, fractions=TRUE)$Product
#' }
#'
#' ## Contrast test (contest) using a custom contrast:
#' ## Here we make the 2-df joint test of the main effects of Gender and Information
#' (L <- diag(length(fixef(fm)))[2:3, ])
#' contest(fm, L = L)
#'
#' ## backward elimination of non-significant effects:
#' step_result <- step(fm)
#'
#' ## Elimination tables for random- and fixed-effect terms:
#' step_result
#'
#' # Extract the model that step found:
#' final_model <- get_model(step_result)
#'
NULL
