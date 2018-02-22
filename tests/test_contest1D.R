# test_contest1D.R
library(lmerTestR)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

data("sleepstudy", package="lme4")

####################################
## Tests of contest1D
####################################

fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
           sleepstudy)
# Basic tests:
L <- c(0, 1, 0)
contest1D(L, fm)
contest1D(L, fm, confint = TRUE)
contest1D(L, fm, confint = TRUE, level=0.99)
contest1D(L, fm, ddf="KR")

# Test too long L
assertError(contest1D(c(0, 1, 1, 1), fm))

# Test too short L
assertError(contest1D(c(0, 1), fm))

# Test matrix L
contest1D(matrix(L, nrow=1), fm)
contest1D(matrix(L, ncol=1), fm)
assertError(contest1D(matrix(c(0, 1), ncol=1), fm))
assertError(contest1D(matrix(c(0, 1, 0, 0), nrow=1), fm))
L <- matrix(numeric(0L), ncol=3)
assertError(contest1D(L, fm)) # "empty" matrix
assertError(contest1D(matrix(1, ncol=3, nrow=2), fm))

# Test list L
assertError(contest1D(list(c(0, 1, 0)), fm))

# Test equivalence to coef(summary(fm)):
Lmat <- diag(length(fixef(fm)))
(coef_mat <- lmerTestR:::rbindall(lapply(1:ncol(Lmat), function(i)
  contest1D(Lmat[i, ], fm))))
(coef_mat_KR <- lmerTestR:::rbindall(lapply(1:ncol(Lmat), function(i)
  contest1D(Lmat[i, ], fm, ddf="KR"))))
(coef_mat_lme4 <- coef(summary(fm, ddf="lme4")))
rownames(coef_mat_KR) <- rownames(coef_mat) <- rownames(coef_mat_lme4)
stopifnot(isTRUE(
  all.equal(as.data.frame(coef_mat_lme4),
            coef_mat[, c("Estimate", "Std. Error", "t value")])
))
stopifnot(isTRUE(
  all.equal(as.data.frame(coef_mat_lme4),
            coef_mat_KR[, c("Estimate", "Std. Error", "t value")], tol=1e-4)
))

# Test of 0-length beta
fm1 <- lmer(Reaction ~ 0 + (1|Subject) + (0+Days|Subject),
            sleepstudy)
stopifnot(length(fixef(fm1)) == 0L)
(ans <- contest1D(numeric(0L), fm1, ddf="KR"))
stopifnot(nrow(ans) == 0L)

# Test using model of wrong class:
fm2 <- lme4::lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
                  sleepstudy)
assertError(contest1D(c(0, 1, 0), fm2)) # fm2 is not of class "lmerModLmerTest"
