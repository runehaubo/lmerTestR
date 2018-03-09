# test_contest1D.R
library(lmerTest)

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
contest1D(fm, L)
contest1D(fm, L, confint = TRUE)
contest1D(fm, L, confint = TRUE, level=0.99)
contest1D(fm, L, ddf="Kenward-Roger")

# Test too long L
assertError(contest1D(fm, c(0, 1, 1, 1)))

# Test too short L
assertError(contest1D(fm, c(0, 1)))

# Test matrix L
contest1D(fm, matrix(L, nrow=1))
contest1D(fm, matrix(L, ncol=1))
assertError(contest1D(fm, matrix(c(0, 1), ncol=1)))
assertError(contest1D(fm, matrix(c(0, 1, 0, 0), nrow=1)))
L <- matrix(numeric(0L), ncol=3)
assertError(contest1D(fm, L)) # "empty" matrix
assertError(contest1D(fm, matrix(1, ncol=3, nrow=2)))

# Test list L
assertError(contest1D(fm, list(c(0, 1, 0))))

# Test equivalence to coef(summary(fm)):
Lmat <- diag(length(fixef(fm)))
(coef_mat <- lmerTest:::rbindall(lapply(1:ncol(Lmat), function(i)
  contest1D(fm, Lmat[i, ]))))
(coef_mat_KR <- lmerTest:::rbindall(lapply(1:ncol(Lmat), function(i)
  contest1D(fm, Lmat[i, ], ddf="Kenward-Roger"))))
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
(ans <- contest1D(fm1, numeric(0L), ddf="Kenward-Roger"))
stopifnot(nrow(ans) == 0L)

## Test rhs argument:
fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
contest1D(fm, L=cbind(0, 1))
contest1D(fm, L=cbind(0, 1), ddf="Kenward-Roger")
contest1D(fm, L=cbind(0, 1), rhs=10)
contest1D(fm, L=cbind(0, 1), ddf="Kenward-Roger", rhs=10)

contest1D(fm, L=c(0, 1), rhs = 10.467)

(ct1 <- contest1D(fm, L=cbind(c(0, 1)), rhs = 10))
(ct2 <- contestMD(fm, L=rbind(c(0, 1)), rhs = 10))
stopifnot(
  isTRUE(all.equal(ct1[, "t value"]^2, ct2[, "F value"]))
)

## Test 'lmerMod' method:
fm <- lme4::lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
                 sleepstudy)
# Basic tests:
L <- c(0, 1, 0)
contest1D(fm, L)
contest1D(fm, L, confint = TRUE)
contest1D(fm, L, confint = TRUE, level=0.99)
contest1D(fm, L, ddf="Kenward-Roger")

