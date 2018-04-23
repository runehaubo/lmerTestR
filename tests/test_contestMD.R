# test_contestMD.R
library(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

# Kenward-Roger only available with pbkrtest and only then validated in R >= 3.3.3
# (faulty results for R < 3.3.3 may be due to unstated dependencies in pbkrtest)
has_pbkrtest <- requireNamespace("pbkrtest", quietly = TRUE) && getRversion() >= "3.3.3"

data("sleepstudy", package="lme4")

####################################
## Tests of contestMD
####################################

fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
           sleepstudy)
# Basic tests:
L <- diag(3L)
contestMD(fm, L)

# Tests of ddf arg:
contestMD(fm, L, ddf="Sat")
if(has_pbkrtest)
  contestMD(fm, L, ddf="Kenward-Roger")
assertError(contestMD(fm, L, ddf="sat")) # Invalid ddf arg.

# Tests of simple 2-df test:
(ans <- contestMD(fm, L[2:3, ], ddf="Sat"))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 2L)
if(has_pbkrtest) {
  (ans <- contestMD(fm, L[2:3, ], ddf="Kenward-Roger"))
  stopifnot(nrow(ans) == 1L,
            ans$NumDF == 2L)
}

# Tests of simple 1-df test:
(ans <- contestMD(fm, L[3, , drop=FALSE], ddf="Sat"))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 1L)
if(has_pbkrtest) {
  (ans <- contestMD(fm, L[3, , drop=FALSE], ddf="Kenward-Roger"))
  stopifnot(nrow(ans) == 1L,
            ans$NumDF == 1L)
}

# Test of vector input:
(ans <- contestMD(fm, L[3, ], ddf="Sat")) # OK since length(L[3, ]) == length(fixef(fm))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 1L)
assertError(contestMD(fm, c(1, 0))) # L is too short
assertError(contestMD(fm, c(1, 0, 1, 1))) # L is too long

# Test of list input:
assertError(contestMD(fm, list(L[3, , drop=FALSE]), ddf="Sat")) # Need L to be a matrix

# zero-row L's are allowed (if ncol(L) is correct):
ans1 <- contestMD(fm, L[0, , drop=FALSE], ddf="Sat")
stopifnot(nrow(ans1) == 0L)
if(has_pbkrtest) {
  ans2 <- contestMD(fm, L[0, , drop=FALSE], ddf="Kenward-Roger")
  stopifnot(nrow(ans2) == 0L)
}

# Test wrong ncol(L):
assertError(contestMD(fm, L[2:3, 2:3])) # need ncol(L) == length(fixef(fm))

# row-rank deficient L are allowed:
L <- rbind(c(1, 0, 1),
           c(0, 1, 0),
           c(1, -1, 1))
ans <- contestMD(fm, L)
stopifnot(nrow(L) == 3L,
          qr(L)$rank == 2,
          ans$NumDF == 2)
if(has_pbkrtest) {
  ans_KR <- contestMD(fm, L, ddf="Kenward-Roger")
  stopifnot(ans_KR$NumDF == 2)
}

# Test of 0-length beta
fm1 <- lmer(Reaction ~ 0 + (1|Subject) + (0+Days|Subject),
            sleepstudy)
stopifnot(length(fixef(fm1)) == 0L)
L <- numeric(0L)
(ans <- contestMD(fm1, L))
stopifnot(nrow(ans) == 0L)
L <- matrix(numeric(0L), ncol=0L)
(ans <- contestMD(fm1, L))
stopifnot(nrow(ans) == 0L)


## rhs argument:
data("cake", package="lme4")
model <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
(L <- diag(length(fixef(model)))[2:3, ])
(an <- anova(model, type="marginal"))

ct <- contestMD(model, L, rhs = 0)
ct2 <- contestMD(model, L, rhs = c(2, 2))
stopifnot(
  isTRUE(all.equal(ct[1, ], an[1, ], check.attributes=FALSE, tolerance=1e-6)),
  ct[, "F value"] < ct2[, "F value"]
)

L2 <- rbind(L, L[1, ] + L[2, ]) # rank deficient!
contestMD(model, L2, rhs = c(0, 0, 0)) # no warning
assertWarning(contestMD(model, L2, rhs = c(2, 2, 2))) # warning since L2 is rank def.
if(has_pbkrtest)
  assertWarning(contestMD(model, L2, rhs = c(2, 2, 2), ddf="Kenward-Roger"))

fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
contestMD(fm, L=cbind(0, 1))
contestMD(fm, L=cbind(0, 1), rhs=10)
if(has_pbkrtest) {
  contestMD(fm, L=cbind(0, 1), ddf="Kenward-Roger")
  contestMD(fm, L=cbind(0, 1), ddf="Kenward-Roger", rhs=10)
}


## Test 'lmerMod' method:
fm <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
contestMD(fm, L=cbind(0, 1))
contestMD(fm, L=cbind(0, 1), rhs=10)
if(has_pbkrtest) {
  contestMD(fm, L=cbind(0, 1), ddf="Kenward-Roger")
  contestMD(fm, L=cbind(0, 1), ddf="Kenward-Roger", rhs=10)
}
