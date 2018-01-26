# test_contestMD.R
library(lmerTestR)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

data("sleepstudy", package="lme4")

####################################
## Tests of contestMD
####################################

fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
           sleepstudy)
# Basic tests:
L <- diag(3L)
contestMD(L, fm)

# Tests of ddf arg:
contestMD(L, fm, ddf="Sat")
contestMD(L, fm, ddf="KR")
assertError(contestMD(L, fm, ddf="sat")) # Invalid ddf arg.

# Tests of simple 2-df test:
(ans <- contestMD(L[2:3, ], fm, ddf="Sat"))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 2L)
(ans <- contestMD(L[2:3, ], fm, ddf="KR"))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 2L)

# Tests of simple 1-df test:
(ans <- contestMD(L[3, , drop=FALSE], fm, ddf="Sat"))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 1L)
(ans <- contestMD(L[3, , drop=FALSE], fm, ddf="KR"))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 1L)

# Test of vector input:
(ans <- contestMD(L[3, ], fm, ddf="Sat")) # OK since length(L[3, ]) == length(fixef(fm))
stopifnot(nrow(ans) == 1L,
          ans$NumDF == 1L)
assertError(contestMD(c(1, 0), fm)) # L is too short
assertError(contestMD(c(1, 0, 1, 1), fm)) # L is too long

# Test of list input:
assertError(contestMD(list(L[3, , drop=FALSE]), fm, ddf="Sat")) # Need L to be a matrix

# zero-row L's are allowed (if ncol(L) is correct):
ans1 <- contestMD(L[0, , drop=FALSE], fm, ddf="Sat")
ans2 <- contestMD(L[0, , drop=FALSE], fm, ddf="KR")
stopifnot(nrow(ans1) == 0L,
          nrow(ans2) == 0L)

# Test wrong ncol(L):
assertError(contestMD(L[2:3, 2:3], fm)) # need ncol(L) == length(fixef(fm))

# row-rank deficient L are allowed:
L <- rbind(c(1, 0, 1),
           c(0, 1, 0),
           c(1, -1, 1))
ans <- contestMD(L, fm)
ans_KR <- contestMD(L, fm, ddf="KR")
stopifnot(nrow(L) == 3L,
          qr(L)$rank == 2,
          ans$NumDF == 2,
          ans_KR$NumDF == 2)

# Test of 0-length beta
fm1 <- lmer(Reaction ~ 0 + (1|Subject) + (0+Days|Subject),
            sleepstudy)
stopifnot(length(fixef(fm1)) == 0L)
L <- numeric(0L)
(ans <- contestMD(L, fm1))
stopifnot(nrow(ans) == 0L)
L <- matrix(numeric(0L), ncol=0L)
(ans <- contestMD(L, fm1))
stopifnot(nrow(ans) == 0L)

# Test using model of wrong class:
fm2 <- lme4::lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
                  sleepstudy)
assertError(contestMD(L[2:3, ], fm2)) # fm2 is not of class "lmerModLmerTest"

