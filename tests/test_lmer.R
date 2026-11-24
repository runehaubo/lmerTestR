# test_lmer.R

if(isNamespaceLoaded("lmerTest")) unloadNamespace("lmerTest")
data("sleepstudy", package="lme4")
lmerTest::lmer(Reaction ~ 1 + (1|Subject), sleepstudy)

if(isNamespaceLoaded("lmerTest")) unloadNamespace("lmerTest")
f <- function(form, data) lmerTest::lmer(form, data=data)
form <- "Reaction ~ Days + (Days|Subject)"
fm <- f(form, data=sleepstudy)
anova(fm)
summary(fm)

# cf. GitHub issue #2:
if(isNamespaceLoaded("lmerTest")) unloadNamespace("lmerTest")
test <- function() {
  tmp <- sleepstudy
  m <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  summary(m)
}
test()

if(isNamespaceLoaded("lmerTest")) unloadNamespace("lmerTest")
if(isNamespaceLoaded("lme4")) unloadNamespace("lme4")
test <- function() {
  tmp <- sleepstudy
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  if(requireNamespace("lmerTest", quietly = TRUE)) {
    summary(lmerTest::as_lmerModLmerTest(m))
  }
}
test()

library(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

TOL <- 1e-4

#####################################################################
# Check that lme4::lmer and lmerTest::lmer have the same arguments

lmer_args <- formals(lme4::lmer)
names(lmer_args)
lmerTest_args <- formals(lmerTest::lmer)
seq_args <- seq_along(lmerTest_args)
if(packageVersion("lme4") > '1.1.21') {
  stopifnot(
    all.equal(names(lmer_args), names(lmerTest_args)),
    all.equal(lmer_args, lmerTest_args)
  )
} else { # Older versions of 'lme4' has a "..." argument:
  stopifnot(
    all.equal(names(lmer_args)[seq_args], names(lmerTest_args[seq_args])),
    all.equal(lmer_args[seq_args], lmerTest_args[seq_args])
  )
}

#####################################################################
# Test evaluation of update inside a function:
myupdate <- function(m, ...) {
  update(m, ...)
}

data("sleepstudy", package="lme4")
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
tmp <- sleepstudy
rm(sleepstudy)
fmA <- update(fm1, data = tmp) # works
fmB <- myupdate(fm1, data = tmp) # also works
# Same except for 'call':
fmB@call <- fmA@call
stopifnot(isTRUE(all.equal(fmA, fmB, tolerance=TOL)))
# Based on bug-report by Henrik Singmann, github issue #3

#####################################################################
# Test update when formula is a character vector:

form <- "Informed.liking ~ Product+Information+
            (1|Consumer) + (1|Product:Consumer) + (1|Information:Consumer)"
m <- lmer(form, data=ham)
class(m)
class(update(m, ~.- Product))
stopifnot(inherits(update(m, ~.- Product), "lmerModLmerTest"))

# In version < 3.0-1.9002 class(update(m, ~.- Product)) was "lmerMod"
#####################################################################
# Coerce fits whose original call cannot be evaluated:
data("sleepstudy", package="lme4")
myfit <- function(formula, data) {
  lme4::lmer(formula = formula, data = data)
}
fm2 <- myfit(Reaction ~ Days + (Days|Subject), sleepstudy)
wrapped <- as_lmerModLmerTest(fm2)
original_data <- sleepstudy
fm2 <- lme4::lmer(Reaction ~ Days + (Days|Subject), data=original_data)
expected <- as_lmerModLmerTest(fm2)
saved <- serialize(fm2, NULL)
rm(original_data)
restored <- unserialize(saved)
fitted_before <- serialize(fitted(restored), NULL)
converted <- as_lmerModLmerTest(restored)
stopifnot(
  inherits(converted, "lmerModLmerTest"),
  all.equal(coef(summary(wrapped)), coef(summary(expected)), tolerance=TOL),
  all.equal(coef(summary(converted)), coef(summary(expected)), tolerance=TOL),
  identical(serialize(fitted(restored), NULL), fitted_before)
)

#####################################################################
# Check that devFunOnly argument works:
data("sleepstudy", package="lme4")
fun <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, devFunOnly = TRUE)
stopifnot(is.function(fun)) # && names(formals(fun)[1]) == "theta")
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fun <- update(fm1, devFunOnly=TRUE)
stopifnot(is.function(fun)) # && names(formals(fun)[1]) == "theta")
# devFunOnly = FALSE:
notfun <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, devFunOnly = FALSE)
stopifnot(inherits(notfun, "lmerModLmerTest"))
# Partial matching:
notfun <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, devFun = FALSE)
stopifnot(inherits(notfun, "lmerModLmerTest"))

#####################################################################
# Use of as_lmerModLmerTest
data("sleepstudy", package="lme4")
m <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
bm <- lmerTest:::as_lmerModLmerTest(m)
stopifnot(
  inherits(bm, "lmerModLmerTest"),
  !inherits(m, "lmerModLmerTest"),
  inherits(bm, "lmerMod"),
  all(c("vcov_varpar", "Jac_list", "vcov_beta", "sigma") %in% slotNames(bm))
)

#####################################################################
# Update method

m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
m1 <- update(m, ~.-Days)
m2 <- lmer(Reaction ~ (Days | Subject), sleepstudy)

stopifnot(
  inherits(m, "lmerModLmerTest"),
  inherits(m1, "lmerModLmerTest"),
  inherits(m2, "lmerModLmerTest"),
  all.equal(m1, m2, tolerance=1e-6)
)

