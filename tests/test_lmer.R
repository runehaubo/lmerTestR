# test_lmer.R

stopifnot(!"lmerTest" %in% .packages()) # ensure that lmerTest is NOT attached
data("sleepstudy", package="lme4")
f <- function(form, data) lmerTest::lmer(form, data=data)
form <- "Reaction ~ Days + (Days|Subject)"
fm <- f(form, data=sleepstudy)
anova(fm)
summary(fm)

# cf. GitHub issue #2:
test <- function() {
  tmp <- sleepstudy
  m <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  summary(m)
}
test()
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


#####################################################################
# Check that lme4::lmer and lmerTest::lmer have the same arguments
lmer_args <- formals(lme4::lmer)
lmerTest_args <- formals(lmerTest::lmer)
stopifnot(
  all.equal(names(lmer_args), names(lmerTest_args)),
  all.equal(lmer_args, lmerTest_args)
)

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
stopifnot(isTRUE(all.equal(fmA, fmB)))
# Based on bug-report by Henrik Singmann, github issue #3

#####################################################################
# Test error message from as_lmerModLmerTest:
data("sleepstudy", package="lme4")
myfit <- function(formula, data) {
  lme4::lmer(formula = formula, data = data)
}
fm2 <- myfit(Reaction ~ Days + (Days|Subject), sleepstudy)
m <- assertError(as_lmerModLmerTest(fm2))
stopifnot(
  grepl("Unable to extract deviance function from model fit", m[[1]], fixed=TRUE)
)

#####################################################################
# Check that devFunOnly argument works:
data("sleepstudy", package="lme4")
fun <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, devFunOnly = TRUE)
stopifnot(is.function(fun) && names(formals(fun)[1]) == "theta")
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fun <- update(fm1, devFunOnly=TRUE)
stopifnot(is.function(fun) && names(formals(fun)[1]) == "theta")
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
# Expect warning when model fails to converge:
m <- assertWarning(
  fm <- lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
               (1|Consumer) + (1 |Consumer:Income), data=carrots)
)
# May not fail to converge on all platforms in all versions:
if(any(grepl("Model failed to converge", fm@optinfo$conv$lme$messages))) {
  messages <- unlist(lapply(m, "[[", "message"))
  stopifnot(
    any(grepl("Model failed to converge", messages))
  )
}

#####################################################################
# Update method

m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
m1 <- update(m, ~.-Days)
m2 <- lmer(Reaction ~ (Days | Subject), sleepstudy)

stopifnot(
  inherits(m, "lmerModLmerTest"),
  inherits(m1, "lmerModLmerTest"),
  inherits(m2, "lmerModLmerTest"),
  all.equal(m1, m2)
)

