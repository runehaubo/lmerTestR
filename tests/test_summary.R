# test_summary.R

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

library(lmerTest)
data("sleepstudy", package="lme4")
data("cake", package="lme4")

# Fit basic model and compute summary:
fm <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
(sfm <- summary(fm))

## Test class:
stopifnot(all(
  class(sfm) == c("summary.lmerModLmerTest", "summary.merMod"),
  all(c("df", "Pr(>|t|)") %in% colnames(coef(sfm)))
))
stopifnot(class(summary(fm, ddf="lme4")) == "summary.merMod")

## Test coefficient table names:
mat <- coef(summary(fm))
stopifnot(all( # colnames
  colnames(mat) == c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")
))
stopifnot(all( # rownames
  names(fixef(fm)) == rownames(mat)
))

## Test pass of 'correlation' argument to lme4:::summary.merMod:
x <- capture.output(summary(fm))
x_nocor <- capture.output(summary(fm, correlation=FALSE))
txt <- "Correlation of Fixed Effects:"
stopifnot(
  any(grep(txt, x)),
  !any(grepl(txt, x_nocor))
)

# Test warning with unrecognized arguments (caught by lme4:::summary.merMod):
assertWarning(summary(fm, false_arg=FALSE))

## Test pass of extra arguments to lme4:::print.summary.merMod:
x <- capture.output(print(summary(fm), signif.stars=TRUE))
x_nocor <- capture.output(print(summary(fm), signif.stars=FALSE))
txt <- "Signif. codes:"
stopifnot(
  any(grep(txt, x)),
  !any(grepl(txt, x_nocor))
)

####### ddf argument:
(an1 <- summary(fm)) # Also testing print method.
(an2 <- summary(fm, ddf="Satterthwaite"))
stopifnot(isTRUE(
  all.equal(an1, an2)
))
(an3 <- summary(fm, ddf="Sat")) ## Abbreviated argument
stopifnot(isTRUE(
  all.equal(an1, an3)
))
(summary(fm, ddf="Kenward-Roger"))
(summary(fm, ddf="lme4"))
assertError(summary(fm, ddf="KR")) ## Error on incorrect arg.

## lme4 method:
an1 <- summary(fm, ddf="lme4")
an2 <- summary(as(fm, "lmerMod"))
stopifnot(isTRUE(
  all.equal(an1, an2)
))


# Test printed output
# - Satterthwaite
x <- capture.output(sfm) # equal to output of 'print(sfm)'
txt <- c("lmerModLmerTest", "t-tests use Satterthwaite's method",
         "df", "t value", "Pr(>|t|)")
stopifnot(all(
  sapply(txt, function(text) any(grepl(text, x)))
))

# Test printed output
# - KR
(sfm <- summary(fm, ddf="Kenward-Roger"))
x <- capture.output(sfm)
txt <- c("lmerModLmerTest", "t-tests use Kenward-Roger's method",
         "df", "t value", "Pr(>|t|)")
stopifnot(all(
  sapply(txt, function(text) any(grepl(text, x)))
))


####################################
## Test 'boundary' fixef structures:
####################################

# Example with no fixef:
m <- lmer(Reaction ~ -1 + (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 0 + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 0L)
stopifnot(
  nrow(coef(summary(m))) == 0L,
  nrow(coef(summary(m, ddf="Kenward-Roger"))) == 0L,
  nrow(coef(summary(m, ddf="lme4"))) == 0L
)

# Example with intercept only:
m <- lmer(Reaction ~ (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 1 + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 1L,
          names(fixef(m)) == "(Intercept)")
stopifnot(
  nrow(coef(summary(m))) == 1L,
  nrow(coef(summary(m, ddf="Kenward-Roger"))) == 1L,
  nrow(coef(summary(m, ddf="lme4"))) == 1L
)

# Example with >1 fixef without intercept:
m <- lmer(Reaction ~ Days - 1 + I(Days^2) + (Days | Subject), sleepstudy)
stopifnot(length(fixef(m)) == 2L,
          names(fixef(m)) == c("Days", "I(Days^2)"))
stopifnot(
  nrow(coef(summary(m))) == 2L,
  nrow(coef(summary(m, ddf="Kenward-Roger"))) == 2L,
  nrow(coef(summary(m, ddf="lme4"))) == 2L
)

