# test_drop1.R

library(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

TOL <- 1e-4
# Kenward-Roger only available with pbkrtest and only then validated in R >= 3.3.3
# (faulty results for R < 3.3.3 may be due to unstated dependencies in pbkrtest)
has_pbkrtest <- requireNamespace("pbkrtest", quietly = TRUE) && getRversion() >= "3.3.3"

data("sleepstudy", package="lme4")

######### Basic usage

data("cake", package="lme4")
cake2 <- cake
cake2$temperature <- factor(cake2$temperature, ordered = FALSE)
fm <- lmer(angle ~ recipe + temperature + (1|recipe:replicate), cake2)
(an1 <- drop1(fm))
(an2 <- drop1(fm, force_get_contrasts = TRUE))
drop1(fm, ddf="lme4", test="Chi")
if(has_pbkrtest)
  drop1(fm, ddf="Kenward-Roger")

tests1 <- show_tests(an1)
tests2 <- show_tests(an2)

stopifnot(
  # Tests are the same:
  isTRUE(all.equal(an1, an2, check.attributes = FALSE, tolerance=TOL)),
  # But contrast matrices are not:
  all(!mapply(function(x, y) isTRUE(all.equal(x, y)), tests1, tests2))
)

fm <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake2)
drop1(fm)
drop1(fm, ddf="lme4")
if(has_pbkrtest)
  drop1(fm, ddf="Kenward-Roger")

# Incorrect arguments:
assertError(drop1(fm, scope="recipe")) # Correct Error
assertError(drop1(fm, scope=3)) # Correct Error
assertError(drop1(fm, scope=list("recipe"))) # Correct Error

# Polynomial terms:

fm <- lmer(Reaction ~ 0 + (Days|Subject), sleepstudy)
(an0 <- drop1(fm)) # No fixef!
fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
(an1 <- drop1(fm))
fm <- lmer(Reaction ~ Days + I(Days^2) + (Days|Subject), sleepstudy)
(an2 <- (drop1(fm)))
fm <- lmer(Reaction ~ poly(Days, 2) + (Days|Subject), sleepstudy)
(an3 <- drop1(fm))
stopifnot(
  nrow(an0) == 0L,
  nrow(an1) == 1L,
  nrow(an2) == 2L,
  nrow(an3) == 1L
)

# Consider a rank-deficient design matrix:
fm <- lmer(angle ~ recipe + temp + temperature + (1|recipe:replicate), cake)
# Here temp accounts for the linear effect of temperature, and
# temperature is an (ordered) factor that accounts for the remaining
# variation between temperatures (4 df).
(an4 <- drop1(fm))
# While temperature is in the model, we cannot test the effect of dropping
# temp. After removing temperature we can test the effect of dropping temp:
(an5 <- drop1(update(fm, ~.-temperature)))

stopifnot(
  nrow(an4) == 3,
  rownames(an4)[2] == "temp",
  all(is.na(an4[2, ])),
  all(!is.na(an4[-2, ])),
  all(rownames(an5) == c("recipe", "temp"))
)



