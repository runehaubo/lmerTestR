# test_lmer.R
library(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

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
  lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
         (1|Consumer) + (1 |Consumer:Income), data=carrots)
)
messages <- unlist(lapply(m, "[[", "message"))
stopifnot(
  any(grepl("Model failed to converge", messages))
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
  all.equal(m1, m2)
)

