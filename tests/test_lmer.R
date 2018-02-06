# test_lmer.R
library(lmerTestR)

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
bm <- lmerTestR:::as_lmerModLmerTest(m)
stopifnot(
  inherits(bm, "lmerModLmerTest"),
  !inherits(m, "lmerModLmerTest"),
  inherits(bm, "lmerMod"),
  all(c("vcov_varpar", "Jac_list", "vcov_beta", "sigma") %in% slotNames(bm))
)

#####################################################################
# Expect warning when model fails to converge:
data("carrots", package="lmerTest")
m <- assertWarning(
  lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
         (1|Consumer) + (1 |Consumer:Income), data=carrots)
)
messages <- unlist(lapply(m, "[[", "message"))
stopifnot(
  any(grepl("Model failed to converge with 1 negative eigenvalue", messages))
)
