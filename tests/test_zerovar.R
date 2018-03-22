# test_zerovar.R

library(lmerTest)
data("sleepstudy", package="lme4")

# Baseline fit:
m0 <- lmer(Reaction ~ Days +  (Days | Subject), sleepstudy)
m0
(an0 <- anova(m0))

# Make a fit with a zero variance estimate:
n <- nrow(sleepstudy)
g <- factor(rep(1:2, c(n - 10, 10)))
m <- lmer(Reaction ~ Days +  (Days | Subject) + (1|g), sleepstudy)
m
(an <- anova(m))

# check that fit has a zero variance
vc <- as.data.frame(VarCorr(m))
stopifnot(isTRUE(
  all.equal(0, vc[vc$grp == "g", "sdcor"], tol=1e-4)
))
# The hessian/vcov is actually positive definite:
stopifnot(isTRUE(
  all(eigen(m@vcov_varpar, only.values = TRUE)$values > 0)
))

# Check that ANOVA tables are the same:
stopifnot(isTRUE(
  all.equal(an0[, 1:5], an[, 1:5], tolerance=1e-4)
))

stopifnot(isTRUE( # Equality of summary tables
  all.equal(coef(summary(m0)), coef(summary(m)), tolerance=1e-4)
))
stopifnot(isTRUE( # Equality of lme4-anova tables
  all.equal(anova(m0, ddf="lme4"), anova(m, ddf="lme4"), tolerance=1e-4)
))

