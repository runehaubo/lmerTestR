# test_legacy.R
library(lmerTest)
TOL <- 1e-4
#####################################################################

# Read in data set
load(system.file("testdata", "legacy_fits.RData", package="lmerTest"))
# Generated with the following code using lmerTest version 2.0-37.9002
#
# library("lmerTest")
# packageVersion("lmerTest")
# fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
# (an1 <- anova(fm1))
# (sfm1 <- summary(fm1))
#
# fm2 <- lmer(Informed.liking ~ Product + Information + Gender +
#               (1|Product:Consumer) , data=ham)
# (an2 <- anova(fm2))
# (sfm2 <- summary(fm2))
#
# save(fm1, an1, sfm1, fm2, an2, sfm2,
#      file="~/GitHub/lmerTestR/package/inst/testdata/legacy_fits.RData")


#######################################
### Check that arguments for merModLmerTest and lmerModLmerTest methods match up:

stopifnot(
  isTRUE(all.equal(formals(lmerTest:::anova.merModLmerTest),
                   formals(lmerTest:::anova.lmerModLmerTest))),
  isTRUE(all.equal(formals(lmerTest:::summary.merModLmerTest),
                   formals(lmerTest:::summary.lmerModLmerTest))),
  isTRUE(all.equal(formals(lmerTest:::drop1.merModLmerTest),
                   formals(lmerTest:::drop1.lmerModLmerTest))),
  isTRUE(all.equal(formals(lmerTest:::step.merModLmerTest),
                   formals(lmerTest:::step.lmerModLmerTest))),
  isTRUE(all.equal(formals(lmerTest:::ls_means.merModLmerTest),
                   formals(lmerTest:::ls_means.lmerModLmerTest))),
  isTRUE(all.equal(formals(lmerTest:::difflsmeans.merModLmerTest),
                   formals(lmerTest:::difflsmeans.lmerModLmerTest))))


#######################################
## Tests for fm1:

an1new <- anova(fm1)
sfm1new <- summary(fm1)

stopifnot(
  isTRUE(all.equal(an1new, an1, check.attributes=FALSE, tol=TOL)),
  isTRUE(all.equal(coef(sfm1new), coef(sfm1), tol=TOL))
)

contest(fm1, c(0, 1))
contest(fm1, c(0, 1), joint=FALSE)
drop1(fm1)
ranova(fm1)
step(fm1)

fm1new <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
               control=lmerControl(optimizer="bobyqa"))
stopifnot(
  isTRUE(all.equal(drop1(fm1), drop1(fm1new), tol=TOL)),
  isTRUE(all.equal(ranova(fm1), ranova(fm1new), tol=TOL)),
  isTRUE(all.equal(contest(fm1, c(0, 1)), contest(fm1new, c(0, 1)), tol=TOL)),
  isTRUE(all.equal(contest(fm1, c(0, 1), joint=FALSE),
                   contest(fm1new, c(0, 1), joint=FALSE), tol=TOL))
)

# Test that lme4 methods work:
coef(fm1)
fixef(fm1)
resid(fm1)

#######################################
## Tests for fm2:
an2new <- anova(fm2)
sfm2new <- summary(fm2)

stopifnot(
  isTRUE(all.equal(an2new, an2, check.attributes=FALSE, tol=TOL)),
  isTRUE(all.equal(coef(sfm2new), coef(sfm2), tol=TOL))
)

drop1(fm2)
ranova(fm2)
ls_means(fm2)
difflsmeans(fm2)
nbeta <- length(fixef(fm2))
L <- diag(nbeta)
L[1:4, ] <- 0
contest(fm2, L)
contest(fm2, diag(nbeta), joint=FALSE)
step(fm2)

fm2new <- lmer(Informed.liking ~ Product + Information + Gender +
                 (1|Product:Consumer), data=ham)
stopifnot(
  isTRUE(all.equal(drop1(fm2), drop1(fm2new), tol=TOL)),
  isTRUE(all.equal(ranova(fm2), ranova(fm2new), tol=TOL)),
  isTRUE(all.equal(ls_means(fm2), ls_means(fm2new), tol=TOL)),
  isTRUE(all.equal(difflsmeans(fm2), difflsmeans(fm2new), tol=TOL))
)

# Test that lme4 methods work:
coef(fm2)
fixef(fm2)
resid(fm2)

