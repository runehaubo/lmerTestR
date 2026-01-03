# test_re_covar_structure.R
library(lmerTest)

# library(devtools)
# r2path <- "~/GitHub/lmerTestR/lmerTest"
# load_all(r2path)


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
TOL <- 1e-4

####################################
## Test that lmerMod objects can be coerced to lmerModLmerTest and 
## that the deviance function can be evaluated with expected results
####################################


####################################
## Basic tests of new (lme4 >= 2.0.0) covariance structure
##  - simple sleepstudy data (random coefficient)
####################################

## From the examples of ?`Covariance-class`:
## Unstructured
fm1.us <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
## Diagional
fm1.diag <- lmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
fm1.diag.hom <- lmer(Reaction ~ Days + diag(Days | Subject, hom = TRUE), 
                     sleepstudy)
## Compound symmetry
fm1.cs <- lmer(Reaction ~ Days + cs(1 + Days | Subject), sleepstudy)
fm1.cs.hom <- lmer(Reaction ~ Days + cs(1 + Days | Subject, hom = TRUE), 
                   sleepstudy)
## Auto-regressive order 1
sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm1.ar1 <- lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), 
                sleepstudy, REML = TRUE)

## Also adding a double-vertical-bar model (though not from 
## 'Covariance-class' examples):
fm1.bv <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)


ltmodels <- namedList(fm1.us,
                      fm1.diag, 
                      fm1.diag.hom,
                      fm1.cs, 
                      fm1.cs.hom, 
                      fm1.ar1,
                      fm1.bv)
## Run various methods on all models:
for(model in ltmodels) {
  # model <- ltmodels[[7]]
  print(model)
  print(summary(model))
  # model
  # summary(model)
  L <- diag(c(0, rep(1, length(fixef(model)) -1)))
  contest(model, L, joint = TRUE)
  contest(model, L[2, ], joint = FALSE)
  (an <- anova(model)) ## ddf is Satterthwaite
  anova(model, ddf = "Ken")
  anova(model, ddf="lme4")
  anova(model, type="I")
  anova(model, type="II")
  anova(model, type="III")
  show_tests(an)
  drop1(model)
  ranova(model)
  ranova(model, reduce.terms = FALSE)
  step(model)
  (lsm <- ls_means(model))
  show_tests(lsm)
  (dlsm <- difflsmeans(model))
  show_tests(dlsm)
}


## FIXME:
fm1.bv <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
ranova(fm1.bv, reduce.terms = FALSE)
## It seems to me that here ranova should drop the whole term. Much like it 
## correctly does with the diag-model:
ranova(fm1.diag, reduce.terms = FALSE)
## Even here the printing is off as 'diag' appears to be dropped.

####################################
## Basic tests of new (lme4 >= 2.0.0) covariance structure
##  - Unbalanced categorical dataset with multiple RE terms
####################################




