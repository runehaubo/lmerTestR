# test_ranova.R

# Test functionality _before_ attaching lmerTest
stopifnot(!"lmerTest" %in% .packages()) # ensure that lmerTest is NOT attached
data("sleepstudy", package="lme4")
f <- function(form, data) lmerTest::lmer(form, data=data)
form <- "Reaction ~ Days + (Days|Subject)"
fm <- f(form, data=sleepstudy)
lmerTest::ranova(fm)
lmerTest::rand(fm)
lmerTest::step(fm)

library(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

TOL <- 1e-4
#####################################################################
data("sleepstudy", package="lme4")

# Test reduction of (Days | Subject) to (1 | Subject):
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
(an <- rand(fm1)) # 2 df test
(an <- ranova(fm1)) # 2 df test
step(fm1)
stopifnot(
  nrow(an) == 2L,
  an[2L, "Df"] == 2L
)

# This test can also be achieved with anova():
fm2 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
(stp <- step(fm2))
get_model(stp)
(ana <- anova(fm1, fm2, refit=FALSE))

stopifnot(
  all.equal(an[2L, "LRT"], ana[2L, "Chisq"], tolerance=TOL)
)

# Illustrate complete.test argument:
# Test removal of (Days | Subject):
(an <- ranova(fm1, reduce.terms = FALSE)) # 3 df test

# The likelihood ratio test statistic is in this case:
fm3 <- lm(Reaction ~ Days, sleepstudy)
LRT <- 2*c(logLik(fm1, REML=TRUE) - logLik(fm3, REML=TRUE)) # LRT
stopifnot(
  nrow(an) == 2L,
  an[2L, "Df"] == 3L,
  all.equal(an[2L, "LRT"], LRT, tolerance=TOL)
)

## _NULL_ model:
fm <- lmer(Reaction ~ -1 + (1|Subject), sleepstudy)
step(fm)
ranova(fm)
lm1 <- lm(Reaction ~ 0, data=sleepstudy)
LRT <- 2*c(logLik(fm, REML=FALSE) - logLik(lm1, REML=FALSE))

## Tests of ML-fits agree with anova():
fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
step(fm1)
lm2 <- lm(Reaction ~ Days, sleepstudy)
(an1 <- ranova(fm1))
(an2 <- anova(fm1, lm2))
j <- grep("Chi Df|Df", colnames(an2))
stopifnot(
  all.equal(an1[2, "LRT"], an2[2, "Chisq"], tolerance=TOL),
  all.equal(an1[2, "Df"], an2[2, j[length(j)]], tolerance=TOL),
  all.equal(an1[1:2, "logLik"], an2[2:1, "logLik"], tolerance=TOL)
)
## Note that lme4 version <1.1-22 use "Chi Df" while >=1.1-22 use "Df"

# Expect warnings when old (version < 3.0-0) arguments are used:
assertWarning(step(fm, reduce.fixed = FALSE, reduce.random = FALSE,
                   type=3, fixed.calc = FALSE, lsmeans.calc = FALSE,
                   difflsmeans.calc = TRUE, test.effs = 42, keep.e="save"))
assertWarning(step(fm, reduce.fixed = FALSE, reduce.random = FALSE,
                   lsmeans=3))


check_nrow <- function(obj, expect_nrow) {
  stopifnot(
    is.numeric(expect_nrow),
    nrow(obj) == expect_nrow
  )
}

# Statistical nonsense, but it works:
fm1 <- lmer(Reaction ~ Days + (1 | Subject) + (0 + Days|Subject), sleepstudy)
step(fm1)
(an <- ranova(fm1))
check_nrow(an, 3)
ranova(fm1, reduce.terms = FALSE)

# Statistical nonsense, but it works:
fm1 <- lmer(Reaction ~ Days + (0 + Days|Subject), sleepstudy)
step(fm1)
(an <- ranova(fm1)) # no test of non-nested models
stopifnot(
  nrow(an) == 2L,
  an[2L, "Df"] == 0,
  all(is.na(an[2L, "Pr(>Chisq)"]))
)
fm0 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
step(fm0)
(an2 <- anova(fm1, fm0, refit=FALSE))
stopifnot(
  (packageVersion("lme4")<="1.1.23" && an2[2L, "Pr(>Chisq)"] == 1) ||
    is.na(an2[2L, "Pr(>Chisq)"])
)
ranova(fm1, reduce.terms = FALSE)

fm1 <- lmer(Reaction ~ Days + (-1 + Days|Subject), sleepstudy)
step(fm1)
(an3 <- ranova(fm1)) # no test of non-nested models
stopifnot(
  all.equal(an, an3, check.attributes=FALSE, tolerance=TOL)
)

# Example where comparison of non-nested models is generated
fm <- lmer(Reaction ~ poly(Days, 2) + (0 + poly(Days, 2) | Subject), sleepstudy)
step(fm)
an <- ranova(fm)
stopifnot(
  nrow(an) == 2L,
  an[2, "Pr(>Chisq)"] == 1
)
ranova(fm, reduce.terms = FALSE) # test of nested models

# These models are nested, though:
fm <- lmer(Reaction ~ poly(Days, 2) + (1 + poly(Days, 2) | Subject), sleepstudy)
step(fm)
ranova(fm)
fm0 <- lmer(Reaction ~ poly(Days, 2) + (1 | Subject), sleepstudy)
step(fm0)
anova(fm0, fm, refit=FALSE)
ranova(fm, reduce.terms = FALSE)

# A model with ||-notation:
fm1 <- lmer(Reaction ~ Days + (Days||Subject), sleepstudy)
step(fm1)
ranova(fm1)

# What about models with nested factors?
fm <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor:TVset) + (1|Assessor),
           data=TVbo)
step(fm)
(an1 <- ranova(fm))

fm <- lmer(Coloursaturation ~ TVset * Picture +
             (1|Assessor/TVset), data=TVbo)
step(fm)
(an2 <- ranova(fm))
stopifnot(
  all.equal(an1, an2, check.attributes=FALSE, tolerance=TOL)
)

#####################################################################
# Test evaluation within functions, i.e. in other environments etc.
attach(sleepstudy)
fm <- lmer(Reaction ~ Days + (Days|Subject))
step(fm)
ranova(fm) # OK
detach(sleepstudy)

# Evaluating in a function works:
f <- function(form, data) lmer(form, data=data)
form <- "Informed.liking ~ Product+Information+
            (1|Consumer) + (1|Product:Consumer) + (1|Information:Consumer)"
fm <- f(form, data=ham)
ranova(fm)
step_res <- step(fm)
stopifnot(
  all(c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F value", "Pr(>F)") %in%
        colnames(step_res$fixed))
)

# Check that step works when form is a character vector
m <- lmer(form, data=ham)
step_res <- step(m)
(drop1_table <- attr(step_res, "drop1"))
stopifnot(
  all(c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F value", "Pr(>F)") %in%
        colnames(drop1_table))
)
# In version < 3.0-1.9002 attr(step_res, "drop1") picked up lme4::drop1.merMod
# and returned an AIC table after the model had been update'd.

#####################################################################
# Model with 2 ranef covarites:

# Model of the form (x1 + x2 | gr):
model <- lmer(Preference ~ sens2 + Homesize + (sens1 + sens2 | Consumer)
              , data=carrots)
step(model)
stopifnot(
  nrow(ranova(model)) == 3L,
  nrow(ranova(model, reduce.terms = FALSE)) == 2L
)

# Model of the form (f1 + f2 | gr):
model <- lmer(Preference ~ sens2 + Homesize + Gender +
                (Gender+Homesize|Consumer), data=carrots)
step(model)
stopifnot(
  nrow(ranova(model)) == 3L,
  nrow(ranova(model, reduce.terms = FALSE)) == 2L
)

# Model of the form (-1 + f2 | gr):
model <- lmer(Preference ~ sens2 + Homesize + Gender +
                (Gender -1 |Consumer), data=carrots)
step(model)
an1 <- ranova(model)
an1b <- ranova(model, reduce.terms = FALSE)

model <- lmer(Preference ~ sens2 + Homesize + Gender +
                (0 + Gender|Consumer), data=carrots)
step(model)
an2 <- ranova(model)
an2b <- ranova(model, reduce.terms = FALSE)

stopifnot(
  all.equal(an1, an2, check.attributes=FALSE, tolerance=TOL),
  all.equal(an1b, an2b, check.attributes=FALSE, tolerance=TOL)
)

####### Polynomial terms:
model <- lmer(Preference ~ sens2 + Gender + (poly(sens2, 2) | Consumer),
              data=carrots)
(an <- ranova(model))
step(model)

model <- lmer(Preference ~ sens2 + Gender + (sens2 + I(sens2^2) | Consumer),
              data=carrots)
(an2 <- ranova(model))
step(model)
stopifnot(
  nrow(an) == 2L,
  an[2L, "Df"] == 5L,
  nrow(an2) == 3L,
  all(an2[2:3, "Df"] == 3L)
)

######## Functions of terms in random effects:
model <- lmer(Preference ~ sens2 + Gender + (log(10+sens2) | Consumer),
              data=carrots)
ranova(model) # Works
step(model)

#####################################################################

# Missing values changes the number of observations in use:
m <- lmer(Preference ~ sens2 + Homesize +
            (1 |Consumer:Income), data=carrots)
assertError(step(m))
ans <- try(ranova(m), silent = TRUE)
stopifnot(
  inherits(ans, "try-error"),
  grepl("number of rows in use has changed", ans)
)

## Removing missing values solves the problem:
m2 <- lmer(Preference ~ sens2 + Homesize +
             (1 |Consumer:Income), data=carrots[complete.cases(carrots), ])
ranova(m2) # Works
step(m2)

## Including the variable with missing values (Income) among the fixed effects
## also solves the problem:
m <- lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
            (1 |Consumer:Income), data=carrots)
ranova(m)
step(m)

# Missing values in a an insignificant fixed effect causes the an error in step:
m0 <- lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
            (1 |Consumer), data=carrots)
ranova(m0)
ans <- try(step(m0), silent = TRUE)
stopifnot(
  inherits(ans, "try-error"),
  grepl("number of rows in use has changed", ans)
)

# Check that step still works for linear models (etc.)
flm <- lm(Coloursaturation ~ TVset * Picture, data=TVbo)
res <- step(flm, trace=0)
stopifnot(
  inherits(res, "lm")
)

##################### Using reduce and keep args:
# Fit a model to the ham dataset:
m <- lmer(Informed.liking ~ Product*Information+
            (1|Consumer) + (1|Product:Consumer)
          + (1|Information:Consumer), data=ham)

# Backward elimination using terms with default alpha-levels:
(step_res <- step(m))

(step_res <- step(m, reduce.random = FALSE))
(step_res <- step(m, reduce.fixed = FALSE))
(step_res <- step(m, reduce.fixed = FALSE, reduce.random = FALSE))

(step_res <- step(m, reduce.random = FALSE, keep="Information"))
(step_res <- step(m, reduce.random = FALSE, keep="Product:Information"))


###########################
## Test that `step` works even if all random terms are reduced away:
set.seed(101)
test <- data.frame(TM = factor(rep(rep(c("org","min"),each=3),3)),
                   dep = runif(18,0,20),
                   ind = runif(18,0,7),
                   dorp = factor(rep(1:3,each=6)))
full.model <- lmer(dep ~ TM + ind + (1 | dorp),  data=test)
res <- step(full.model)
# res$random
# res$fixed
# attr(res, "model")
# attr(res, "drop1")


