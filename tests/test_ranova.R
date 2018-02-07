# test_ranova.R
library(lmerTestR)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

#####################################################################
data("sleepstudy", package="lme4")

# Test reduction of (Days | Subject) to (1 | Subject):
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
(an <- ranova(fm1)) # 2 df test
stopifnot(
  nrow(an) == 2L,
  an[2L, "Df"] == 2L
)

# This test can also be achieved with anova():
fm2 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
(ana <- anova(fm1, fm2, refit=FALSE))

stopifnot(
  all.equal(an[2L, "LRT"], ana[2L, "Chisq"])
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
  all.equal(an[2L, "LRT"], LRT)
)

## _NULL_ model:
fm <- lmer(Reaction ~ -1 + (1|Subject), sleepstudy)
ranova(fm)
lm1 <- lm(Reaction ~ 0, data=sleepstudy)
LRT <- 2*c(logLik(fm, REML=FALSE) - logLik(lm1, REML=FALSE))

## Tests of ML-fits agree with anova():
fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
lm2 <- lm(Reaction ~ Days, sleepstudy)
(an1 <- ranova(fm1))
(an2 <- anova(fm1, lm2))
stopifnot(
  all.equal(an1[2, "LRT"], an2[2, "Chisq"]),
  all.equal(an1[2, "Df"], an2[2, "Chi Df"]),
  all.equal(an1[1:2, "logLik"], an2[2:1, "logLik"])
)

check_nrow <- function(obj, expect_nrow) {
  stopifnot(
    is.numeric(expect_nrow),
    nrow(obj) == expect_nrow
  )
}

# Statistical nonsense, but it works:
fm1 <- lmer(Reaction ~ Days + (1 | Subject) + (0 + Days|Subject), sleepstudy)
(an <- ranova(fm1))
check_nrow(an, 3)
ranova(fm1, reduce.terms = FALSE)

# Statistical nonsense, but it works:
fm1 <- lmer(Reaction ~ Days + (0 + Days|Subject), sleepstudy)
(an <- ranova(fm1)) # no test of non-nested models
stopifnot(
  nrow(an) == 2L,
  an[2L, "Df"] == 0,
  all(is.na(an[2L, "Pr(>Chisq)"]))
)
fm0 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
(an2 <- anova(fm1, fm0, refit=FALSE))
stopifnot(
  an2[2L, "Pr(>Chisq)"] == 1
)
ranova(fm1, reduce.terms = FALSE)

fm1 <- lmer(Reaction ~ Days + (-1 + Days|Subject), sleepstudy)
(an3 <- ranova(fm1)) # no test of non-nested models
stopifnot(
  all.equal(an, an3, check.attributes=FALSE)
)

# Example where comparison of non-nested models is generated
fm <- lmer(Reaction ~ poly(Days, 2) + (0 + poly(Days, 2) | Subject), sleepstudy)
an <- ranova(fm)
stopifnot(
  nrow(an) == 2L,
  an[2, "Pr(>Chisq)"] == 1
)
ranova(fm, reduce.terms = FALSE) # test of nested models

# These models are nested, though:
fm <- lmer(Reaction ~ poly(Days, 2) + (1 + poly(Days, 2) | Subject), sleepstudy)
ranova(fm)
fm0 <- lmer(Reaction ~ poly(Days, 2) + (1 | Subject), sleepstudy)
anova(fm0, fm, refit=FALSE)
ranova(fm, reduce.terms = FALSE)

# A model with ||-notation:
fm1 <- lmer(Reaction ~ Days + (Days||Subject), sleepstudy)
ranova(fm1)

# What about models with nested factors?
data("TVbo", package="lmerTest")
fm <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor:TVset) + (1|Assessor),
           data=TVbo)
(an1 <- ranova(fm))

fm <- lmer(Coloursaturation ~ TVset * Picture +
             (1|Assessor/TVset), data=TVbo)
(an2 <- ranova(fm))
stopifnot(
  all.equal(an1, an2, check.attributes=FALSE)
)

#####################################################################
# Test evaluation within functions, i.e. in other environments etc.
attach(sleepstudy)
fm <- lmer(Reaction ~ Days + (Days|Subject))
ranova(fm) # OK
detach(sleepstudy)

# Evaluating in a function leads to an error:
f <- function(form, data) lmer(form, data=data)
form <- formula("Reaction ~ Days + (Days|Subject)")
fm <- f(form, data=sleepstudy)
assertError(ranova(fm))

# Evaluating in parent environment does not help:
f <- function(form, data) eval.parent(lmer(form, data=data))
form <- formula("Reaction ~ Days + (Days|Subject)")
fm <- f(form, data=sleepstudy)
assertError(ranova(fm))

#####################################################################
# Model with 2 ranef covarites:

# Model of the form (x1 + x2 | gr):
data("carrots", package="lmerTest")
model <- lmer(Preference ~ sens2 + Homesize + (sens1 + sens2 | Consumer)
              , data=carrots)
stopifnot(
  nrow(ranova(model)) == 3L,
  nrow(ranova(model, reduce.terms = FALSE)) == 2L
)

# Model of the form (f1 + f2 | gr):
carrots$gender <- factor(carrots$Gender)
model <- lmer(Preference ~ sens2 + Homesize + gender +
                (gender+Homesize|Consumer), data=carrots)
stopifnot(
  nrow(ranova(model)) == 3L,
  nrow(ranova(model, reduce.terms = FALSE)) == 2L
)

# Model of the form (-1 + f2 | gr):
model <- lmer(Preference ~ sens2 + Homesize + gender +
                (gender -1 |Consumer), data=carrots)
an1 <- ranova(model)
an1b <- ranova(model, reduce.terms = FALSE)

model <- lmer(Preference ~ sens2 + Homesize + gender +
                (0 + gender|Consumer), data=carrots)
an2 <- ranova(model)
an2b <- ranova(model, reduce.terms = FALSE)

stopifnot(
  all.equal(an1, an2, check.attributes=FALSE),
  all.equal(an1b, an2b, check.attributes=FALSE)
)

####### Polynomial terms:
model <- lmer(Preference ~ sens2 + gender + (poly(sens2, 2) | Consumer),
              data=carrots)
(an <- ranova(model))

model <- lmer(Preference ~ sens2 + gender + (sens2 + I(sens2^2) | Consumer),
              data=carrots)
(an2 <- ranova(model))
stopifnot(
  nrow(an) == 2L,
  an[2L, "Df"] == 5L,
  nrow(an2) == 3L,
  all(an2[2:3, "Df"] == 3L)
)

######## Functions of terms in random effects:
model <- lmer(Preference ~ sens2 + gender + (log(10+sens2) | Consumer),
              data=carrots)
ranova(model) # Works

#####################################################################

# Missing values changes the number of observations in use:
m <- lmer(Preference ~ sens2 + Homesize +
            (1 |Consumer:Income), data=carrots)
ans <- try(ranova(m), silent = TRUE)
stopifnot(
  inherits(ans, "try-error"),
  grepl("number of rows in use has changed", ans)
)

## Removing missing values solves the problem:
m2 <- lmer(Preference ~ sens2 + Homesize +
             (1 |Consumer:Income), data=carrots[complete.cases(carrots), ])
ranova(m2) # Works

## Including the variable with missing values (Income) among the fixed effects
## also solves the problem:
m <- lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
            (1 |Consumer:Income), data=carrots)
ranova(m)

