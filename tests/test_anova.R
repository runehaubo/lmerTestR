# test_anova.R
library(lmerTestR)
stopifnot(require(tools)) # For assertError and assertWarning
data("sleepstudy", package="lme4")
data("cake", package="lme4")

####################################
## Basic anova tests
####################################

m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

####### ddf argument:
(an1 <- anova(m)) # Also testing print method.
(an2 <- anova(m, ddf="Satterthwaite"))
stopifnot(isTRUE(
  all.equal(an1, an2)
))
(an3 <- anova(m, ddf="Sat")) ## Abbreviated argument
stopifnot(isTRUE(
  all.equal(an1, an3)
))
(anova(m, ddf="KR"))
(anova(m, ddf="lme4"))
assertError(anova(m, ddf="Ken")) ## Error on incorrect arg.

## lme4 method:
an1 <- anova(m, ddf="lme4")
an2 <- anova(as(m, "lmerMod"))
stopifnot(isTRUE(
  all.equal(an1, an2)
))

###### type argument:
(an1 <- anova(m, type="1")) # valid type arg.
(an2 <- anova(m, type="I")) # same
stopifnot(isTRUE(
  all.equal(an1, an2)
))
(an3 <- anova(m, type=1)) # Not strictly valid, but accepted
stopifnot(isTRUE(
  all.equal(an1, an3)
))
assertError(anova(m, type=0)) # Not valid arg.
assertError(anova(m, type="i")) # Not valid arg.

####################################
## Example with factor fixef:
####################################

## 'temp' is continuous, 'temperature' an ordered factor with 6 levels
m <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
(an <- anova(m))
(an_KR <- anova(m, ddf="KR"))
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")])
stopifnot(isTRUE(res))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_KR[, c("Sum Sq", "Mean Sq", "F value")])
stopifnot(isTRUE(res))

stopifnot(all.equal(c(2, 5, 10), an$NumDF),
          all.equal(c(42, 210, 210), an$DenDF))

# No intercept:
m <- lmer(angle ~ 0 + recipe * temperature + (1|recipe:replicate), cake)
(an <- anova(m))
(an_KR <- anova(m, ddf="KR"))
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")])
stopifnot(isTRUE(res))

# ML-fit:
m <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake, REML=FALSE)
(an <- anova(m))
assertError(an <- anova(m, ddf="KR")) # KR fits should be REML
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")])
stopifnot(isTRUE(res))


####################################
## Example with continuous fixef:
####################################

# Example with no fixef:
m <- lmer(Reaction ~ -1 + (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 0 + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 0L)
an <- anova(m)
an_KR <- anova(m, ddf="KR")
stopifnot(nrow(an) == 0L,
          nrow(an_KR) == 0L)
# anova(m, ddf="lme4") # Bug in lme4 it seems

# Example with intercept only:
m <- lmer(Reaction ~ (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 1 + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 1L,
          names(fixef(m)) == "(Intercept)")
(an <- anova(m))
(an_KR <- anova(m, ddf="KR"))
(an_lme4 <- anova(m, ddf="lme4"))
stopifnot(nrow(an) == 0L,
          nrow(an_lme4) == 0L,
          nrow(an_KR) == 0L)

# Example with 1 fixef without intercept:
m <- lmer(Reaction ~ Days - 1 + (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 0 + Days + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 1L,
          names(fixef(m)) == "Days")
(an <- anova(m))
(an_KR <- anova(m, ddf="KR"))
(an_lme4 <- anova(m, ddf="lme4"))

stopifnot(nrow(an) == 1L,
          nrow(an_KR) == 1L,
          nrow(an_lme4) == 1L)
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")])
stopifnot(isTRUE(res))
stopifnot(isTRUE(all.equal(
  c(1, 17), unname(unlist(an[, c("NumDF", "DenDF")])),
  tolerance=1e-4
)))

# Example with >1 fixef without intercept:
m <- lmer(Reaction ~ Days - 1 + I(Days^2) + (Days | Subject), sleepstudy)
stopifnot(length(fixef(m)) == 2L,
          names(fixef(m)) == c("Days", "I(Days^2)"))
(an <- anova(m))
(an_KR <- anova(m, ddf="KR"))
(an_lme4 <- anova(m, ddf="lme4"))
stopifnot(nrow(an) == 2L,
          nrow(an_KR) == 2L,
          nrow(an_lme4) == 2L)
# Here is a diff in SSQ which doesn't seem well-defined anyway...
# SSQ for I(Days^2) agree though.
# t-statistics also agree:
coef(summary(m))
Lmat <- diag(length(fixef(m)))
lmerTestR:::rbindall(lapply(1:nrow(Lmat), function(i) contest1D(Lmat[i, ], m)))

# Example with >1 fixef and intercept:
m <- lmer(Reaction ~ Days + I(Days^2) + (Days | Subject), sleepstudy)
stopifnot(length(fixef(m)) == 3L)
(an <- anova(m))
(an_KR <- anova(m, ddf="KR"))
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")])
stopifnot(isTRUE(res))

res <- all.equal(an[, c("Sum Sq", "Mean Sq", "DenDF", "F value")],
                 an_KR[, c("Sum Sq", "Mean Sq", "DenDF", "F value")], tol=1e-6)
stopifnot(isTRUE(res))
