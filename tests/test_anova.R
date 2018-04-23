# test_anova.R
library(lmerTest)

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
## Basic anova tests
####################################

m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

####### ddf argument:
(an1 <- anova(m)) # Also testing print method.
(an2 <- anova(m, ddf="Satterthwaite"))
(an2b <- anova(m, ddf="Satterthwaite", type=3))
(an2c <- anova(m, ddf="Satterthwaite", type=2))
stopifnot(isTRUE(
  all.equal(an1, an2, tolerance=TOL)
))
(an3 <- anova(m, ddf="Sat")) ## Abbreviated argument
stopifnot(isTRUE(
  all.equal(an1, an3, tolerance=TOL)
))
if(has_pbkrtest) {
  (anova(m, ddf="Kenward-Roger"))
  (anova(m, ddf="Kenward-Roger", type=3))
}
(an1 <- anova(m, ddf="lme4"))
(an2 <- anova(m, ddf="lme4", type=3)) # 'type' is ignored with ddf="lme4"
stopifnot(isTRUE(
  all.equal(an1, an2, tolerance=TOL)
))
res <- assertError(anova(m, ddf="KR")) ## Error on incorrect arg.
stopifnot(
  grepl("'arg' should be one of ", unlist(res[[1]])$message)
)

## lme4 method:
an1 <- anova(m, ddf="lme4")
an2 <- anova(as(m, "lmerMod"))
stopifnot(isTRUE(
  all.equal(an1, an2, tolerance=TOL)
))

###### type argument:
(an1 <- anova(m, type="1")) # valid type arg.
(an2 <- anova(m, type="I")) # same
stopifnot(isTRUE(
  all.equal(an1, an2, tolerance=TOL)
))
(an3 <- anova(m, type=1)) # Not strictly valid, but accepted
stopifnot(isTRUE(
  all.equal(an1, an3, tolerance=TOL)
))

(an1 <- anova(m, type="2")) # valid type arg.
(an2 <- anova(m, type="II")) # same
stopifnot(isTRUE(
  all.equal(an1, an2, tolerance=TOL)
))
(an3 <- anova(m, type=3)) # Not strictly valid, but accepted
stopifnot(isTRUE(
  all.equal(an1, an3, check.attributes=FALSE, tolerance=TOL)
))

(an1 <- anova(m, type="3")) # valid type arg.
(an2 <- anova(m, type="III")) # same
stopifnot(isTRUE(
  all.equal(an1, an2, tolerance=TOL)
))
(an3 <- anova(m, type=3)) # Not strictly valid, but accepted
stopifnot(isTRUE(
  all.equal(an1, an3, tolerance=TOL)
))
assertError(anova(m, type=0)) # Not valid arg.
assertError(anova(m, type="i")) # Not valid arg.

####### Model comparison:
fm <- lm(Reaction ~ Days, sleepstudy)
(an <- anova(m, fm))
stopifnot(
  nrow(an) == 2L,
  rownames(an)[2] == "m"
)

m2 <- lmer(Reaction ~ Days + I(Days^2) + (Days | Subject), sleepstudy)
(an <- anova(m, m2, refit=FALSE))
stopifnot(
  nrow(an) == 2L,
  rownames(an)[1] == "m"
)


####################################
## Example with factor fixef:
####################################

## 'temp' is continuous, 'temperature' an ordered factor with 6 levels
data("cake", package="lme4")
m <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
(an <- anova(m))
(an_lme4 <- anova(m, ddf="lme4"))

if(has_pbkrtest) {
  (an_KR <- anova(m, ddf="Kenward-Roger"))
  # res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
  #                  an_lme4[, c("Sum Sq", "Mean Sq", "F value")])
  # stopifnot(isTRUE(res))
  res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                   an_KR[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
  stopifnot(isTRUE(res))
}
stopifnot(all.equal(c(2, 1, 2), an$NumDF, tol=1e-6),
          all.equal(c(254.0157612, 222, 222), an$DenDF, tol=TOL))

an3 <- anova(m, type=3)
an2 <- anova(m, type=2)
an1 <- anova(m, type=1)

## Data is balanced, so Type II and III should be identical:
## One variable is continuous, so Type I and II/III are different:
stopifnot(
  isTRUE(all.equal(an3, an2, check.attributes=FALSE, tolerance=TOL)),
  !isTRUE(all.equal(an1, an2, check.attributes=FALSE, tolerance=1e-8))
)

# Using an ordered factor:
m <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
(an1 <- anova(m, type=1))
(an2 <- anova(m, type=2))
# Type 3 is also available with ordered factors:
(an3 <- anova(m, type=3))
## Balanced data and only factors: Type I, II and III should be the same:
stopifnot(
  isTRUE(all.equal(an1, an2, check.attributes=FALSE, tolerance=TOL)),
  isTRUE(all.equal(an1, an3, check.attributes=FALSE, tolerance=TOL))
)

(an <- anova(m, type=1))
(an_lme4 <- anova(m, type=1, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
stopifnot(isTRUE(res))
if(has_pbkrtest) {
  (an_KR <- anova(m, type=1, ddf="Kenward-Roger"))
  res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                   an_KR[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
  stopifnot(isTRUE(res))
}
stopifnot(all.equal(c(2, 5, 10), an$NumDF, tolerance=TOL),
          all.equal(c(42, 210, 210), an$DenDF, tolerance=TOL))

########
## Make case with balanced unordered factors:
cake2 <- cake
cake2$temperature <- factor(cake2$temperature, ordered = FALSE)
# str(cake2)
stopifnot(
  !is.ordered(cake2$temperature)
)
m <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake2)
(an1 <- anova(m, type=1))
(an2 <- anova(m, type=2))
(an3 <- anova(m, type=3))
## Balanced data and only factors: Type I, II, and III should be the same:
stopifnot(
  isTRUE(all.equal(an1, an2, check.attributes=FALSE, tolerance=TOL)),
  isTRUE(all.equal(an3, an2, check.attributes=FALSE, tolerance=TOL))
)
########

# No intercept:
m <- lmer(angle ~ 0 + recipe * temp + (1|recipe:replicate), cake)
(an <- anova(m, type=1))
(an2 <- anova(m, type=2))
(an2 <- anova(m, type=3))
if(has_pbkrtest)
  (an_KR <- anova(m, ddf="Kenward-Roger"))
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
stopifnot(isTRUE(res))

# ML-fit:
m <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake, REML=FALSE)
(an <- anova(m, type=1))
if(has_pbkrtest)
  assertError(an <- anova(m, ddf="Kenward-Roger")) # KR fits should be REML
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
stopifnot(isTRUE(res))

####################################
## Using contr.sum:
####################################

m <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake,
          contrasts = list('recipe' = "contr.sum"))
(an <- anova(m, type=1))
(an2 <- anova(m, type=2))
(an3 <- anova(m, type=3))
stopifnot(
  isTRUE(all.equal(an2, an3, check.attributes=FALSE, tolerance=TOL))
)
if(has_pbkrtest)
  (an_KR <- anova(m, type=1, ddf="Kenward-Roger"))
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
stopifnot(isTRUE(res))


####################################
## Example with continuous fixef:
####################################

# Example with no fixef:
m <- lmer(Reaction ~ -1 + (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 0 + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 0L)
(an <- anova(m, type=1))
(an_2 <- anova(m, type=2))
(an_3 <- anova(m, type=3))
stopifnot(nrow(an) == 0L,
          nrow(an_2) == 0L,
          nrow(an_3) == 0L)
# anova(m, ddf="lme4") # Bug in lme4 it seems
if(has_pbkrtest) {
  (an_KR <- anova(m, ddf="Kenward-Roger"))
  stopifnot(
    nrow(an_KR) == 0L
  )
}

# Example with intercept only:
m <- lmer(Reaction ~ (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 1 + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 1L,
          names(fixef(m)) == "(Intercept)")
(an <- anova(m))
(an_2 <- anova(m, type=2))
(an_3 <- anova(m, type=3))
(an_lme4 <- anova(m, ddf="lme4"))
stopifnot(nrow(an) == 0L,
          nrow(an_2) == 0L,
          nrow(an_3) == 0L,
          nrow(an_lme4) == 0L)
if(has_pbkrtest) {
  (an_KR <- anova(m, ddf="Kenward-Roger"))
  stopifnot(
    nrow(an_KR) == 0L
  )
}

# Example with 1 fixef without intercept:
m <- lmer(Reaction ~ Days - 1 + (Days | Subject), sleepstudy)
# m <- lmer(Reaction ~ 0 + Days + (Days | Subject), sleepstudy) # alternative
stopifnot(length(fixef(m)) == 1L,
          names(fixef(m)) == "Days")
(an <- anova(m))
(an_2 <- anova(m, type=2))
(an_3 <- anova(m, type=3))
(an_lme4 <- anova(m, ddf="lme4"))
stopifnot(nrow(an) == 1L,
          nrow(an_2) == 1L,
          nrow(an_3) == 1L,
          nrow(an_lme4) == 1L)
if(has_pbkrtest) {
  (an_KR <- anova(m, ddf="Kenward-Roger"))
  stopifnot(
    nrow(an_KR) == 1L
  )
}

res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
stopifnot(isTRUE(res))
stopifnot(isTRUE(all.equal(
  c(1, 17), unname(unlist(an[, c("NumDF", "DenDF")])), tolerance=TOL
)))

# Example with >1 fixef without intercept:
m <- lmer(Reaction ~ Days - 1 + I(Days^2) + (Days | Subject), sleepstudy)
stopifnot(length(fixef(m)) == 2L,
          names(fixef(m)) == c("Days", "I(Days^2)"))
(an <- anova(m))
(an_2 <- anova(m, type=2))
(an_3 <- anova(m, type=3))
(an_lme4 <- anova(m, ddf="lme4"))
stopifnot(nrow(an) == 2L,
          nrow(an_3) == 2L,
          nrow(an_lme4) == 2L)
if(has_pbkrtest) {
  (an_KR <- anova(m, ddf="Kenward-Roger"))
  stopifnot(
    nrow(an_KR) == 2L
  )
}
# Here is a diff in SSQ which doesn't seem well-defined anyway...
# SSQ for I(Days^2) agree though.
# t-statistics also agree:
coef(summary(m))
Lmat <- diag(length(fixef(m)))
lmerTest:::rbindall(lapply(1:nrow(Lmat), function(i) contest1D(m, Lmat[i, ])))

# Example with >1 fixef and intercept:
m <- lmer(Reaction ~ Days + I(Days^2) + (Days | Subject), sleepstudy)
stopifnot(length(fixef(m)) == 3L)
(an <- anova(m, type=1))
(an_2 <- anova(m, type=2))
(an_3 <- anova(m, type=3))
(an_lme4 <- anova(m, ddf="lme4"))
res <- all.equal(an[, c("Sum Sq", "Mean Sq", "F value")],
                 an_lme4[, c("Sum Sq", "Mean Sq", "F value")], tolerance=TOL)
stopifnot(isTRUE(res))

if(has_pbkrtest) {
  (an_KR <- anova(m, ddf="Kenward-Roger"))
  res <- all.equal(an_3[, c("Sum Sq", "Mean Sq", "DenDF", "F value")],
                   an_KR[, c("Sum Sq", "Mean Sq", "DenDF", "F value")], tolerance=TOL)
  stopifnot(isTRUE(res))
}

## FIXME: Test the use of refit arg to lme4:::anova.merMod

##############################
# Test that type III anova is the same regardless of contrast coding:
# 3 x 3 factorial with missing diagonal
data("cake", package="lme4")
cake4 <- cake
cake4$temperature <- factor(cake4$temperature, ordered=FALSE)
cake4 <- droplevels(subset(cake4, temperature %in% levels(cake4$temperature)[1:3]))
cake4 <- droplevels(subset(cake4, !((recipe == "A" & temperature == "175") |
                                      (recipe == "B" & temperature == "185") |
                                      (recipe == "C" & temperature == "195") )))
str(cake4)
with(cake4, table(recipe, temperature))
# load_all(r2path)

fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4)
fm2 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4,
            contrasts=list(recipe="contr.sum", temperature="contr.SAS"))
fm3 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4,
            contrasts=list(recipe="contr.sum", temperature="contr.poly"))
fm4 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4,
            contrasts=list(recipe=contr.helmert, temperature="contr.poly"))
(an1 <- anova(fm1))
(an2 <- anova(fm2))
(an3 <- anova(fm3))
(an4 <- anova(fm4))
options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
fm5 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4)
(an5 <- anova(fm5))
options(contrasts = c("contr.treatment", "contr.poly"))
options("contrasts")
stopifnot(
  isTRUE(all.equal(an1, an2, check.attributes=FALSE, tolerance=TOL)),
  isTRUE(all.equal(an1, an3, check.attributes=FALSE, tolerance=TOL)),
  isTRUE(all.equal(an1, an4, check.attributes=FALSE, tolerance=TOL)),
  isTRUE(all.equal(an1, an5, check.attributes=FALSE, tolerance=TOL))
)
