# test_compare_sas.R
library(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

#####################################################################


# Use contrasts to get particular estimates for the summary table:
l <- list(Frequency="contr.SAS", Income="contr.SAS")
m.carrots <- lmer(Preference ~ sens2*Frequency*Income
                  +(1+sens2|Consumer), data=carrots, contrasts=l)
an.m <- anova(m.carrots)

TOL <- 1e-4
# with 4 decimals should agree with SAS output
# numbers before decimals should agree with SAS output
stopifnot(
  all.equal(an.m[,"Pr(>F)"],
            c(2e-5, 0.15512,  0.06939, 0.08223, 0.52459, 0.03119, 0.48344),
            tol = TOL),
  all.equal(round(an.m$DenDF), c(83, 83, 83, 83, 83, 83, 83))
)

sm <- summary(m.carrots)
stopifnot(
  isTRUE(all.equal(sm$coefficients[,"Pr(>|t|)"],
                   c(1e-10, 0.005061, 0.6865554, 0.342613, 0.129157,
                     0.088231, 0.846000, 0.354472, 0.526318, 0.020646, 0.010188,
                     0.031242, 0.055356, 0.694689, 0.099382, 0.28547,
                     0.977774, 0.855653, 0.427737, 0.321086, 0.417465 , 0.204385, 0.784437,
                     0.681434, 0.106180, 0.149122, 0.390870, 0.273686), tol=TOL,
                   check.attributes = FALSE))
)

# Takes too long to run:
# if(requireNamespace("pbkrtest", quietly = TRUE)) {
#   sm.kr <- summary(m.carrots, ddf = "Kenward-Roger")
#
#   ## coefficients for Sat and KR agree in this example
#   #   cbind(sm$coefficients[,"Pr(>|t|)"], sm.kr$coefficients[,"Pr(>|t|)"])
#   all.equal(sm$coefficients[,"Pr(>|t|)"], sm.kr$coefficients[,"Pr(>|t|)"],
#             tol=TOL)
# }

################################################################################
## checking lsmeans and difflsmeans
## compare with SAS output
m <- lmer(Informed.liking ~ Product*Information*Gender
          + (1|Product:Consumer) + (1|Consumer) , data=ham)


lsm <- lsmeansLT(m, which = "Product")
head(lsm)

TOL <- 1e-5
stopifnot(
  isTRUE(all.equal(lsm[, "Estimate"], c(5.8084, 5.1012, 6.0909, 5.9256),
                   tol=TOL, check.attributes = FALSE)),
  isTRUE(all.equal(round(lsm[, "t value"], 2), c(24.93, 21.89, 26.14, 25.43), tol=TOL,
                   check.attributes = FALSE)),
  isTRUE(all.equal(lsm[, "lower"], c(5.3499, 4.6428, 5.6324, 5.4672), tol=TOL,
                   check.attributes = FALSE)),
  isTRUE(all.equal(lsm[, "upper"], c(6.2668, 5.5597, 6.5493, 6.3840), tol=TOL,
                   check.attributes = FALSE))
)

################################################################################
# Not actually 'hard-coded' tests versus SAS results...

m.carrots <- lmer(Preference ~ 0 + sens2 + Homesize +
                    (1+sens2 | Consumer), data=carrots)
summary(m.carrots)

(an.1 <- anova(m.carrots, type=1))
(an.3 <- anova(m.carrots))
(an.lme4 <- anova(m.carrots, ddf = "lme4")) # difference in SSQ MS and F-values
# Is this a problem with lme4?
# fm <- lm(Preference ~ 0 + sens2 + Homesize, data=carrots)
# anova(fm)
# coef(summary(fm))
# Here the F value is a little greater than the squared t-value (as expected)

TOL <- 1e-5
stopifnot(all.equal(an.1[, "F value"], c(56.5394, 4169.87), tol = TOL),
          all.equal(an.3[, "F value"], c(54.8206, 4169.87), tol = TOL))


################################################################################
# Check exmaple from GLM SAS report

### example from the paper GLM SAS 101 report
a <- factor(c(1,1,1,2,2,2,2,2,1,2))
b <- factor(c(1,1,2,1,2,2,2,2,2,1))
f=factor(c(1,2,1,2,1,2,1,2,1,2))
y <- c(12,14,11,20,17,23,35,46,15,16)
dd <- data.frame(a=a, b=b, y=y, f=f)

## check type 2 is order independent
model <- lmer(y ~ a*b + (1|f), data=dd)
model2 <- lmer(y ~ b*a + (1|f), data=dd)
(an <- anova(model, type=2))
(an2 <- anova(model2, type=2))
TOL <- 1e-5
stopifnot(
  isTRUE(all.equal(an,an2[c(2,1,3),], check.attributes = FALSE, tol=TOL))
)

## check the results are the same as from SAS proc mixed
stopifnot(
  isTRUE(all.equal(an[,"F value"], c(3.90131, 1.32753, 0.99565), tol = 1e-5))
)
################################################################################
## Check type II and III anova tables versus SAS

m.carrots <- lmer(Preference ~ sens2*Homesize
                  +(1+sens2|Consumer), data=carrots)
(ancar <- anova(m.carrots, type=2))

stopifnot(
  isTRUE(all.equal(ancar[,"F value"], c(54.8361, 5.16138, 1.03035), tol = 1e-4))
)

m <- lmer(Informed.liking ~ Product*Age
          + (1|Consumer) , data=ham)
(an <- anova(m, type=2))

stopifnot(
  isTRUE(all.equal(an[,"F value"], c(2.48135, .005387, 1.48451), tol = 1e-5))
)


fm <- lmer(Preference ~ sens2*Homesize*sens1 + (1|Product),
           data=carrots)
(ant2 <- anova(fm, type=2))
(ant3 <- anova(fm, type=3))

stopifnot(
  isTRUE(all.equal(ant2[,"F value"],
                   c(16.4842, 14.0010, .526076, 1.18144,
                     .107570, .335177, 1.05946), tol = 1e-4)),
  isTRUE(all.equal(ant3[,"F value"],
                   c(16.9140, 14.0010,.481148, 1.18144,
                     .074201, .335177, 1.05946), tol = 1e-4))
)

################################################################################


