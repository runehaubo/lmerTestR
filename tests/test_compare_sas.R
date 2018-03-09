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
  isTRUE(all.equal(sm$coefficients[,"Pr(>|t|)"], c(1e-10, 0.005061, 0.6865554, 0.342613, 0.129157,
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
################################################################################


