library(lmerTest)

# Read in data set
load(system.file("testdata","potdata.RData", package="lmerTest"))

# Mixed model
lmerout <- lmer(biomass ~ CO2*nutrients + (1|chamber),data=potdata)
summary(lmerout)

an.sat <- anova(lmerout)
anova(lmerout, ddf="lme4")
TOL <- 1e-5
stopifnot(isTRUE(all.equal(
  an.sat[,"DenDF"], c(2, 10, 10), tol=TOL
)))

stopifnot(isTRUE(
  all.equal(an.sat[,"Pr(>F)"], c(0.0224955602, 1e-11, 0.020905569), tol=TOL)
))

# if(require(pbkrtest))
#   an.kr <- anova(lmerout, ddf="Kenward-Roger")
#
# TOL <- 1e-7
# stopifnot(all.equal(an.kr[,"Pr(>F)"], c(0.0224955602, 1e-11, 0.020905569) ,
#                     tol=TOL),
#           all.equal(an.kr[,"DenDF"],
#                     c(2, 10, 10) , tol=TOL),
#           TRUE)
