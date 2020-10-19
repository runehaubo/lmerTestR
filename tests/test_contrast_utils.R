# test_contrast_utils.R

library(lmerTest)

##########
# Test that a message is printed if some cells have zero data:
# Missing a single cell:
data("cake", package="lme4")
cake4 <- cake
cake4$temperature <- factor(cake4$temperature, ordered=FALSE)
cake4 <- droplevels(subset(cake4, !(recipe == "A" & temperature == "175") ))
with(cake4, table(recipe, temperature))

fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4)
an <- anova(fm1)
txt <- capture.output(an <- anova(fm1), type = "message")
stopifnot(length(grep("Missing cells for:", txt)) > 0,
          length(grep("Interpret type III hypotheses with care.", txt)) > 0)

##########
# Test that a message is printed if some cells have zero data:
# Missing diagonal:
cake4 <- cake
cake4$temperature <- factor(cake4$temperature, ordered=FALSE)
cake4 <- droplevels(subset(cake4, temperature %in% levels(cake4$temperature)[1:3]))
cake4 <- droplevels(subset(cake4, !((recipe == "A" & temperature == "175") |
                                      (recipe == "B" & temperature == "185") |
                                      (recipe == "C" & temperature == "195") )))
cake4$temp0 <- cake4$temp - mean(cake4$temp)
with(cake4, table(recipe, temperature))

fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4)
an <- anova(fm1)
txt <- capture.output(an <- anova(fm1), type = "message")
stopifnot(length(grep("Missing cells for:", txt)) > 0,
          length(grep("Interpret type III hypotheses with care.", txt)) > 0)

##########
# Test that a message is NOT printed with centered covariates:
fm1 <- lmer(angle ~ recipe * temp0 + (1|recipe:replicate), cake4)
an <- anova(fm1)
txt <- capture.output(an <- anova(fm1), type = "message")
stopifnot(length(grep("Missing cells for:", txt)) == 0,
          length(grep("Interpret type III hypotheses with care.", txt)) == 0)
# Note: in many cases a message would not be printed anyway because the
# columns sums in the rdX design matrix would not be exactly zero but just a
# small number very close to zero.



