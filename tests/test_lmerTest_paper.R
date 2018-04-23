# test_lmerTest_paper.R

library(lmerTest)

# Kenward-Roger only available with pbkrtest and only then validated in R >= 3.3.3
# (faulty results for R < 3.3.3 may be due to unstated dependencies in pbkrtest)
has_pbkrtest <- requireNamespace("pbkrtest", quietly = TRUE) && getRversion() >= "3.3.3"

# Read in data set
load(system.file("testdata","test_paper_objects.RData", package="lmerTest"))

# Evaluate code from paper:
## Section 8.2:
tv <- lmer(Sharpnessofmovement ~ TVset * Picture + (1 | Assessor) +
             (1 | Assessor:TVset) + (1 | Assessor:Picture), data = TVbo)
(an8.2 <- anova(tv))

if(has_pbkrtest)
  (ankr8.2 <- anova(tv, type=2, ddf="Kenward-Roger"))

## Section 8.3:
m.carrots <- lmer(Preference ~ sens1 + sens2 + (1 + sens1 + sens2 | Consumer) +
                    (1 | Product), data=carrots)
(sum8.3 <- coef(summary(m.carrots)))

## Section 8.4:
tv <- lmer(Sharpnessofmovement ~ TVset * Picture +
             (1 | Assessor:TVset) + (1 | Assessor:Picture) +
             (1 | Assessor:Picture:TVset) + (1 | Repeat) + (1 | Repeat:Picture) +
             (1 | Repeat:TVset) + (1 | Repeat:TVset:Picture) + (1 | Assessor),
           data = TVbo)
st <- step(tv)
names(st)
(elim_tab_random8.4 <- st$random)
(elim_tab_fixed8.4 <- st$fixed)
(an8.4 <- anova(get_model(st)))

## Section 8.5:
# L <- matrix(0, ncol = 12, nrow = 6)
# L[1, 7] <- L[2, 8] <- L[3, 9] <- L[4, 10] <- L[5, 11] <- L[6, 12] <- 1
L <- cbind(array(0, dim=c(6, 6)), diag(6))
(con1_8.5 <- calcSatterth(tv, L))
(con2_8.5 <- contest(tv, L))

## Section C:
# m.carrots <- lmer(Preference ~ sens1 + sens2 + (1 + sens1 + sens2 | Consumer) +
#                     (1 | product), data = carrots)
# step(m.carrots, reduce.fixed = FALSE)
(ran_C <- ranova(m.carrots))

# Compare to validated outputs:
TOL <- 1e-4
stopifnot(
  isTRUE(all.equal(an8.2_save, an8.2, check.attributes = FALSE, tolerance=TOL)),
  isTRUE(all.equal(sum8.3_save, sum8.3, check.attributes = FALSE, tolerance=TOL)),
  isTRUE(all.equal(elim_tab_random8.4_save, elim_tab_random8.4,
                   check.attributes = FALSE, tolerance=TOL)),
  isTRUE(all.equal(elim_tab_fixed8.4_save, elim_tab_fixed8.4,
                   check.attributes = FALSE, tolerance=TOL)),
  isTRUE(all.equal(an8.4_save, an8.4, check.attributes = FALSE, tolerance=TOL)),
  isTRUE(all.equal(con1_8.5_save, con1_8.5, check.attributes = FALSE, tolerance=TOL)),
  isTRUE(all.equal(con2_8.5_save, con2_8.5, check.attributes = FALSE, tolerance=TOL))
)
if(has_pbkrtest) {
  stopifnot(
    isTRUE(all.equal(ankr8.2_save, ankr8.2, check.attributes = FALSE, tolerance=TOL))
  )
}

