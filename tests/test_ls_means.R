# test_lsmeans.R

library(lmerTest)

TOL <- 1e-4
# Kenward-Roger only available with pbkrtest and only then validated in R >= 3.3.3
# (faulty results for R < 3.3.3 may be due to unstated dependencies in pbkrtest)
has_pbkrtest <- requireNamespace("pbkrtest", quietly = TRUE) && getRversion() >= "3.3.3"

########### Basic model structures:

# Factor * covariate:
data("cake", package="lme4")
model <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
(lsm <- ls_means(model))
stopifnot(
  nrow(lsm) == 3L,
  ncol(lsm) == 7L,
  # Balanced, so LS-means equal raw means:
  isTRUE(all.equal(c(with(cake, tapply(angle, recipe, mean))), lsm[, "Estimate"],
                   check.attributes=FALSE, tolerance=TOL))
)

# Pairwise differences of LS-means:
plsm <- ls_means(model, pairwise = TRUE)
plsm2 <- difflsmeans(model)
C <- as.matrix(lmerTest:::get_pairs(rownames(lsm)))
stopifnot(
  isTRUE(all.equal(plsm, plsm2, tolerance=TOL)),
  isTRUE(all.equal(plsm[, "Estimate"], c(lsm[, "Estimate"] %*% C),
                   check.attributes=FALSE, tolerance=TOL))
)

# Contrasts vectors:
show_tests(lsm)
show_tests(plsm)

# Factor * Ordered:
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
(lsm2 <- ls_means(model))
stopifnot(
  nrow(lsm2) == 3 + 6 + 3*6,
  ncol(lsm) == 7L,
  # Balanced, so LS-means equal raw means:
  isTRUE(all.equal(lsm[1:3, ], lsm2[1:3, ],
                   check.attributes=FALSE, tolerance=TOL))
)


# Factor * Factor:
cake2 <- cake
cake2$temperature <- factor(cake2$temperature, ordered = FALSE)
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake2)
(lsm3 <- ls_means(model))
stopifnot(
  isTRUE(all.equal(lsm2, lsm3, check.attributes=FALSE, tolerance=TOL))
)

# Covariate (only):
data("sleepstudy", package="lme4")
m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
(lsm <- ls_means(m))
stopifnot(
  nrow(lsm) == 0L,
  ncol(lsm) == 7L
)

# No fixef:
m <- lmer(Reaction ~ 0 + (Days | Subject), sleepstudy)
(lsm <- ls_means(m))
stopifnot(
  nrow(lsm) == 0L,
  ncol(lsm) == 7L
)

########### Arguments and options:

# which
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake2)
(lsm4 <- ls_means(model, which = "recipe"))
stopifnot(
  nrow(lsm4) == 3L,
  ncol(lsm4) == 7L,
  isTRUE(all.equal(lsm3[1:3, ], lsm4, check.attributes=FALSE, tolerance=TOL))
)

# KR:
if(has_pbkrtest)
  (lsm5 <- ls_means(model, which = "recipe", ddf = "Kenward-Roger"))

# level:
(lsm6 <- ls_means(model, which = "recipe", level=0.99))

stopifnot(
  all(lsm6[, "lower"] < lsm4[, "lower"]),
  all(lsm6[, "upper"] > lsm4[, "upper"])
)



########### Missing cels -> unestimable contrasts:

# Missing cell:
cake3 <- cake
cake3$temperature <- factor(cake3$temperature, ordered=FALSE)
cake3 <- droplevels(subset(cake3, temperature %in% levels(cake3$temperature)[1:3]))
cake3 <- droplevels(subset(cake3, !(recipe == "C" & temperature == "195") ))
str(cake3)
with(cake3, table(recipe, temperature))
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake3)
(lsm7 <- ls_means(model))

# Using show_tests with options:
show_tests(lsm7, fractions = TRUE)
show_tests(lsm7, fractions = TRUE, names = FALSE)

# Missing diagonal:
cake4 <- cake
cake4$temperature <- factor(cake4$temperature, ordered=FALSE)
cake4 <- droplevels(subset(cake4, temperature %in% levels(cake4$temperature)[1:3]))
cake4 <- droplevels(subset(cake4, !((recipe == "A" & temperature == "175") |
                                      (recipe == "B" & temperature == "185") |
                                      (recipe == "C" & temperature == "195") )))
# str(cake4)
with(cake4, table(recipe, temperature))
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4)
ls_means(model)


########### Various contrasts codings:

model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake3,
              contrasts = list(recipe="contr.sum", temperature="contr.helmert"))
(lsm8 <- ls_means(model))
# show_tests(lsm7)
# show_tests(lsm8)
stopifnot(
  isTRUE(all.equal(lsm7, lsm8, check.attributes=FALSE, tolerance=TOL))
)

# ambient contrasts not contr.treatment:
options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake3)
(lsm9 <- ls_means(model))
options(contrasts = c("contr.treatment", "contr.poly"))
options("contrasts")
stopifnot(
  isTRUE(all.equal(lsm7, lsm9, check.attributes=FALSE, tolerance=TOL))
)


