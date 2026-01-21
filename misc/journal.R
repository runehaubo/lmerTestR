# misc.R

library(devtools)
# has_devel()
r2path <- "C:/Users/RCHR0155/OneDrive - Region Hovedstaden/github/lmerTestR/lmerTest" # "~/GitHub/lmerTestR/lmerTest"
# document(pkg=r2path, roclets = c("namespace", "rd"))
document(pkg=r2path)
load_all(r2path)
# system.time(check(r2path, cran = FALSE))
# check(r2path, cran = TRUE)
# install(r2path)
# library(lmerTest)
# sessionInfo()

#####################################################################
## Install development version of lme4:


# devtools::install_github("lme4/lme4",dependencies=TRUE, force=TRUE)
# library(lme4)
# sessionInfo()

# remove.packages("lme4")
# install.packages("lme4")

#####################################################################
## Run reverse dependency checks:

## Install the revdepcheck package:

## Restart R!
# install.packages("pak")
# pak::pkg_install("r-lib/revdepcheck")

## In a terminal (eg. in RStudio), do the following:
## 1. cd to the root of the package
## 2. Start R by writing "R" -> Enter
## 3. Run the following code:
##      library(revdepcheck)
##      revdep_check(dependencies = "Depends", num_workers = 4)
##      revdep_check(dependencies = c("Depends", "Imports"), num_workers = 4)
##      revdep_check(num_workers = 4)
##    This will run reverse dependency check only on packages that depends on 
##    pkg = "lmerTest". Extend to 'dependencies = c("Depends", "Imports")' to 
##    run on more packages or even 'revdep_check()' to run all.
## 4. In 'this' (the usual) R-session run 
##      r2path <- "~/GitHub/lmerTestR/lmerTest"
##      revdep_summary(pkg=r2path) 
##    to get a status on the reverse dependency checks, and, say,
##      revdep_details(pkg=r2path, revdep = "AssumpSure")
##    to get details on a revdep-package.
## 5. When the revdep_check() is done use
##      revdep_report(pkg=r2path)
##    to generate humanly readable .md files under lmerTest/revdep.


library(revdepcheck)
?revdep_check()

r2path <- "~/GitHub/lmerTestR/lmerTest"
revdep_summary(pkg=r2path)  
revdep_report(pkg=r2path)
revdep_details(pkg=r2path, revdep = "sosta")
revdepcheck::revdep_details(, "variancePartition")

#####################################################################
## 2026-01-11

##  --no-manual --no-build-vignettes
install.packages("NCC", dependencies = TRUE)

library(devtools)

r2path <- "~/GitHub/lmerTestR/lmerTest"
check(pkg="~/GitHub/lmerTestR/lmerTest_3.2-0.tar.gz")
check(pkg="~/GitHub/lmerTestR/pkgs_revdep/variancePartition", 
      build_args = "--no-build-vignettes",
      args="--no-manual --no-build-vignettes")

install.packages("variancePartition", dependencies = TRUE)
?BiocManager::install("variancePartition")

pak::pkg_install("variancePartition", dependencies = TRUE)

# Your NAMESPACE would do

importFrom(utils, packageVersion)
if (utils::packageVersion("lme4") >= "2.0-0")
  importFrom(lme4, forceNewMerMod)

# and your '.onLoad' would do:
  
ns <- parent.env(environment())
if (packageVersion("lme4") < "2.0-0")
  assign("forceNewMerMod", envir = ns, inherits = FALSE,
         function(object, reference) object)


# ##############################################
# ######## lmerModLmerTest_to_lmerMod
# ##############################################
# lmerModLmerTest_to_lmerMod <- function(object) {
#   ## This replaces a simple "as(object, "lmerMod")"
#   res <- as(object, "lmerMod")
#   if(!is.null(attr(object, "upper")))
#     attr(res, "upper") <- attr(object, "upper")
#   if(!is.null(attr(object, "reCovs")))
#     attr(res, "reCovs") <- attr(object, "reCovs")
#   res
# }



#####################################################################
## 2026-01-08
model <- lmer(Preference ~ sens2 + Homesize + Gender +
                (Gender+Homesize|Consumer), data=carrots)
ranova(model)
drop1(model)
step(model)

forceNewMerMod
isNewMerMod
lme4:::.anyStructured

#####################################################################
## 2026-01-07
library(reformulas)

model <- lmer(Depth ~ Picture + diag(TVset | Assessor), TVbo)
model
ranova(model, reduce.terms = FALSE)
ranova(model, reduce.terms = TRUE)

findbars(orig_rhs)
splitForm(orig_rhs)$reTrmClasses

reforms[[1]]
full_form

get_newforms(reforms[[1]], special = "diag", full_formula = full_form)

?splitForm(Depth ~ Picture + cs(1 + Assessor | TVset, hom=TRUE) + (1 | Assessor), 
          defaultTerm = "", )

###################
orig_form <- Depth ~ Picture + diag(TVset | Assessor) + (1 | Assessor)
splitForm(orig_form)
splitForm(orig_form)$reTrmFormulas
(specials <- splitForm(orig_form)$reTrmClasses)
.known_specials <- c("us", "cs", "diag", "ar1")
(has_specials <- any(!splitForm(orig_form)$reTrmClasses %in% c("", "us")))

(orig_rhs <- orig_form[[length(orig_form)]])
(orig_lhs <- orig_form[[2]])
# Reconstruct formula - needed for terms like (1 | g1 / g2):
(fe_rhs <- deparse2(nobars(orig_rhs)))
findbars_x(orig_form, debug = FALSE, default.special = "us",
           specials=c("us", "diag"))
findbars_x(orig_form, debug = FALSE, default.special = "us",
           specials=.known_specials)

(reforms <- lapply(findbars_x(orig_rhs, default.special = "us", 
                              specials = .known_specials), deparse2)) 

(re_rhs <- lapply(reforms, function(rf) paste0("(", rf, ")")))
(full_rhs <- paste(c(list(fe_rhs), re_rhs), collapse=" + "))
(full_rhs <- paste(c(list(fe_rhs), reforms), collapse=" + "))
(full_form <- update.formula(orig_form, paste0(". ~", full_rhs)))
(full_form <- as.formula(paste0(orig_lhs, "~", full_rhs)))

rm_complete_terms(reforms, full_form)

(reform <- reforms[[1]])

#####################
## Pruned code:

####
# Test reduction of (Days | Subject) to (1 | Subject):
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
ranova(fm1) # 2 df test

# This test can also be achieved with anova():
fm2 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
anova(fm1, fm2, refit=FALSE)

# Illustrate reduce.test argument:
# Test removal of (Days | Subject):
ranova(fm1, reduce.terms = FALSE) # 3 df test

# The likelihood ratio test statistic is in this case:
fm3 <- lm(Reaction ~ Days, sleepstudy)
2*c(logLik(fm1, REML=TRUE) - logLik(fm3, REML=TRUE)) # LRT

# anova() is not always able to perform the same tests as ranova(),
# for example:
anova(fm1, fm3, refit=FALSE) # compares REML with ML and should not be used
anova(fm1, fm3, refit=TRUE) # is a test of ML fits and not what we seek

# Also note that the lmer-fit needs to come first - not an lm-fit:
# anova(fm3, fm1) # does not work and gives an error

# ranova() may not generate all relevant test:
# For the following model ranova() indicates that we should not reduce
# (TVset | Assessor):
fm <- lmer(Coloursaturation ~ TVset * Picture + (TVset | Assessor), data=TVbo)
ranova(fm)
# However, a more appropriate model is:
fm2 <- lmer(Coloursaturation ~ TVset * Picture + (1 | TVset:Assessor), data=TVbo)
anova(fm, fm2, refit=FALSE)
# fm and fm2 has essentially the same fit to data but fm uses 5 parameters
# more than fm2.

####

orig_form <- Reaction ~ Days + (Days || Subject)
orig_form <- Depth ~ Picture + diag(TVset | Assessor, hom=TRUE) + (1 | Assessor)
orig_form <- Depth ~ Picture + diag(TVset | Assessor) + (1 | Assessor)
orig_form <- Depth ~ Picture + diag(1 | Assessor / TVset) + (1 | TVset)
orig_form <- Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE)


sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
model <- lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), data=sleepstudy)
ranova(model)
model <- lmer(Reaction ~ Days + diag(Days | Subject), data=sleepstudy)
ranova(model)
model <- lmer(Depth ~ Picture + diag(TVset | Assessor, hom=TRUE) + (1 | Assessor), 
              data=TVbo)
ranova(model)
model <- lmer(Depth ~ Picture + (0 +  TVset | Assessor), 
              data=TVbo)
ranova(model)
ranova(model, reduce.terms = FALSE)

model <- lmer(Depth ~ Picture + diag(TVset | Assessor) + (1 | Assessor), 
              data=TVbo)
ranova(model)
model <- lmer(Depth ~ Picture + diag(1 | Assessor / TVset) + (1 | TVset), 
              data=TVbo)
ranova(model)
model <- lmer(Depth ~ Picture + (1 | Assessor / TVset) + (1 | TVset), 
              data=TVbo)
ranova(model)
model <- lmer(Depth ~ Picture + us(1 | Assessor / TVset) + (1 | TVset), 
              data=TVbo)
ranova(model)
model <- lmer(Reaction ~ Days + (Days | Subject), data=sleepstudy)
ranova(model)
model <- lmer(Reaction ~ Days + us(Days | Subject), data=sleepstudy)
ranova(model)
model <- lmer(Reaction ~ Days + (Days || Subject), data=sleepstudy)
ranova(model)
ranova(model, reduce.terms = FALSE)

ranova(model)
summary(model)
slotNames(model)
model@call
str(model, max.level = 2)


(orig_rhs <- orig_form[[length(orig_form)]])
(orig_lhs <- orig_form[[2]])
# Reconstruct formula - needed for terms like (1 | g1 / g2):
(fe_rhs <- deparse2(nobars(orig_rhs)))
# findbars_x(orig_form, default.special = "us", specials=.known_specials,
#            , expand_doublevert_method = "split")

has_double_bar <- grepl("||", deparse2(orig_rhs), fixed = TRUE)
(reforms <- findbars_x(orig_rhs, default.special = "us", 
                       specials = .known_specials,
                       expand_doublevert_method = "split")) 
(reforms_chr <- lapply(reforms, deparse2))
# expandDoubleVerts()
# splitForm(orig_form)
# splitForm(orig_form)$reTrmFormulas
(specials <- vapply(reforms, function(ref) splitForm(ref)$reTrmClasses, 
                    character(1L)))
.known_specials <- c("us", "cs", "diag", "ar1")
(has_specials <- any(!specials %in% c("", "us")))
(has_unknown_specials <- has_specials && !all(specials %in% .known_specials))

# (re_rhs <- lapply(reforms, function(rf) paste0("(", rf, ")")))
# (full_rhs <- paste(c(list(fe_rhs), re_rhs), collapse=" + "))
(full_rhs <- paste(c(list(fe_rhs), reforms), collapse=" + "))
# (full_form <- update.formula(orig_form, paste0(". ~", full_rhs)))
(full_form <- as.formula(paste0(orig_lhs, "~", full_rhs)))

(res_forms <- rm_complete_terms(reforms_chr, full_form))

get_newforms()

res_forms[[2]]
lmer(res_forms[[2]], data=TVbo)

# (reform <- reforms[[1]])
rm_complete_terms <- function(terms, full_formula, random=TRUE) {
  # Remove random-effect formula terms from original model formula (full_formula)
  forms <- lapply(terms, function(reform) {
    form <- update.formula(full_formula, paste0("~.- ", reform))
    environment(form) <- environment(full_formula)
    form
  })
  names(forms) <- terms
  forms
}

#####################################################################
## 2026-01-07

?TVbo
data("TVbo")
all(complete.cases(TVbo))

fm1 <- lmer(Depth ~ Picture + (1|Assessor:TVset) + 
              (1 | Assessor) + (1 | TVset), TVbo)
fm1
summary(fm1)

fm2 <- lmer(Coloursaturation ~ Picture + us(1 + TVset | Assessor), TVbo)
summary(fm2)

fm3 <- lmer(Coloursaturation ~ Picture + cs(1 + TVset | Assessor, hom=TRUE), TVbo)
fm3

fm4 <- lmer(Coloursaturation ~ Picture + cs(1 + TVset | Assessor, hom=FALSE), TVbo)
summary(fm4)

fm5 <- lme4::lmer(Depth ~ Picture + cs(1 + Assessor | TVset, hom=TRUE) + (1 | Assessor), TVbo)
fm5
summary(fm5)

fm6 <- lmer(Coloursaturation ~ Picture + cs(1 + Assessor | TVset, hom=FALSE), TVbo)
summary(fm6)

library(data.table)
fm1
fm5
varcor <- VarCorr(fm1) |> as.data.frame() |> as.data.table()
vars <- varcor[1:2, vcov]
sum(vars) |> sqrt()


VarCorr(fm5) |> as.data.frame()
REMLcrit(fm1)
REMLcrit(fm5)

sqrt(0.69631660)
s1 <- 0.69631660
s2 <- 0.05747056

s1 / (s1 + s2) |> sqrt()
s2 / (s1 + s2)

######################################
## There appears to be an error here:
remotes::install_github("bbolker/reformulas")

library(lme4)
data("TVbo", package="lmerTest")
fm5 <- lmer(Depth ~ Picture + cs(Assessor | TVset, hom=TRUE) + (1 | Assessor), TVbo)
fm5

######################################

#####################################################################
## 2026-01-06
## Fix scoping issue in ranova():

library(lmerTest)
f <- function(data) {
  lmer(Petal.Length ~ Sepal.Length + (1|Species), data=data)
}

res <- f(iris)
res
ranova(res) # fails

data <- iris
ranova(res) # now it works

# A model with ||-notation:
f2 <- function(data) {
  lmer(Reaction ~ Days + (Days||Subject), data)
}
res <- f2(sleepstudy)
res
ranova(res)

# A model with multiple RE terms:
f3 <- function(data) {
  lmer(Coloursaturation ~ TVset*Picture + (1|Assessor:TVset) + (1|Assessor), data)
}
res <- f3(TVbo)
res
ranova(res)

fm <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor:TVset) + (1|Assessor),
           data=TVbo)
step(fm)
(an1 <- ranova(fm))

#####################################################################

library(reformulas)
?findbars
RHSForm(y ~ x + (1|g))
anySpecial(~diag(1))
anySpecial(~diag)
anySpecial(~diag[[1]])
anySpecial(~diag[1])
anySpecial(~s)
anySpecial(~s(hello+goodbye,whatever))
reformulas:::findReTrmClasses()


f <- y ~ 1 + x
RHSForm(f) <- quote(2+x^2)
print(f)
reOnly(~ 1 + x + y + (1|f) + (1|g))
addForm0(y~x,~1)
addForm0(~x,~y)
ff <- findbars_x(y~1+(x|f/g))
expandAllGrpVar(ff)
expandAllGrpVar(quote(1|(f/g)/h))
expandAllGrpVar(quote(1|f/g/h))
expandAllGrpVar(quote(1|f*g))
expandAllGrpVar(quote(1|f+g))
expandAllGrpVar(quote(a+b|f+g+h*i))
expandAllGrpVar(quote(s(log(d), k = 4)))
expandAllGrpVar(quote(s(log(d+1))))
splitForm(quote(us(x,n=2)))
findbars_x(~ 1 + (x + y || g), expand_doublevert_method = "diag_special")
findbars_x(~ 1 + (x + y || g), expand_doublevert_method = "split")
findbars_x(~ 1 + (1 | f) + (1 | g))
findbars_x(~ 1 + (1 | f) + (1 | g))
findbars_x(~ 1 + (1|h) + (x + y || g), expand_doublevert_method = "split")
findbars_x(~ 1 + (1|Subject))
findbars_x(~ (1||Subject))
findbars_x(~ (1|Subject) + diag(1 | g))
findbars_x(~ diag(1|Subject), default.special = NULL, specials = c("diag"))
findbars_x(~ 1 + x)
findbars_x(~ s(x, bs = "tp"))
findbars_x(y ~ a + log(b) + s(x, bs = "tp") + s(y, bs = "gp"),
           target = "s", default.special = NULL)
findbars(~ 1 + (x + y || g))
findbars(~ 1 + diag(1 | f) + (1 | g))
findbars_x(~ 1 + (1|h) + (x + y || g), expand_doublevert_method = "split")
findbars(~ 1 + (1|Subject))
findbars(~  1 + (1||Subject))
findbars(~ (1||Subject))
findbars(~ 1 + x)
inForm(z~.,quote(.))
inForm(z~y,quote(.))
inForm(z~a+b+c,quote(c))
inForm(z~a+b+(d+e),quote(c))
f <- ~ a + offset(x)
f2 <- z ~ a
inForm(f,quote(offset))
inForm(f2,quote(offset))
extractForm(~a+offset(b),quote(offset))
extractForm(~c,quote(offset))
extractForm(~a+offset(b)+offset(c),quote(offset))
extractForm(~offset(x),quote(offset))
dropHead(~a+offset(b),quote(offset))
dropHead(~a+poly(x+z,3)+offset(b),quote(offset))
drop.special(x~a + b+ offset(z))

replaceForm(quote(a(b+x*c(y,z))),quote(y),quote(R))

ss <- ~(1 | cask:batch) + (1 | batch)
replaceForm(ss,quote(cask:batch),quote(batch:cask))
replaceForm(ss, quote(`:`), quote(`%:%`))

replaceForm(quote(Depth ~ Picture + diag(1 | TVset:Assessor) + diag(1 | Assessor)), 
            quote(diag(1 | TVset:Assessor)), as.name(" "))

substitute(A ~ B + C, list(C = quote(D)))
?quote
parse(text=quote(D))

full_form


removeForm <- function(term, target) {
  if (identical(term, target)) 
    return(repl)
  if (!inForm(term, target)) 
    return(term)
  if (length(term) == 2) {
    return(substitute(OP(x), list(OP = replaceForm(term[[1]], target, repl), 
                                  x = replaceForm(term[[2]], target, repl))))
  }
  return(substitute(OP(x, y), list(OP = replaceForm(term[[1]], target, repl), 
                                   x = replaceForm(term[[2]], target, repl), 
                                   y = replaceForm(term[[3]], target, repl))))
  
}

formula(reform)
as.name(reform)
?quote

expandGrpVar(quote(x*y))
expandGrpVar(quote(x/y))
findReTrmClasses()

## (silly/impractical formula, for illustration only)
form <- mpg ~ 1 + (1|gear) + (factor(cyl)|gear) + (1 + hp | carb)
fr <- model.frame(subbars(form), data = mtcars)
findbars(form)
rterms <- mkReTrms(findbars(form), fr)
names(rterms)
## block sizes (latent variables per block) of each term
(nperblock <- lengths(rterms$cnms))
## latent variables per term
(nperterm <- diff(rterms$Gp))
with(rterms, identical(unname(nl*nperblock), nperterm))
## illustrate reordering of terms
dd <- expand.grid(a = 1:7, b = 1:3, c = 1:5, d = 1:9)
dd$y <- 1
form2 <- y ~ 1 + (1|a) + (1|b) + (1|c) + (1|d)
rterms2 <- mkReTrms(findbars(form2), dd, reorder.terms = TRUE)
## reorder elements into original formula order
with(rterms2, cnms[order(ord)])
## reorder splitForm output to match mkReTrms components
ss <- splitForm(form2)
ss$reTrmFormulas[rterms2$ord]

findbars_x(~ 1 + s(x) + (f|g) + diag(x|y))
no_specials(findbars_x(~ 1 + s(x) + (f|g) + diag(x|y)))
no_specials(~us(f|g))
no_specials(~us(f|g, extra_arg))

nobars(Reaction ~ Days + (Days|Subject))
nobars_(Reaction ~ Days + (Days|Subject))

f <- ~ 1 + a  + b + (1 | f) +  (0 + a | f) + (1 + a | g) + (a + b | h ) + (1 + a + b | i)
randint(f)


?splitForm(~x+y)                     ## no specials or RE
splitForm(~x+y+(f|g))               ## no specials
splitForm(~x+y+diag(f|g))           ## one special
splitForm(~x+y+(diag(f|g)))         ## 'hidden' special
splitForm(~x+y+(f|g)+cs(1|g))       ## combination
splitForm(~x+y+(1|f/g))             ## 'slash'; term
splitForm(~x+y+(1|f/g/h))             ## 'slash'; term
splitForm(~x+y+(1|(f/g)/h))             ## 'slash'; term
splitForm(~x+y+(f|g)+cs(1|g)+cs(a|b,stuff))  ## complex special
splitForm(~(((x+y))))               ## lots of parentheses
splitForm(~1+kr(f|g,n=2))
splitForm(~1+s(x, bs = "tp"))

noSpecials(y~1+us(1|f))
noSpecials(y~1+us(1|f),delete=FALSE)
noSpecials(y~us(1|f))
noSpecials(y~us(1|f), delete=FALSE)
noSpecials(y~us(1|f), debug=TRUE)
noSpecials(y~us+1)  ## should *not* delete unless head of a function
noSpecials(~us(1|f)+1)   ## should work on a one-sided formula!
noSpecials(~s(stuff) + a + b, specials = "s")
noSpecials(cbind(b1, 20-b1) ~ s(x, bs = "tp"))

sub_specials( ~ s(a, k=4))
sub_specials( ~ (1|x) + (a + b || y) + s(a, k=4))
sub_specials(Reaction ~ s(Days) + (1 + Subject))
sub_specials(~ s(cos((y^2*3)/2), bs = "tp"))

subbars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days + (Days + Subject)

#####################################################################
## 2026-01-03

########################
## Does something get stuck in the evaluation of devfun?
## This appears to happen if illegal parameter values are applied.
model <- lme4::lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
devfun <- update(model, devFunOnly=TRUE)
(optpar <- model@optinfo$val)
getME(model, "par")

attr(model, "upper")
model@lower
devfun(optpar)
optpar3 <- c(1, .2, 1.1)
devfun(optpar3) # NaN understandable
devfun(optpar) # NaN ??

########################


has_new_covs <- function(x) {
  covs <- lme4:::getReCovs(x)
  any(vapply(covs, function(x) !inherits(x, "Covariance.us"),
             FUN.VALUE = logical(1)))
}


fm0 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm1 <- lmer(Reaction ~ Days + diag(Days|Subject), sleepstudy)

has_new_covs(fm0)
has_new_covs(fm1)

lme4:::getReCovs(fm0)
lme4:::getReCovs(fm1)

?`Covariance-class`

getPar(fm0)


## Check out KR in comparison to Satterthwaite:
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
anova(fm, ddf="Sat")
anova(fm, ddf="Ken")

fm <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
anova(fm, ddf="Sat")
anova(fm, ddf="Ken")

fm <- lmer(Reaction ~ Days + diag(Days | Subject, hom = TRUE), sleepstudy)
anova(fm, ddf="Sat")
anova(fm, ddf="Ken")

fm <- lmer(Reaction ~ Days + ar1(Days | Subject), sleepstudy)
lme4:::getReCovs(fm)
getME(fm, "theta")
fm@theta


anova(fm, ddf="Sat")
an <- anova(fm, ddf="Ken")
show_tests(an)

L <- c(0, 1)
model <- fm
?getME

pbkrtest::KRmodcomp(model, L, betaH=0)$test
pbkrtest::SATmodcomp(model, L, betaH=0)
contestMD(model, drop(L), rhs=0, confint=FALSE)

library(pbkrtest)
(fm0 <- lmer(Reaction ~ diag(Days|Subject, hom=TRUE), sleepstudy))
(fm1 <- lmer(Reaction ~ Days + diag(Days|Subject, hom=TRUE), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + I(Days^2) + diag(Days|Subject, hom=TRUE), sleepstudy))

## Test for no effect of Days in fm1, i.e. test fm0 under fm1
SATmodcomp(fm1, "Days")
SATmodcomp(fm1, ~.-Days)
L1 <- cbind(0, 1) 
contestMD(fm1, L1)
## SATmodcomp(fm1, L1) ## FIXME
SATmodcomp(fm1, fm0)

KRmodcomp(fm1, "Days")
KRmodcomp(fm1, ~.-Days)
L1 <- cbind(0, 1) 
KRmodcomp(fm1, L1) ## FIXME
KRmodcomp(fm1, fm0)

## Could we walk through all the lme4models with specials and check the 
## ddf computed by KR and Sat? 
## Also check using another categorical re-model possibly with unbalanced data.





mk_Ftable(Fvalue=x["FtestU", "stat"], ndf=x[1L, "ndf"],
          ddf=x[1L, "ddf"], sigma=sigma(model),
          Fscale=x["Ftest", "F.scaling"])
#####################################################################
## 2026-01-02

fm <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
fm
summary(fm)
anova(fm, type="I")
anova(fm, ddf="lme4")

object@optinfo$val

# library(lme4)
sessionInfo()

sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm1.ar1 <- lme4::lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), 
                      sleepstudy)
fm1.ar1
fm1.ar1@optinfo$val
fm1.ar1@theta
getME(fm1.ar1, "theta")
lme4:::getReCovs(fm1.ar1)
lme4:::getPar(fm1.ar1)
devFun <- update(fm1.ar1, devFunOnly=TRUE)
devFun(fm1.ar1@optinfo$val) 
fm1.ar1@devcomp$cmp["REML"]
varpar <- getVarPar(fm1.ar1)
devfun_vp(varpar, devFun, reml=TRUE)


fm1.us <- lme4::lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
model <- fm1.us
model@optinfo$val
model@theta
devfun <- update(model, devFunOnly=TRUE)
devfun(model@optinfo$val)
varpar_opt <- getVarPar(model)
varpar <- getVarPar(model)

devfun_vp(varpar_opt, devFun, reml=TRUE)

as_lmerModLT(model, devFun)

fm1.cs <- lme4::lmer(Reaction ~ Days + cs(Days | Subject), 
                         sleepstudy)
fm1.cs.lt <- lmer(Reaction ~ Days + cs(Days | Subject), 
                  sleepstudy)
model <- fm1.cs
model <- fm1.ar1
model <- fm1.cs.lt
devFun <- update(model, devFunOnly=TRUE)
devFun(model@optinfo$val) 
model@devcomp$cmp["REML"]

varpar_opt <- unname(c(model@optinfo$val, sigma(model)))
devfun_vp(varpar_opt, devFun, reml=TRUE)

?lme4::getME()
?lme4::modular
?`Covariance-class`
methods(class = "merMod")

## Unstructured
fm1.us <- lme4::lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
## Diagional
fm1.diag <- lme4::lmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
fm1.diag.hom <- lme4::lmer(Reaction ~ Days + diag(Days | Subject, hom = TRUE), 
                     sleepstudy)
## Compound symmetry
fm1.cs <- lme4::lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
fm1.cs.hom <- lme4::lmer(Reaction ~ Days + cs(Days | Subject, hom = TRUE), 
                   sleepstudy)
## Auto-regressive order 1
sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm1.ar1 <- lme4::lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), 
                sleepstudy, REML = TRUE)

lme4models <- namedList(fm1.us,
                        fm1.diag, 
                        fm1.diag.hom,
                        fm1.cs, 
                        fm1.cs.hom, 
                        fm1.ar1)

model <- lme4models[[1]]
## Native devfun:
devfun <- update(model, devFunOnly=TRUE)
## Evaluate native devfun at optimum:
devfun(getOptPar(model))
## Check that devfun returns the same value as that saved in the model object:
stopifnot(
  all.equal(unname(getME(model, "devcomp")$cmp["REML"]), 
            devfun(getOptPar(model)), tolerance=1e-6) # TRUE
)
## Get varpar (including residual SD):
(varpar <- getVarPar(model))
## Evaluate devfun_vp at the optimum:
devfun_vp(varpar, devfun, reml=TRUE)
## Check that devfun_vp returns the same value as native devfun:
stopifnot(
  all.equal(unname(getME(model, "devcomp")$cmp["REML"]), 
            devfun(getOptPar(model))) # TRUE
)
## Here we also want to check that devfun and and devfun_vp returns the same 
## value at non-optimum values of varpar.
## Hmm. Because sigma is profiled out of devfun this cannot be done right away.
## Maybe we can optimize over the sigma param to get the right thing? 
## No, we would need to optimize over all parameters to get equivalence.
## Though we should be able to set one of the parameters in common to a 
## particular value and optimize the rest. 

model <- lme4models[[6]]
devfun <- update(model, devFunOnly=TRUE)
(optpar <- getOptPar(model))
(varpar <- getVarPar(model))
(optpar2 <- optpar * 1.1)
(varpar2 <- c(optpar2, varpar[length(varpar)]))

## Evaluate gradients to ensure that devfun and devfun_vp are both 
## functions with optima at optpar and varpar respectively:
library(numDeriv)
(g_devfun <- grad(devfun, optpar))
(g_devfun_vp <- grad(devfun_vp, varpar, devfun=devfun, reml=TRUE))
stopifnot(
  all(abs(g_devfun) < 1e-3),
  all(abs(g_devfun_vp) < 1e-3)
)
## These are not zero as expected:
(g_devfun <- grad(devfun, optpar2))
(g_devfun_vp <- grad(devfun_vp, varpar2, devfun=devfun, reml=TRUE)) 
## Try optimizing devfun_vp:
x <- nlminb(start=varpar2, objective = devfun_vp, devfun=devfun, reml=TRUE,
            control=list(trace=1))
## Check that the optimum is re-achieved:
stopifnot(
  all(abs(varpar - x$par) < 1e-4),
  abs(devfun_vp(varpar, devfun=devfun, reml=TRUE) - 
        devfun_vp(x$par, devfun=devfun, reml=TRUE)) < 1e-6
)

## Optimize devfun and devfun_vp over all but one of the parameters in turn to
## check that devfun and devfun_vp gives the same deviance and parameter values
## for settings away from the REML optimum. This is to build confidence that
## devfun_vp is a valid implementation of the deviance function from LLMs.
for(j in seq_along(optpar)) { # j <- 1
  ## Check that all parameters are within bounds:
  stopifnot(
    model@lower < optpar,
    optpar < attr(model, "upper"),
    model@lower < optpar2,
    optpar2 < attr(model, "upper")
  )
  ## Evaluate deviance function at optimum (for safety):
  devfun(optpar)
  ## Optimize devfun over all but the j'th parameter:
  (startpar <- optpar[-j])
  res <- nlminb(start=startpar, objective = function(p) {
    (Par <- optpar2)
    (Par[-j] <- p)
    devfun(Par) 
    }, control = list(trace=1),
    lower = model@lower[-j], 
    upper = attr(model, "upper")[-j])
  ## Evaluate devfun_vp:
  devfun_vp(varpar, devfun=devfun, reml=TRUE)
  ## Optimize devfun_vp over all but the j'th parameter:
  (startpar_vp <- varpar[-j])
  (np <- length(startpar_vp))
  res_vp <- nlminb(start=startpar_vp, objective = function(p) {
    (Par <- optpar2)
    (Par[-j] <- p[-np])
    Par <- c(Par, p[np])
    devfun_vp(Par, devfun=devfun, reml=TRUE) 
  }, control = list(trace=1),
  lower = c(model@lower[-j], 0), 
  upper = c(attr(model, "upper")[-j], Inf))
  ## Compare parameter estimates (except for sigma):
  res$objective - res_vp$objective
  res$par - res_vp$par[seq_along(res$par)]
  ## Check that parameter estimates and deviance values agree:
  stopifnot(
    abs(res$objective - res_vp$objective) < 1e-8,
    all(abs(res$par - res_vp$par[seq_along(res$par)]) < 1e-4)
  )
}

########################
## Does something get stuck in the evaluation of devfun?
## This appears to happen if illegal parameter values are applied.
model <- lme4::lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
devfun <- update(model, devFunOnly=TRUE)
(optpar <- model@optinfo$val)
attr(model, "upper")
model@lower
devfun(optpar)
optpar3 <- c(1, .2, 1.1)
devfun(optpar3) # NaN
devfun(optpar) # NaN ??

########################

nlminb(start=varpar[-(1:2)], objective = function(p) {
  devfun_vp(c(1, .02, p), devfun=devfun, reml=TRUE) })
nlminb(start=optpar[-(1:2)], objective = function(p) {
  devfun(c(1, .02, p)) })


hessian(devfun, optpar)
hessian(devfun_vp, varpar, devfun=devfun, reml=TRUE)

optimize(function(sig) -devfun_vp(c(allpar[, 1], sig), devfun = devfun, reml=TRUE), 
         interval=c(20, 30)) 

i <- 2
devfun(allpar[-length(varpar), i])
devfun_vp(allpar[, i], devfun, reml = TRUE)



model@optinfo$val
model@theta
getOptPar(model)
(varpar <- getVarPar(model))
devfun <- update(model, devFunOnly=TRUE)
devfun(getOptPar(model))
unname(model@devcomp$cmp["REML"]) - devfun(getOptPar(model)) # 0
all.equal(unname(model@devcomp$cmp["REML"]), 
          devfun(getOptPar(model))) # TRUE
devfun_vp(varpar, devfun, reml=TRUE)
unname(model@devcomp$cmp["REML"]) - devfun_vp(varpar, devfun, reml=TRUE) # 0
all.equal(unname(model@devcomp$cmp["REML"]), 
          devfun_vp(varpar, devfun, reml=TRUE)) # TRUE

fm_lt <- as_lmerModLT(model, devfun)
print(fm_lt)
summary(fm_lt)
summary(fm_lt, ddf="Ken")
anova(fm_lt)
anova(fm_lt, ddf = "Ken")
lme4::VarCorr(model)
VarCorr

lapply(lme4models, VarCorr)

fm <- lme4::lmer(Reaction ~ Daysf + diag(1 + Daysf | Subject, hom=TRUE), 
                 sleepstudy, REML = FALSE)
fm <- lme4::lmer(Reaction ~ Daysf + diag(1 | Subject, hom=TRUE), 
                 sleepstudy, REML = FALSE)
VarCorr(fm)
logLik(fm, REML=TRUE)
deviance(fm)

#####################################################################
## 2026-01-01

## Try out this version:
forceNewMerMod <-
  function (object, reference = object) {
    if (is.null(attr(object, "upper")))
      attr(object, "upper") <- getUpper(reference)
    if (is.null(attr(object, "reCovs")))
      attr(object, "reCovs") <- getReCovs(reference)
    object
  }


## From: ?`Covariance-class`
## Unstructured
fm1.us <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
fm1.us
## Diagional
fm1.diag <- lmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
fm1.diag
fm1.diag.hom <- lmer(Reaction ~ Days + diag(Days | Subject, hom = TRUE), 
                     sleepstudy)
fm1.diag.hom

## Compound symmetry
fm1.cs <- lme4::lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)

eigen(fm1.cs@optinfo$derivs$Hessian)


fm1.cs.hom <- lmer(Reaction ~ Days + cs(Days | Subject, hom = TRUE), 
                   sleepstudy)
fm1.cs.hom

## Auto-regressive order 1
sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm1.ar1 <- lme4::lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), 
                      sleepstudy, REML = FALSE)
fm1.ar1
devFun <- update(fm1.ar1, devFunOnly=TRUE)
devFun(fm1.ar1@theta) # Why does this not work?
# [1] NaN
# Warning message:
#   In sqrt(1 - rho2) : NaNs produced

## Works fine for other covariance-structures, for example:
fm1.cs.hom <- lme4::lmer(Reaction ~ Days + cs(Days | Subject, hom = TRUE), 
                         sleepstudy)
fm1.cs.hom
devFun <- update(fm1.cs.hom, devFunOnly=TRUE)
devFun(fm1.cs.hom@theta) 


model <- fm1.cs
devFun <- update(model, devFunOnly=TRUE)
devFun(model@theta)
model@devcomp$cmp["REML"]

varpar_opt <- unname(c(model@theta, sigma(model)))
devfun_vp(varpar_opt, devFun, reml=TRUE)
model@devcomp$cmp["REML"]

# Compute Hessian:
h <- numDeriv::hessian(func=devfun_vp, x=varpar_opt, devfun=devfun,
                       reml=is_reml)





m <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
bm <- as_lmerModLmerTest(m)
slotNames(bm)
attr(m, "reCovs")
attr(m, "upper")

lme4:::getReCovs(m)
lme4:::getReCovs(bm)
lme4:::getUpper(m)
lme4:::getUpper(bm)

attr(bm, "reCovs")
attr(bm, "upper")

fm0 <- lme4::lmer(Reaction ~ Days + diag(Days|Subject), sleepstudy)
fm0
(fm0_reCovs <- attr(fm0, "reCovs"))
attr(fm0, "reCovs") <- NULL
attr(fm0, "reCovs") # now NULL
fm0 # Fails (expected)
lme4:::getReCovs(fm0) # Fails - unexpected
lme4:::getUpper(fm0) # OK
forceNewMerMod(fm0) # Fails -> unexpected.

attr(fm0, "reCovs") <- fm0_reCovs
attr(fm0, "reCovs") # non-NULL
fm0 # Now works

lme4:::getReCovs(fm0)
devfun <- update(fm0, devfunonly)

as_lmerModLmerTest(fm0)

fm1 <- lmer(Reaction ~ Days + diag(Days|Subject), sleepstudy)
fm1  ## Fail - now OK :-)
summary(fm1) ## Fail - now OK :-)
summary(fm1, ddf="lme4")

summary(fm1, ddf="lme4")
lme4:::getReCovs(fm1)
lme4:::upReCovs(lme4:::getReCovs(fm1), fm1@theta)
getAnywhere("summary.merMod")
VarCorr(fm1)

attr(fm1, "reCovs")


#####################################################################
## 2025-12-30

#################################
## Look into ranova behavior with lme4-2.0-0
library(lmerTest)
library(lme4)
packageVersion("lmerTest")
packageVersion("lme4") == "2.0.0"

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
ranova(fm1, reduce.terms = FALSE)
fm1 <- lmer(Reaction ~ Days + us(Days|Subject), sleepstudy)
ranova(fm1, reduce.terms = FALSE)
fm1 <- lmer(Reaction ~ Days + cs(Days|Subject), sleepstudy)
fm1
ranova(fm1, reduce.terms = FALSE)
fm1 <- lmer(Reaction ~ Days + diag(Days|Subject), sleepstudy)
attr(fm1, "reCovs")
fm1  ## Fail - now ok
summary(fm1) ## Fail - now ok
ranova(fm1, reduce.terms = FALSE)
fm1 <- lme4::lmer(Reaction ~ Days + diag(Days|Subject), sleepstudy)
attr(fm1, "reCovs")
fm1 ## OK
summary(fm1) ## OK
ranova(fm1, reduce.terms = FALSE)
ranova(fm1, reduce.terms = TRUE)

fm7 <- lmer(Reaction ~ Days + ar1(Days|Subject), sleepstudy)
fm7
summary(fm7)
ranova(fm7)
str(fm7, max.level = 2)
fm7@reCovs
fm7@upper
fm7@beta
attr(fm7, "reCovs")
attr(fm7, "beta")

fm2 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
anova(fm1, fm2, refit=FALSE)

ranova(fm1, reduce.terms = FALSE)
logLik(fm1)
fm1 <- lmer(Reaction ~ Days + (Days||Subject), sleepstudy)
fm1
summary(fm1)


drop1(fm1)
ranova(fm1) # 2 df test

# Test reduction of (Days | Subject) to (1 | Subject):
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
ranova(fm1) # 2 df test

# This test can also be achieved with anova():
fm2 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
anova(fm1, fm2, refit=FALSE)

# Illustrate reduce.test argument:

# Test removal of (Days | Subject):
ranova(fm1, reduce.terms = FALSE) # 3 df test
# The likelihood ratio test statistic is in this case:
fm3 <- lm(Reaction ~ Days, sleepstudy)
2*c(logLik(fm1, REML=TRUE) - logLik(fm3, REML=TRUE)) # LRT

# anova() is not always able to perform the same tests as ranova(),
# for example:
anova(fm1, fm3, refit=FALSE) # compares REML with ML and should not be used
anova(fm1, fm3, refit=TRUE) # is a test of ML fits and not what we seek

# Also note that the lmer-fit needs to come first - not an lm-fit:
# anova(fm3, fm1) # does not work and gives an error

# ranova() may not generate all relevant test:
# For the following model ranova() indicates that we should not reduce
# (TVset | Assessor):
fm <- lmer(Coloursaturation ~ TVset * Picture + (TVset | Assessor), data=TVbo)
ranova(fm)
# However, a more appropriate model is:
fm2 <- lmer(Coloursaturation ~ TVset * Picture + (1 | TVset:Assessor), data=TVbo)
anova(fm, fm2, refit=FALSE)
# fm and fm2 has essentially the same fit to data but fm2 uses 5 parameters
# more than fm.


#####################################################################
## 2025-12-29

library(lmerTest)
# Test reduction of (Days | Subject) to (1 | Subject):
fm1 <- lmer(Reaction ~ Days + (Days||Subject), sleepstudy)
summary(fm1)
ranova(fm1) # 2 df test


# devtools::install_github("lme4/lme4",dependencies=TRUE)
library(lme4)
sessionInfo()

m <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

str(m, max.level = 2)
attr(m, "upper")
attr(m, "reCovs")
attr(m, "call")
attr_m <- attributes(m)
str(attr_m, max.level = 1)
bm <- as_lmerModLmerTest(m)
slotNames(bm)
slotNames(m)

ff <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, devFunOnly = TRUE)
ff

#####################################################################
## 23-Oct-2020

install.packages("panelr", dependencies = TRUE)
install.packages("lme4")
install.packages("gsDesign")
library(gsDesign)

nEvents(hr=.5, beta=.2, tbl = TRUE, alpha = .05, sided = 2)

LRPower(total.sample.size=160, effect.size=0.5)

(coef <- log(.5))^2
(d <- (qnorm(.975) + qnorm(.80))^2 / (.5 *.5 * coef^2)) # 66

0.02*5000

4636/2 # 2318 i hver arm
.02*2318 + 0.01*2318 # => ~70 cases

nEvents(n = 66, beta=.2, tbl = TRUE, alpha = .05, sided = 2, hr=0.5)

# What would the size of CI be if 3 events were observed in each arm?
d1 <- 40; d2 <- 40
(lambda <- d2/d1)
(theta <- log(lambda))
(se.theta <- sqrt(1/d1 + 1/d2))
(CI.theta <- theta + c(-1, 1) * 1.96 * se.theta)
(CI.lambda <- exp(CI.theta))

ci_rate_ratio <- function(d1, d2) {
  if(missing(d2)) d2 <- d1
  lambda <- d2/d1
  theta <- log(lambda)
  se.theta <- sqrt(1/d1 + 1/d2)
  CI.theta <- theta + c(-1, 1) * 1.96 * se.theta
  CI.lambda <- exp(CI.theta)
  data.frame(rate_ratio=lambda, lower=CI.lambda[1], upper=CI.lambda[2])
}
d1 <- d2 <- 50
ci_rate_ratio(d1, d2)

x <- seq(20, 160, 2)
X <- do.call(rbind, lapply(x, function(y) ci_rate_ratio(y, y)))
X$n <- x
library(ggplot2)
# dev.new()
ggplot(X, aes(n, upper)) + geom_line() +
  # geom_line(data=X, aes(x, lower)) +
  geom_hline(yintercept = 1)


ci_rate_ratio(33, 33)

ci_rate_ratio(47)
ci_rate_ratio(155)

# Risk reduction:
1 - ci_rate_ratio(42, 66-42)
ci_rate_ratio(24, 42)
1 - ci_rate_ratio(44, 22)
ci_rate_ratio(22, 44)
42/24
24/42
44+22
24 * .75


#####################################################################
## 19-Oct-2020

set.seed(101)
test <- data.frame(TM = factor(rep(rep(c("org","min"),each=3),3)),
                   dep = runif(18,0,20),
                   ind = runif(18,0,7),
                   dorp = factor(rep(1:3,each=6)))
# library(lmerTest)
full.model <- lmer(dep ~ TM + ind + (1 | dorp),  data=test)
x <- reduce_random(full.model)
str(x, max.level=1)
fm <- attr(x, "model")
fm

res <- step(full.model)
str(res, max.level = 1)
res$random
res$fixed

attr(res, "model")
attr(res, "drop1")

library(lmerTest)
multipar_test <- function(mod, pars = "all", ddf = "Satterthwaite"){
  cfs <- attributes(terms(mod))$term.labels
  pars <- gsub(" ", "", pars)
  cfs <- gsub(":", "*", cfs)
  if (length(pars) == 1){
    if (pars == "all"){
      pars <- cfs
    }
  }
  if (all(pars %in% cfs)){
    contrast_mat <- matrix(0, nrow = length(pars),
                           ncol = length(cfs) + 1)
    for (i in 1:length(pars)){
      contrast_mat[i, which(cfs %in% pars[i]) + 1] <- 1
    }
    outpt <- lmerTest::contestMD(mod, contrast_mat, ddf = ddf)
    return(outpt)
  }
  stop("Error: pars do not match model pars")
}

fm <- lmer(Informed.liking ~ Product*Information + (1|Consumer) , data=ham)
an <- anova(fm)
show_tests(an)

contestMD(fm, L=diag(length(fixef(fm))))

contestMD(fm, L=diag(x = c(0, rep(1, length(fixef(fm)) - 1))))
multipar_test(fm)

data("cake", package="lme4")
?ham
cake$temperature <- factor(cake$temperature, ordered=FALSE)
with(cake, table(recipe, temperature))
fm1 <- lmer(angle ~ temperature + (1|recipe:replicate), cake)

contestMD(fm1, L=diag(length(fixef(fm1))))



multipar_test(fm1)

library("devtools");
library(lmerTest)
install.packages("lme4")
install_github("lme4/lme4")
update.packages(ask = FALSE, checkBuilt = TRUE)
remove.packages(pkgs = "lme4")
sessionInfo()

library(lmerTest)
data("Orthodont", package = "nlme")
Orthodont[,"age"] <- Orthodont[,"age"] - 11 ## center the time covariate
model <- lmer(distance ~ Sex * age + (age | Subject), data = Orthodont)
anova(fit)

######

data("cake", package="lme4")
cake4 <- cake
cake4$temperature <- factor(cake4$temperature, ordered=FALSE)
cake4 <- droplevels(subset(cake4, temperature %in% levels(cake4$temperature)[1:3]))
cake4 <- droplevels(subset(cake4, !((recipe == "A" & temperature == "175") |
                                      (recipe == "B" & temperature == "185") |
                                      (recipe == "C" & temperature == "195") )))
str(cake4)
with(cake4, table(recipe, temperature))
fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake4)
(an1 <- anova(fm1))
model <- fm1
XX <- lmerTest:::get_rdX(fm1)

#####################################################################
## 08-Apr-2020
sessionInfo()
?devtools::document

install.packages("Rcpp")
install.packages(c("pbkrtest", "lme4", "numDeriv", "ggplot2"))
update.packages(checkBuilt = TRUE, ask=FALSE)

update.packages(ask = FALSE, checkBuilt = TRUE)  # update R packages
remove.packages("tinytex")
install.packages("tinytex")

#####################################################################
test <- data.frame(TM = as.factor(c(rep("org", 3), rep("min", 3),rep("org", 3), rep("min", 3),rep("org", 3), rep("min", 3))),
                   dep = runif(18,0,20),
                   ind = runif(18,0,7),
                   dorp = as.factor(c(rep(1,6),rep(2,6),rep(3,6))))
test

library(lmerTest)
full.model <- lme4::lmer(dep ~ TM + ind + (1 | dorp),  data=test)  #lmerTest:: give same outcome
ranova(full.model)

step.model <- lmerTest::step(full.model)   # n=16
step.model<- lmerTest::step(full.model, direction="both",k=log(16))   # n=16
traceback()

summary(step.model)

full.model <- lmer(dep ~ TM + ind + (1 | dorp),  data=test)  #lmerTest:: give same outcome
ranova(full.model)

step.model <- step(full.model)   # n=16

traceback()

BIC(step.model)
#####################################################################
data("TVbo", package="lmerTest")
library(lme4)

fm <- lmerTest::lmer(Coloursaturation ~ TVset + Picture + (1|Assessor:TVset) + (1|Assessor),
           data=TVbo)
# summary(fm)
class(fm)
fm2 <- update(fm, ~.-TVset)
class(fm2)
# summary(fm2)
class(as(fm2, "lmerMod"))

sessionInfo()

#####################################################################
# install.packages("lmerTest") # current CRAN version
data(soup, package="ordinal")
soup$response <- as.numeric(as.character(soup$SURENESS))
fm <- lmer(response ~ PRODID * DAY + (1 | RESP), data=soup)
anova(fm)

##############
# Load model fits fm1 and fm2 generated with lmerTest version 2.3-37:
load(system.file("testdata","legacy_fits.RData", package="lmerTest"))

# Apply some methods defined by lmerTest:
anova(fm1)
summary(fm1)
contest(fm1, c(0, 1))
contest(fm1, c(0, 1), joint=FALSE)
drop1(fm1)
ranova(fm1)
step(fm1, ddf="Ken")

# lme4-methods also work:
fixef(fm1)

# Ditto for second model fit:
anova(fm2)
summary(fm2)
ls_means(fm2)
difflsmeans(fm2)
step(fm2)

fm2b <- fm2
class(fm2b) <- "lmerMod"
fm2c <- as_lmerModLmerTest(fm2b)
getCall(fm2c)
summary(fm2c)
drop1(fm2c)
ranova(fm2c)
step(fm2c)

object <- fm2c
alpha.random=0.1
alpha.fixed=0.05
reduce.fixed=TRUE
reduce.random=TRUE
keep <- character(0L)
ddf=c("Satterthwaite")



##############

?anova.merModLmerTest
?summary.merModLmerTest
?legacy


fm1_new <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
(an1_new <- anova(fm1))
(sfm1_new <- summary(fm1))

fm2_new <- lmer(Informed.liking ~ Product*Information*Gender
            + (1|Product:Consumer) , data=ham)
(an2_new <- anova(fm2_new))
(sfm2_new <- summary(fm2_new))
step(fm2_new)

#########

library("lmerTest")
packageVersion("lmerTest")
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
(an1 <- anova(fm1))
(sfm1 <- summary(fm1))

fm2 <- lmer(Informed.liking ~ Product + Information + Gender +
              (1|Product:Consumer) , data=ham)
(an2 <- anova(fm2))
(sfm2 <- summary(fm2))
step(fm2)

save(fm1, an1, sfm1, fm2, an2, sfm2,
     file="~/GitHub/lmerTestR/package/inst/testdata/legacy_fits.RData")

load(file="~/GitHub/lmerTestR/package/inst/testdata/legacy_fits.RData")

data("sleepstudy", package="lme4")
m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
anova(m) # with p-values from F-tests using Satterthwaite's denominator df


?summary.lmerModLmerTest
?anova.lmerModLmerTest
install.packages("pbkrtest")
library(pbkrtest)
#####################################################################

m <- lmer(Reaction ~ Days + I(Days^2) + (Days | Subject), sleepstudy)
(an_3 <- anova(m, type=3))
(an_KR <- anova(m, ddf="Kenward-Roger"))

show_tests(an_3)

library(pbkrtest)
m <- lme4::lmer(Reaction ~ Days + I(Days^2) + (Days | Subject), sleepstudy)
KRmodcomp(m, matrix(c(0, 1, 0), nrow=1))
sessionInfo()


update.packages("Matrix")

nrow(sleepstudy)
samp <- sort(sample(180, 170))
d <- sleepstudy[samp, ]

models <- lapply(1:10, function(i) {
  samp <- sort(sample(180, 170))
  d <- sleepstudy[samp, ]
  lmer(Reaction ~ Days + (Days|Subject), data = d)
})

str(models, max.level = 1)
for(m in models) print(anova(m))

library(lmerTest)
fm <- lmer(Informed.liking ~ Gender + Information * Product + (1 | Consumer) +
             (1 | Consumer:Product), data=ham)
lmerTest::summary(fm)
## Error: 'summary' is not an exported object from 'namespace:lmerTest'
summary(fm)
summary
packageVersion("lmerTest")

lmerTest:::summary.lmerModLmerTest(fm)
sessionInfo()

getNamespaceExports("lmerTest")

#####################################################################
data("sleepstudy", package="lme4")
myfit <- function(formula, data) {
  lme4::lmer(formula = formula, data = data)
}
fm2 <- myfit(Reaction ~ Days + (Days|Subject), sleepstudy)
as_lmerModLmerTest(fm2)


#####################################################################

install.packages("lmerTest") # current CRAN version
library("lmerTest")
packageVersion("lmerTest")
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
anova(fm1)
# Analysis of Variance Table of type III  with  Satterthwaite
# approximation for degrees of freedom
#      Sum Sq Mean Sq NumDF DenDF F.value    Pr(>F)
# Days  30031   30031     1    17  45.853 3.264e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
getwd()
save(fm1, file = "../misc/cranversion.rda")
# restart R
devtools::install_github("runehaubo/lmerTestR")
library("lmerTest")
load("../misc/cranversion.rda")
fm0 <- fm1
# fm1 ## no print method
anova(fm1) # Error
class(fm1) <- "lmerMod"
fm1 ## print.merMod
anova(fm1) # lme4::anova.merMod
emmeans(fm1, "Days") # works

fm2 <- as_lmerModLmerTest(fm1)
anova(fm2)   # works
summary(fm2)  # works
emmeans(fm2, "Days")  # works



lmerTest:::anova.lmerModLmerTest
anova.merModLmerTest <- function(object, ..., type = c("III", "II", "I", "3", "2", "1"),
                                 ddf = c("Satterthwaite", "Kenward-Roger", "lme4")) {
  class(object) <- "lmerMod"
  dots <- list(...)
  models <- if (length(dots))
    sapply(dots, is, "lmerModLmerTest") | sapply(dots, is,
                                                 "merMod") | sapply(dots, is, "lm")
  else logical(0)
  if(any(models)) return(NextMethod())
  df <- match.arg(ddf)
  if (df == "lme4")
    return(anova(object, ...))

  object <- as_lmerModLmerTest(object)
  anova(object, type=type, ddf=ddf)
}
anova(fm0)

library(afex)
load(url("http://singmann.org/download/r/m_machines_lmerTest-pre3.0.rda"))
class(m_machines)
emmeans(m_machines, "Machine", lmer.df = "asymptotic")

anova(m_machines)


get_anova <- function(m) anova(m)
get_anova(fm0)

#####################################################################

data("sleepstudy", package="lme4")
myfit <- function(formula, data) {
  fm2 <- lme4::lmer(formula = formula, data = data)
  fm2
}
fm2 <- myfit(Reaction ~ Days + (Days|Subject), sleepstudy)
getCall(fm2)
environment(myfit)
environment(formula(fm2))
environment(fm2)
eval(getCall(fm2))
fm3 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
eval(getCall(fm3))


library("lme4")

data("sleepstudy", package="lme4")
fm <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
anova(lmerTest::as_lmerModLmerTest(fm))
lmerTest_anova(fm)

lmerTest_summary(fm)

myfun <- function(formula, data) {
  fm <- lme4::lmer(formula = formula, data = data)
  fm2 <- lmerTest::as_lmerModLmerTest(fm)
  anova(fm2)
}

myfun(Reaction ~ Days + (Days|Subject), sleepstudy)

myfun2 <- function(formula, data) {
  fm <- lme4::lmer(formula = formula, data = data)
  lmerTest_anova(fm)
}

myfun2(Reaction ~ Days + (Days|Subject), sleepstudy)




myfit <- function(formula, data) {
  fm <- lme4::lmer(formula = formula, data = data)
  lmerTest::as_lmerModLmerTest(fm)
}

fm2 <- myfit(Reaction ~ Days + (Days|Subject), sleepstudy)
getCall(fm2)
lmerTest::as_lmerModLmerTest(fm2)

myfit2 <- function(formula, data) {
  lmerTest::lmer(formula = formula, data = data)
}
fm22 <- myfit2(Reaction ~ Days + (Days|Subject), sleepstudy)
class(fm22)

#####################################################################

library("lmerTest")

# Test evaluation of update inside a function:
myupdate <- function(m, ...) {
  update(m, ...)
}

data("sleepstudy", package="lme4")
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
tmp <- sleepstudy
rm(sleepstudy)
fmA <- update(fm1, data = tmp) # works
fmB <- myupdate(fm1, data = tmp) # also works
# Same except for 'call':
fmB@call <- fmA@call
stopifnot(isTRUE(all.equal(fmA, fmB)))




library("lmerTest")
sessionInfo()

myupdate <- function(m, ...) {
  update(m, ...)
}



myupdateB <- function(m, ...) {
  # browser()
  m <- update(m, ...)
  print(getCall(m))
  lmerTest::as_lmerModLmerTest(m)
}


data("sleepstudy", package="lme4")
(fm1 <- lmerTest::lmer(Reaction ~ Days + (Days|Subject), sleepstudy))


devf <- update(fm2, devFunOnly=TRUE)
devf
names(formals(devf)[1]) == "theta"

getCall(fm1)
tmp <- sleepstudy
rm(sleepstudy)
update(fm1, data = tmp) # works
myupdate(fm1, data = tmp) # fails
update(fm1, data = sleepstudy) # works
myupdate(fm1, data = sleepstudy) # fails

str(fm1, max.level = 1)

fm2 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
getCall(fm2)
update(fm2, data = sleepstudy) # works
myupdate(fm2, data = sleepstudy) # works
myupdateB(fm2)
myupdateB(fm2, data=sleepstudy)

myupdate2 <- function(m, ...) {
  eval.parent(update(m, ...))
}

myupdate2(fm1, data = sleepstudy) # fails

m <- fm1

myupdate3 <- function(m, data=NULL) {
  mc <- as.list(getCall(m))
  mc <- if(!is.null(data)) mc$data <- data
  Call <- as.call(mc)
  eval.parent(Call)
}

myupdate3(fm1)
myupdate3(fm1, data=sleepstudy)
tmp <- sleepstudy

myupdate3(fm1, data=tmp)

myupdate4 <- function(m, ...) {
  m <- update(m, ...)
  m
}

myupdate4(fm2)
myupdate4(fm2, data=sleepstudy)


#####################################################################

library(lmerTest)
system.time(fm1 <- lme4::lmer(Preference ~ sens1 + sens2 + (Age + Homesize + Income)^2 +
                                (1 + sens1 + sens2 | Consumer) + (1 | Product),
                              data=carrots)) # time to fit model
system.time(fm2 <- as_lmerModLmerTest(fm1)) # time to compute derivatives
system.time(an <- anova(fm2, type = "2"))  # time to compute tests


# Perhaps I need to check that the right subset, weights etc. are included
# correctly in the model.frame.
data("sleepstudy", package="lme4")
fm <- lme4::lmer(Reaction ~ poly(Days, 2) + (0 + poly(Days, 2) | Subject), sleepstudy)
fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
head(model.frame(fm))
?na.action(fm)
head(recover_data(getCall(fm), terms(fm), na.action(fm)))
.all.vars <- emmeans:::.all.vars
.reformulate <- emmeans:::.reformulate

get_all_vars


#####################################################################
data("sleepstudy", package="lme4")
fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
fm

tmp <- sleepstudy
m <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
m

args <- as.list(getCall(m))[-1]
foo <- do.call(lmerTest::lmer, args)
foo

args <- c(as.list(getCall(m))[-1], devFunOnly=TRUE)
foo <- do.call(lme4::lmer, args, envir = parent.frame())

foo <- as.call(list(quote(lmerTest::lmer), formula=m[[1]],
                    data=imageMat, ...))

mf <- model.frame(m)
args <- c(as.list(getCall(m)), devFunOnly=TRUE)
args <- as.list(getCall(m))
args$data <- mf
Call <- as.call(c(list(quote(lme4::lmer)), args[-1]))
eval.parent(Call)

# foo <- do.call(lmerTest::lmer, list(formula = m[[1]], data=imageMat, ...))

#####################################################################
# library(lmerTest)

data("sleepstudy", package="lme4")
test <- function() {
  tmp <- sleepstudy
  m <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  summary(m)
}
test()

test <- function() {
  tmp <- sleepstudy
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  summary(m)
}
test()

test <- function() {
  tmp <- sleepstudy
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  if(requireNamespace("lmerTest", quietly = TRUE)) {
    summary(lmerTest::as_lmerModLmerTest(m))
  }
}
test()

sessionInfo()

#######

rm(tmp)
data("sleepstudy", package="lme4")
test <- function() {
  tmp <- sleepstudy
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  if(requireNamespace("lmerTest", quietly = TRUE)) {
    summary(lmerTest::as_lmerModLmerTest(m))
  }
}

# error
test()

######
data("sleepstudy", package="lme4")
test <- function() {
  tmp <- sleepstudy
  m <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), data = tmp)
  summary(m)
}

test()
sessionInfo()




#####################################################################
data("cake", package="lme4")
cake$Temp <- factor(cake$temperature, ordered = FALSE)

fm1 <- lmer(angle ~ recipe * Temp + (1|recipe:replicate), cake)
fm2 <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
terms(fm2)
containment(fm1)
containment(fm2)

m <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
containment(m)
old_containment(m)

object <- m
terms <- terms(object)
data_classes <- attr(terms(object, fixed.only=FALSE), "dataClasses")
# Note: need fixed.only for merMod objects to get dataClasses
term_names <- attr(terms, "term.labels")
factor_mat <- attr(terms, "factors")

relatives(data_classes, "temp", term_names, factor_mat)

object <- fm1
terms <- terms(object)
data_classes <- attr(terms(object, fixed.only=FALSE), "dataClasses")
# Note: need fixed.only for merMod objects to get dataClasses
term_names <- attr(terms, "term.labels")
factor_mat <- attr(terms, "factors")
# setNames(lapply(term_names, function(term) {
#   term_names[relatives(data_classes, term, term_names, factor_mat)]
# }), term_names)

term_names[term_contain("recipe", factor_mat, data_classes)]

(factors <- factor_mat)
(classes.term <- data_classes)

(term1 <- term <- term_names[1])
(term2 <- term_names[3])
is.relative(term, term2)
factors[, term]
factors[, term2]
all(!(factors[, term1] & (!factors[, term2])))
!(factors[, term1] & (!factors[, term2]))
factors[, term1] & (!factors[, term2])

# relatives <- function(classes.term, term, term_names, factors) {
  ## checks if the terms have the same number of covariates (if any)
checkCovContain <- function(term1, term2) {
  num.numeric <- which(classes.term=="numeric")
  num.numeric.term1 <- which((num.numeric %in% which(factors[,term1]!=0))==TRUE)
  num.numeric.term2 <- which((num.numeric %in% which(factors[,term2]!=0))==TRUE)

  if((length(num.numeric.term1)>0 && length(num.numeric.term2)>0)||
     (length(num.numeric.term1)==0 && length(num.numeric.term2)==0))
    return(all(num.numeric.term2 == num.numeric.term1))
  else
    return(FALSE)
}
is.relative <- function(term1, term2) {
  all(!(factors[, term1] & (!factors[, term2]))) && checkCovContain(term1, term2)
}
if(length(term_names) == 1) return(NULL)
which.term <- which(term == term_names)
(1:length(term_names))[-which.term][sapply(term_names[-which.term],
                                           function(term2) is.relative(term, term2))]
# }

# Containment:
# F1 (A) is contained in F2 (A:B) [F2 (A:B) contains F1 (A)] if
# - F1 and F2 involve the same X's (if any)
# - F2 involve more factors (say B) than F1
# - All factors in F1 (if any) are part of F2

# Check if all variables in F1 are also in F2:
get_vars <- function(term, factors) rownames(factors)[factors[, term] == 1]

(var_names <- rownames(factors))

vars_F1 <- get_vars(term, factors)
vars_F2 <- get_vars(term2, factors)
all(vars_F1 %in% vars_F2) # all variables in F1 are also in F2
length(setdiff(vars_F2, vars_F1)) > 0L # F2 involve more variables than F1

data_classes[vars_F1]
data_classes[vars_F2]

numerics_F1 <- vars_F1[match(data_classes[vars_F1], "numeric", 0L)]
numerics_F2 <- vars_F2[match(data_classes[vars_F2], "numeric", 0L)]
setequal(numerics_F1, numerics_F2)


contain <- function(F1, F2, factors, dataClasses) {
  get_vars <- function(term, factors) rownames(factors)[factors[, term] == 1]
  vars_F1 <- get_vars(F1, factors)
  vars_F2 <- get_vars(F2, factors)
  numerics_F1 <- vars_F1[match(data_classes[vars_F1], "numeric", 0L)]
  numerics_F2 <- vars_F2[match(data_classes[vars_F2], "numeric", 0L)]
  all(vars_F1 %in% vars_F2) && # all variables in F1 are also in F2
    length(setdiff(vars_F2, vars_F1)) > 0L && # F2 involve more variables than F1
    setequal(numerics_F1, numerics_F2) # F1 and F2 involve the same covariates (if any)
}

contain(term, "Temp", factors, data_classes)
contain(term, "recipe:Temp", factors, data_classes)

term_names

sapply(term_names, function(term_nm)
  contain(term, term_nm, factors, data_classes))

term_contain <- function(term, term_names, factors, dataClasses) {
  contain <- function(F1, F2) {
    get_vars <- function(term, factors) rownames(factors)[factors[, term] == 1]
    vars_F1 <- get_vars(F1, factors)
    vars_F2 <- get_vars(F2, factors)
    numerics_F1 <- vars_F1[match(data_classes[vars_F1], "numeric", 0L)]
    numerics_F2 <- vars_F2[match(data_classes[vars_F2], "numeric", 0L)]
    all(vars_F1 %in% vars_F2) && # all variables in F1 are also in F2
      length(setdiff(vars_F2, vars_F1)) > 0L && # F2 involve more variables than F1
      setequal(numerics_F1, numerics_F2) # F1 and F2 involve the same covariates (if any)
  }
  sapply(term_names, function(term_nm) contain(term, term_nm))
}

dataClasses <- data_classes

term_contain <- function(object, ...) UseMethod("term_contain")

term_contain.default <- function(object, factors, dataClasses) {
  # Compute logical vector indicating if term is contained in each of the terms
  # in the model (term_names)
  get_vars <- function(term) {
    # Extract vector of names of all variables in a term
    # term: name of a term
    # factors: attr(terms_object, "factors")
    rownames(factors)[factors[, term] == 1]
  }
  contain <- function(F1, F2) {
    # Returns TRUE if F1 is contained in F2 (i.e. if F2 contains F1)
    # F1, F2: Names of terms, i.e. attr(terms_object, "term.labels")
    all(vars[[F1]] %in% vars[[F2]]) && # all variables in F1 are also in F2
      length(setdiff(vars[[F2]], vars[[F1]])) > 0L && # F2 involve more variables than F1
      setequal(numerics[[F1]], numerics[[F2]]) # F1 and F2 involve the same covariates (if any)
  }
  term_names <- colnames(factors)
  # Get (named) list of all variables in terms:
  vars <- lapply(setNames(term_names, term_names), get_vars)
  # Get (named) list of all _numeric_ variables in all terms:
  numerics <- lapply(vars, function(varnms)
    varnms[match(dataClasses[varnms], "numeric", 0L)])
  # Check if object=term is contained in each model term:
  sapply(term_names, function(term_nm) contain(object, term_nm))
}

term_contain.terms <- function(object, dataClasses) {
  term_names <- attr(object, "term.labels")
  factor_mat <- attr(object, "factors")
  lapply(setNames(term_names, term_names), term_contain.default,
         factors=factor_mat, dataClasses=dataClasses)
}

class(terms)
term_contain("recipe", factor_mat, data_classes)
term_contain(terms, data_classes)



## A + B + C + A:B
##
## A:B - The conventional test of A:B (last) - R(A:B | A, B, C)
## C - Adjusted for A, B and A:B - R(C | A, B, AB)
## A - Adjusted for B and C, but not A:B - R(A | B, C)
## B - Adjusted for A and C, but not A:B - R(B | A, C)

## If a term ('T'; e.g. C or A:B) IS NOT CONTAINED in any other term, move the
## columns in X corresponding to T last, and obtain the Type I contrast for T
## using doolittle.
##
## If a term ('T'; e.g. A or B) IS CONTAINED in another term (such as A:B),
## remove the columns in X corresponding to the containing terms, move the
## columns in X corresponding to T last, and obtain the Type I contrast for T
## using doolittle.
##
## A Type II test is a test for a term T that adjusts for all terms that do not
## contain T.

## Containment:
## F1 (A) is contained in F2 (A:B) [F2 (A:B) contains F1 (A)] if
## - F1 and F2 involve the same X's (if any)
## - F2 involve more factors (say B) than F1
## - All factors in F1 (if any) are part of F2

## A + X + A:X
## - A:X is not contained in any term
## - A is not contained in A:X since they do not involve the same X's
## - X is contained in A:X
##
## R(A:X | A, X)
## R(A | X, A:X)
## R(X | A)

## A + B + C + X + A:B + A:X
##
## R(A | B, C, X)
## R(B | A, C)
## R(C | A, B, AB)
## R(X | A, B, C, A:B) - not A:X
## R(A:B | A, B, C, X, A:X)
## R(A:X | A, B, C, X, A:X)
##


## For each term
## -


#####################################################################
library(ggplot2)
data("cake", package="lme4")
cake$Temp <- factor(cake$temperature, ordered = FALSE)

m1 <- lmer(angle ~ recipe + temp + (1|recipe:replicate), cake)
m2 <- lmer(angle ~ recipe + temp + Temp + (1|recipe:replicate), cake)

as.data.frame(VarCorr(m1))
as.data.frame(VarCorr(m2))

X1 <- model.matrix(m1)
X2 <- model.matrix(m2)

ncol(X1)
ncol(X2)

(LL <- zapsmall(get_Ldiffmat(X1, X2)))
(LL <- get_Ldiffmat2(X1, X2))
nicify_cols(LL)
contestMD(m2, LL)
drop1(m2)

fixef(m2)


model <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
model <- lmer(angle ~ recipe * Temp + (1|recipe:replicate), cake)
model <- lmer(angle ~ recipe + Temp + (1|recipe:replicate), cake)
head(as.data.frame(ls_means(model)))
head(lsm <- difflsmeans(model))
head(as.data.frame(difflsmeans(model)))
tail(as.data.frame(difflsmeans(model)))
plot(ls_means(model))
plot(difflsmeans(model), which=c("recipe", "Temp"))
tail(lsm)
x <- lsm

########

#' ## load lmerTest package
#' library(lmerTest)
#'
#' ## Fit linear mixed model to the ham data:
#' fm <- lmer(Informed.liking ~ Gender + Information * Product + (1 | Consumer) +
#'              (1 | Consumer:Product), data=ham)
#'
#' ## Summary including coefficient table with p-values for t-statistics using
#' ## Satterthwaite's method for denominator degrees of freedom:
#' summary(fm)
#'
#' ## Type III anova table with p-values for F-tests based on Satterthwaite's
#' ## method:
#' (aov <- anova(fm))
#'
#' ## Inspect the contrast matrix for the Type III test of Product:
#' show_tests(aov, fractions = TRUE)$Product
#'
#' ## Choose type II anova table with Kenward-Roger method for the F-test:
#' \dontrun{
#' if(requireNamespace("pbkrtest", quietly = TRUE))
#'   anova(fm, type=2, ddf="Kenward-Roger")
#' }
#'
#' ## Anova-like table of random-effect terms using likelihood ratio tests:
#' ranova(fm)
#'
#' ## F-tests of 'single term deletions' for all marginal terms:
#' drop1(fm)
#'
#' ## Least-Square means and pairwise differences:
#' (lsm <- ls_means(fm))
#' ls_means(fm, which = "Product", pairwise = TRUE)
#'
#' ## ls_means also have plot and as.data.frame methods:
#' \dontrun{
#' plot(lsm, which=c("Product", "Information"))
#' as.data.frame(lsm)
#' ## Inspect the LS-means contrasts:
#' show_tests(lsm, fractions=TRUE)$Product
#' }
#'
#' ## Contrast test (contest) using a custom contrast:
#' ## Here we make the 2-df joint test of the main effects of Gender and Information
#' (L <- diag(length(fixef(fm)))[2:3, ])
#' contest(fm, L = L)
#'
#' ## backward elimination of non-significant effects:
#' step_result <- step(fm)
#'
#' ## Elimination tables for random- and fixed-effect terms:
#' step_result
#'
#' # Extract the model that step found:
#' final_model <- get_model(step_result)

#####################################################################

prop <- structure(list(ArcSin = c(0.767445588347834, 0.523598775598299,
                          0.698443658152744, 0.777637319623078, 0.955316618124509, 0.507098504392337,
                          0.982059733068725, 0.973655486815118, 0.680739483658771, 0.63613206280501,
                          0.69048963272509, 0.553691179473802, 1.02887367459657, 1.08499756903507,
                          0.890333459528817, 1.06236790283103, 1.00958333472834, 0.478839982231845,
                          1.25400123962673, 0.771035875105878, 1.05476152712072, 0.765177018212134,
                          0.756281180193236, 1.20522277649318, 0.935491702817486, 0.700211555534291,
                          0.47524724767499, 0.68239850277138, 0.830915552416156, 0.655835164072781,
                          0.854584252090002, 0.858473568847598, 0.803807230236046, 0.276571799805924,
                          1.5707963267949, 1.09201210800869, 0.841856088093851, 1.5707963267949,
                          0.879845286975491, 1.5707963267949, 1.5707963267949, 1.06122504401673,
                          0.84986232200633, 1.14580754383249, 0.958753378292407, 0.817237859525958,
                          1.10714871779409), Part = c("First Block", "First Block", "First Block",
                                                      "First Block", "First Block", "First Block", "First Block", "First Block",
                                                      "First Block", "First Block", "First Block", "First Block", "Last Block",
                                                      "Last Block", "Last Block", "Last Block", "Last Block", "Last Block",
                                                      "Last Block", "Last Block", "Last Block", "Last Block", "Last Block",
                                                      "Last Block", "First Block", "First Block", "First Block", "First Block",
                                                      "First Block", "First Block", "First Block", "First Block", "First Block",
                                                      "First Block", "First Block", "First Block", "Last Block", "Last Block",
                                                      "Last Block", "Last Block", "Last Block", "Last Block", "Last Block",
                                                      "Last Block", "Last Block", "Last Block", "Last Block"), Participant = structure(c(1L,
                                                                                                                                         1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                                                                         1L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
                                                                                                                                         3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L), .Label = c("1",
                                                                                                                                                                                                             "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                                                                                                                                                                                             "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
                                                                                                                                                                                                             "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
                                                                                                                                                                                                             "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46",
                                                                                                                                                                                                             "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
                                                                                                                                                                                                             "58", "59", "60"), class = "factor"), Condition = c("Label",
                                                                                                                                                                                                                                                                 "Label", "Label", "Label", "Label", "Label", "Label", "Label",
                                                                                                                                                                                                                                                                 "Label", "Label", "Label", "Label", "Label", "Label", "Label",
                                                                                                                                                                                                                                                                 "Label", "Label", "Label", "Label", "Label", "Label", "Label",
                                                                                                                                                                                                                                                                 "Label", "Label", "NoLabel", "NoLabel", "NoLabel", "NoLabel",
                                                                                                                                                                                                                                                                 "NoLabel", "NoLabel", "NoLabel", "NoLabel", "NoLabel", "NoLabel",
                                                                                                                                                                                                                                                                 "NoLabel", "NoLabel", "NoLabel", "NoLabel", "NoLabel", "NoLabel",
                                                                                                                                                                                                                                                                 "NoLabel", "NoLabel", "NoLabel", "NoLabel", "NoLabel", "NoLabel",
                                                                                                                                                                                                                                                                 "NoLabel")), .Names = c("ArcSin", "Part", "Participant", "Condition"
                                                                                                                                                                                                                                                                 ), row.names = c(NA, -47L), class = c("grouped_df", "tbl_df",
                                                                                                                                                                                                                                                                                                       "tbl", "data.frame"), vars = "Participant", drop = TRUE, indices = list(
                                                                                                                                                                                                                                                                                                         0:23, 24:46), group_sizes = c(24L, 23L), biggest_group_size = 24L, labels = structure(list(
                                                                                                                                                                                                                                                                                                           Participant = structure(c(1L, 3L), .Label = c("1", "2", "3",
                                                                                                                                                                                                                                                                                                                                                         "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                                                                                                                                                                                                                                                                                                                                                         "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
                                                                                                                                                                                                                                                                                                                                                         "25", "26", "27", "28", "29", "30", "31", "32", "33", "34",
                                                                                                                                                                                                                                                                                                                                                         "35", "36", "37", "38", "39", "40", "41", "42", "43", "44",
                                                                                                                                                                                                                                                                                                                                                         "45", "46", "47", "48", "49", "50", "51", "52", "53", "54",
                                                                                                                                                                                                                                                                                                                                                         "55", "56", "57", "58", "59", "60"), class = "factor")), row.names = c(NA,
library(lmerTest)                                                                                                                                                                                                                                                                                                                                                                                                                                -2L), class = "data.frame", vars = "Participant", drop = TRUE, .Names = "Participant"))
prop.lmer.0 <- lmer(ArcSin ~ Part*Condition + (1 | Participant), data = prop)
summary(prop.lmer.0)
anova(prop.lmer.0, type=1)
drop1(prop.lmer.0)

#####################################################################
# Moving 'hard-coded' tests from old lmerTest to new lmerTest

model <- lmer(Preference ~ (Age + Homesize + Income)^3 +
                (1 + sens2 | Consumer), data=carrots)
ranova(model)
drop1(model) # FIX this
anova(model, type=1)
anova(model, type=2)
anova(model, type=3)
anova(model, type="yates")
anova(model, type="marginal")

show_tests(anova(model, type=3))


model <- lmer(Preference ~ (Age + Homesize + Income)^2 +
                (1 + sens2 | Consumer), data=carrots)
ranova(model)
drop1(model)
anova(model, type=1)
anova(model, type=2)
anova(model, type=3)
anova(model, type="yates")
anova(model, type="marginal")

summary(m.carrots)


#####################################################################


fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
class(fm)

contest1D(fm, c(0, 1))
contest1D(fm, c(0, 1), ddf="Ken")

get_KR1D(c(0, 1), fm)

model <- fm

data("soup", package="ordinal")
str(soup)
soup$sure <- as.numeric(soup$SURENESS)

model <- lme4::lmer(sure ~ PRODID * DAY + (1 | RESP), data=soup)
model@call

fm <- lmer(sure ~ PRODID * DAY + (1 | RESP), data=soup)


## update doesn't quite works as it should:
update(fm, control=lmerControl(calc.derivs = TRUE), evaluate=FALSE)
update(fm, control=lmerControl(calc.derivs = TRUE), evaluate=TRUE)


library(ordinal)
m0 <- clm(rating ~ temp * contact, data = wine)
m1 <- clm(rating ~ temp * contact, data = wine, doFit = FALSE)
names(m1) # model environment including functions to compute log-likelihood etc.
# ls.str(m1) # more information
m1$clm.nll(m1) # evaluate log-likelihood at starting values
m1$clm.nll(m1, par=coef(m0)) # evaluate (negative) log-likelihood at MLE
-logLik(m0) # same as above
m1$clm.grad(m1) # gradient vector
m1$clm.hess(m1) # Hessian matrix



## Consider a function 'means' which computes the marginal means
## in addition, define a function, 'effects' which computes the
## marginal means with the overall mean subtracted - potentially subtracting
## the mean/effect for all terms which it containes? Will this be similar to the
## effects-function for aov objects?
#####################################################################
sessionInfo()
Sys.getenv("R_LIBS_USER")
.libPaths()

model <- lmer(Reaction ~ 0 + (Days|Subject), sleepstudy)
terms(fm)
library(ordinal)
str(wine)

data("wine", package="ordinal")

fm1 <- lmer(response ~ temp + contact + (1|judge), data=wine)
anova(fm1)
lmerTest::anova(fm1)
lmerTest::summary.lmerModLmerTest(fm1)

emmeans(fm1, ~ temp)

fm1 <- lme4::lmer(response ~ temp + contact + (1|judge), data=wine)
emmeans(fm1, ~ temp, options=list(df="Satterthwaite"))
emmeans(fm1, ~ temp, options=list(dfargs="Satterthwaite"))

emmeans(fm1, ~ temp, options=list(disable.pbkrtest=TRUE))
?emmeans

library(emmeans)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

grepl("merModLmerTest", class(fm1))
grepl("merModLmerTest", "lmerModLmerTest")

anova_lmerTest <- function(object, ...) {
  # Dispatch the right anova method across lmerTest versions
  if(is_lmerTest_class(object) && requireNamespace("lmerTest", quietly = TRUE)) {
    if(packageVersion("lmerTest") < "2.0.37.90012")
      return(lmerTest::anova(object, ...)) else return(anova(object, ...))
  } else if(inherits(object, "merMod") && requireNamespace("lmerTest", quietly = TRUE)) {
    if(packageVersion("lmerTest") < "2.0.37.90012")
      return(lmerTest::anova(as(object, "merModLmerTest"), ...)) else
        return(anova(lmerTest::as_lmerModLmerTest(object), ...))
  } # Default:
  anova(object, ...)
}
is_lmerTest_class <- function(object)
  # Check if an object is of class merModLmerTest or lmerModLmerTest
  # Bridges across versions of lmerTest
  any(grepl("merModLmerTest", class(object)))
summary_lmerTest <- function(object, ...) {
  # Dispatch the right summary method across lmerTest versions
  if(is_lmerTest_class(object) && requireNamespace("lmerTest", quietly = TRUE)) {
    if(packageVersion("lmerTest") < "2.0.37.90012")
      return(lmerTest::summary(object, ...)) else return(summary(object, ...))
  } else if(inherits(object, "merMod") && requireNamespace("lmerTest", quietly = TRUE)) {
    if(packageVersion("lmerTest") < "2.0.37.90012")
      return(lmerTest::summary(as(object, "merModLmerTest"), ...)) else
        return(summary(lmerTest::as_lmerModLmerTest(object), ...))
  } # Default:
  summary(object, ...)
}


# anova_lmerTest <- function(object, ...) {
#   if(!inherits(object, "merMod")) anova(object, ...) else
#     if(requireNamespace("lmerTest") && packageVersion("lmerTest") < "2.0.37.90012")
#       lmerTest::anova(as(object, "merModLmerTest"), ...) else
#         anova(lmerTest::as_lmerModLmerTest(object), ...)
# }
# summary_lmerTest <- function(object, ...) {
#   if(!inherits(object, "merMod")) summary(object, ...) else
#     if(requireNamespace("lmerTest") && packageVersion("lmerTest") < "2.0.37.90012")
#       summary(as(object, "merModLmerTest"), ...) else
#         summary(lmerTest::as_lmerModLmerTest(object), ...)
# }

library(lme4)
requireNamespace("lmerTest")
sessionInfo()

fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
summary(fm1) # lme4 summary
lme4::summary(fm1) # Error: 'summary' is not an exported object from 'namespace:lme4'
summary(as(fm1, "lmerModLmerTest")) # Fake result
summary(as_lmerModLmerTest(fm1)) # could not find function "as_lmerModLmerTest"
summary(lmerTest::as_lmerModLmerTest(fm1)) # lmerTest summary

lmerTest::summary.lmerModLmerTest(fm1) # Error: 'summary.lmerModLmerTest' is not an exported object from 'namespace:lmerTest'
lmerTest::summary(fm1) # Error: 'summary' is not an exported object from 'namespace:lmerTest'

summary_lmerTest(fm1) # lmerTest summary
summary(as(fm1, "merModLmerTest"))

anova(fm1) # lme4 anova
lme4::anova(fm1) # Error: 'anova' is not an exported object from 'namespace:lme4'
anova(as(fm1, "lmerModLmerTest")) # Error
lmerTest::anova(as(fm1, "merModLmerTest")) # lmerTest anova
lmerTest::anova(fm1) # lme4 anova
anova(as_lmerModLmerTest(fm1)) # could not find function "as_lmerModLmerTest"
anova(lmerTest::as_lmerModLmerTest(fm1)) # lmerTest anova

lmerTest::anova.lmerModLmerTest(fm1) # Error: 'anova.lmerModLmerTest' is not an exported object from 'namespace:lmerTest'
lmerTest::anova.merModLmerTest(fm1) # Error: 'anova.lmerModLmerTest' is not an exported object from 'namespace:lmerTest'
lmerTest::anova(fm1) # Error: 'anova' is not an exported object from 'namespace:lmerTest'


fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
anova_lmerTest(fm1) # lmerTest anova
fm2 <- lmerTest::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
anova_lmerTest(fm2) # lmerTest anova

fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
summary_lmerTest(fm1) # lmerTest anova
fm2 <- lmerTest::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
summary_lmerTest(fm2) # lmerTest anova


summary(fm1)
anova(fm1)
anova.lmerModLmerTest(fm1)
anova.lmerModLmerTest(as_lmerModLmerTest(fm1))

anova_lmerTest(fm1)
summary_lmerTest(fm1)
calcSatterth(fm1, c(0, 1))

fm2 <- lm(Reaction ~ Days, sleepstudy)
summary.lmerModLmerTest(fm2)
summary(fm2)
anova.lmerModLmerTest(fm2)
anova(fm2)

anova_lmerTest(fm2)
summary_lmerTest(fm2)
calcSatterth(fm2, c(0, 1))

(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
summary.lmerModLmerTest(gm1)
summary(gm1)
anova.lmerModLmerTest(gm1)
anova(gm1)

anova_lmerTest(gm1)
summary_lmerTest(gm1)
calcSatterth(gm1, c(0, 1))

fm3 <- glm(cbind(incidence, size - incidence) ~ period,
           data = cbpp, family = binomial)
summary.lmerModLmerTest(fm3)
summary(fm3)
anova.lmerModLmerTest(fm3)
anova(fm3)

anova_lmerTest(fm3)
summary_lmerTest(fm3)
calcSatterth(fm3, c(0, 1))

library(hamlet)
library(lmerTest)
summary(fm1)
lmerTest::summary(fm1)
showMethods("summary")
methods("summary")

lme4::summary(fm1)

#####################################################################
library(lmerTest)
sessionInfo()
library(afex)

### replicate results from Table 15.2 to 15.6 (Maxwell & Delaney, 2004, pp. 774)
data(md_15.1)
### ANOVA results (Table 15.2)
aov_4(iq ~ timecat + (timecat|id),data=md_15.1, anova_table=list(correction = "none"))
### Table 15.3 (random intercept only)
# we need to set the base level on the last level:
contrasts(md_15.1$timecat) <- contr.treatment(4, base = 4)
# "Type 3 Tests of Fixed Effects"
(t15.3 <- mixed(iq ~ timecat + (1|id),data=md_15.1, check.contrasts=FALSE))
# "Solution for Fixed Effects" and "Covariance Parameter Estimates"
summary(t15.3$full.model)
### make Figure 15.2
plot(NULL, NULL, ylim = c(80, 140), xlim = c(30, 48), ylab = "iq", xlab = "time")

library(hamlet)
# install.packages("hamlet")

# Use the VCaP ARN data as an example
data(vcaplong)
arn <- vcaplong[vcaplong[,"Group"] == "Vehicle" | vcaplong[,"Group"] == "ARN",]
# lme4 is required for mixed-effects models
library(lme4)
# Fit an example fixed effects model
fit <- lmer(PSA ~ 1 + DrugWeek + ARN:DrugWeek + (1|ID) + (0 + DrugWeek|ID), data = arn)
# For reproducibility, set a seed
set.seed(123)
# Run a brief power analysis with only a few selected N values and a limited number of bootstrapping
# Balance strata over the ARN and non-ARN (=Vehicle) so that both contain equal count of individuals
power <- mem.powersimu(fit, N=c(4, 7, 10), boot=20, strata="ARN", plot=TRUE)
# Power curves are plotted, along with returning the power matrix at:
power
# Notice that each column corresponds to a specified fixed effects at the formula
# "1 + DrugWeek + ARN:DrugWeek"


s <- lmerTest::summary(fit)
anova(fit)
lmerTest::anova(fit)
sessionInfo()

#####################################################################
data(ham, package="lmerTest")
str(ham)
summary(ham)

m <- lmer(Informed.liking ~ Product + Information
          + (1|Consumer) + (1|Consumer:Product) , data=ham)
# anova table with p-values with Satterthwaite's approximation for denominator
# degrees of freedom
anova(m)


drop1(fm)
ranova(fm)

data("carrots")
str(carrots)
fm <- lmer(Preference ~ sens2 + Homesize + (1 + sens2 | Consumer), data=carrots)
anova(fm)
carr <- within(carrots, {
  Frequency <- factor(Frequency)
  Gender <- factor(Gender)
  Work <- factor(Work)
  Income <- factor(Income)
  Bitter <- as.integer(as.character(carrots$BITTER))
  Product <- product
  product <- BITTER <- NULL
})
str(carr)
nm <- c("Consumer", "Frequency", "Gender", "Age", "Homesize", "Work",
        "Income", "Preference", "Sweetness", "Bitter", "Crisp", "sens1",
        "sens2", "Product")
str(carr[nm])
carrots <- carr[nm]
getwd()
save(carrots, file="data/carrots.rda")

library(tools)
?checkRdaFiles()

files <- paste0("R/", list.files("R/"))
file <- files[2]
for(file in files) {
  fileLines <- readLines(file)
  # ind <- grep('lmerTestR', fileLines, fixed=TRUE)
  # fileLines[ind]
  tmp <- gsub('lmerTestR', 'lmerTest', fileLines, fixed=TRUE)
  # tmp[ind]
  writeLines(tmp, con=file)
}

files <- paste0("tests/", list.files("tests/"))
for(file in files) {
  fileLines <- readLines(file)
  ind <- grep("3.3.3", fileLines, fixed=TRUE)
  print(fileLines[ind])
  tmp <- gsub('"3.3.3"', '"3.3.2"', fileLines, fixed=TRUE)
  tmp[ind]
  writeLines(tmp, con=file)
}


fileLines <- readLines("R/lsmeans.R")
head(fileLines)

ind <- grep('"KR"', fileLines, fixed=TRUE)
fileLines[ind]

tmp <- gsub('"KR"', '"Kenward-Roger"', fileLines, fixed=TRUE)
tmp[ind]
?readLines
writeLines(tmp, con="R/lsmeans.R")
read

# Fit a model to the ham dataset:
data("ham", package="lmerTest")
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

# Extract final model:
final <- get_model(step_res)
# Proceed with anova table and model summary:
anova(final)
summary(final)


data("TVbo", package="lmerTest")
str(TVbo)
fm <- lmer(Coloursaturation ~ TVset + Picture + (1|Assessor:TVset) + (1|Assessor),
           data=TVbo)
ranova(fm)
drop1(fm)

(an1 <- ranova(fm))
(step_res <- step(fm))
reduce_fixed(fm, keep=c("Pic", "TV"))
model <- fm

fm <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor:TVset) + (1|Assessor) +
             (1 | Assessor:Picture) + (1 | Assessor:Picture:TVset),
           data=TVbo)
terms(fm)
model <- fm
(step_res <- step(fm))

fm <- lmer(Preference ~ sens2 + gender + (1 | product) +
             (sens2 | Consumer),
           data=carrots)
step(fm)

set.seed(1)
g <- sample(c("1", "2"), nrow(TVbo), replace=TRUE)
fm <- lmer(Coloursaturation ~ TVset*Picture + (1|g), data=TVbo)
(step_res <- step(fm))

model <- fm

flm <- lm(Coloursaturation ~ TVset * Picture, data=TVbo)
drop1(flm, test="F")
args(stats::step)
res <- step(flm, trace=0)

stopifnot(
  inherits(res, "lm")
)

#####################################################################
Rank <- function(X) qr(X)$rank
A <- cbind(c(1, 1, 1, 1),
           c(2, 1, 1, 2),
           c(1, 2, 3, 4))
Rank(A) # = 3, OK.
A0 <- A[, c(1, 3)]
get_Ldiffmat2(A0, A)

R <- resid(.lm.fit(x=A0, y=A))
keep <- order(colSums(abs(R)), decreasing = TRUE)[Rank(A) - Rank(A0)]
keep <- order(colSums(abs(R)), decreasing = TRUE)[ncol(A) - ncol(A0)]
Rh <- R[, keep, drop=FALSE]
zapsmall(t(A) %*% R)
Lt <- crossprod(A, R)
Lt[] <- zapsmall(qr.Q(qr(Lt))) # orthonormalize contrasts
t(Lt)

P <- t(A) %*% R

A %*% P
zapsmall(t(A) %*% R)
zapsmall(MASS::ginv(A) %*% R)

solve(crossprod(A)) %*% t(A) %*% R

t(A0) %*% R
t(A) %*% R

#####################################################################
X <- cbind(c(1, 1, 1, 1, 1),
           c(2, 1, 1, 2, 1),
           c(1, 2, 3, 4, 1))
(X0 <- X[, 1L, drop=FALSE])
(X0 <- X[, c(1, 3)])
Rank(X); Rank(X0); Rank(cbind(X0, X))
(Lt <- zapsmall(nullspace(cbind(X0, X))))
zapsmall(nullspace(t(cbind(X0, X))))

zapsmall(t(X) %*% zapsmall(nullspace(cbind(X0, X))))
zapsmall(t(X0) %*% zapsmall(nullspace(cbind(X0, X))))


t(X) %*% zapsmall(nullspace(t(cbind(X0, X))))
t(X0) %*% zapsmall(nullspace(t(cbind(X0, X))))

B <- rbind(c(4/5, -2/5),
           c(-2/5, 1/5))
Pc <- rbind(c(1/5, 2/5),
            c(2/5, 4/5))
B %*% Pc %*% c(3, 2)
fractions(B %*% c(3, 2))
#####################################################################
library(MASS)

proj_col <- function(X)
  # Compute the matrix that projects onto the columns space of X
  X %*% MASS::ginv(X)
orthonorm <- function(X)
  # Orthonormalize the columns of X
  qr.Q(qr(X))
iszero <- function(x, tol=1e-8) all(abs(x) < tol)
norm <- function(x) sqrt(sum(x^2)) # function to compute norm/length of a vector
nicify_cols <- function(X, standardize=TRUE, tol=1e-8) {
  # Only operate on non-zero columns:
  zero_col <- apply(X, 2, function(x) all(abs(x) < 1e-8))
  x <- X[, !zero_col, drop=FALSE]
  #
  x[, 1L] <- c(1L, integer(nrow(X) - 1L))
  x <- gram_schmidt(x, normalize = FALSE)
  if(standardize) {
    for(j in 1:min(dim(x))) {
      x[, j] <- x[, j] / ifelse(abs(x[j, j]) < tol, 1, x[j, j])
    }
  }
  # Sweep-out non-unit diagonals:
  # d <- diag(x)
  # d[abs(d) < 1e-8] <- 1
  # x <- sweep(x, 2, d, FUN = "/")
  # if(ncol(x) > 1L) for(j in 2:ncol(x)) {
  #   x[, j] <- x[, j] - proj(x[, j-1L], x[, j])
  #   x[, j] <- x[, j] / x[j, j] #  norm(x[, j])
  # }
  X[, !zero_col] <- x
  X
}
# nicify_cols(L)

gram_schmidt <- function(A, normalize=TRUE) {
  if((nc <- ncol(A)) > 1) for(j in 2:nc) for(i in seq_len(j-1))
    A[, j] <- A[, j] - proj(A[, i], A[, j])
  if(normalize) sweep(A, 2, apply(A, 2, norm), FUN="/") else A
}
proj <- function(e, a, tol=1e-8) {# Projection of a onto e
  e * sum(e*a) / sum(e^2)
}
proj <- function(e, a, tol=1e-8) {# Projection of a onto e
  stopifnot(length(e) == length(a))
  if(abs(el <- sum(e^2)) < tol) rep(0, length(e)) else e * sum(e*a) / sum(e^2)
}


X <- cbind(c(1, 1, 1, 1),
           c(2, 1, 1, 2),
           c(1, 2, 3, 4))
(X0 <- X[, 1L, drop=FALSE])
(X0 <- X[, c(1, 3)])
(X0 <- cbind(X[, 1], X[, 2] + X[, 3]))

(P <- proj_col(X))
(P0 <- proj_col(X0))
nL <- ncol(X) - ncol(X0) # assuming full rank of X and X0

(Po <- (diag(4) - P0) %*% P)
(Po <- P - P0)

(L <- L_all[seq_len(Rank(L_all)), , drop=FALSE])
# select nL linearly independent rows from L:
(L <- t(qr.Q(qr(t(L_all)))[, nL, drop=FALSE]))
(L <- L_all[seq_len(nL), , drop=FALSE])
zapsmall(L)

L %*% t(X)
X %*% t(L)

orthonorm(t(L))

t(X) %*% zapsmall((P - P0) %*% X)
zapsmall(orthonorm(t(X) %*% zapsmall((P - P0) %*% X)))

# Could we write a small function that makes a contrast look nice.
# Instead of normalizing the (row or col?) vectors, it starts with a 1
# and orthogonalized the rest onto that in a sequential manor.

# Get full and null design matrices:
X <- cbind(c(1, 1, 1, 1),
           c(2, 1, 1, 2),
           c(1, 2, 3, 4))
(X0 <- X[, 1L, drop=FALSE])
# Get Projection matrix onto col-space of X:
(P <- proj_col(X))
# Get Projection matrix onto col-space of X0:
(P0 <- proj_col(X0))
# Get Projection matrix onto orthogonal complement of X0 in X:
(Po <- P - P0)
# Compute full L-matrix:
(L_all <- zapsmall(Po %*% X))
(L_all <- R <- resid(fm <- .lm.fit(x=X0, y=X)))
# nL is the number of linearly independent rows of L:
(nL <- ncol(X) - ncol(X0)) # assuming full rank of X and X0
# Select linearly independent rows of L:
(L <- L_all[seq_len(nL), , drop=FALSE])
(L <- L_all[seq_len(Rank(L_all)), , drop=FALSE])

X %*% t(L)
Rank(L)
R
names(fm)
fm$pivot


# Orthogonalize rows of R:
(Ro <- t(gram_schmidt(t(R), normalize = FALSE)))
(Ro <- zapsmall(t(qr.Q(qr(t(R))))))
# this ensures that the first two rows are linearly independent (also orthogonal)
# Select the first nL rows of Ro:
(L <- Ro[seq_len(nL), , drop=FALSE])
(L <- Ro[seq_len(Rank(Ro)), , drop=FALSE])
# Bring L on a more humanly readable form:
(LL <- nicify_cols(L))

gram_schmidt(L, norm=FALSE)
zapsmall(t(orthonorm(t(L))))
X %*% t(LL)

(K <- (diag(4) - P0) %*% P %*% X)
(K <- (diag(4) - P0) %*% X)
(K <- (P - P0) %*% X)
(set <- with(qr(t(K))[-1], pivot[seq_len(rank)]))
(L <- K[set, , drop=FALSE])

t(K) %*% X
t(K) %*% X0
zapsmall(ginv(K) %*% X)
(bhat <- zapsmall(ginv(X) %*% K))
(bhat <- zapsmall(ginv(X) %*% K[, 2:3]))
X %*% bhat
bhat <- zapsmall(ginv(X) %*% X0)

P %*% (diag(4) - P0) %*% X
(diag(4) - P0) %*% X
(P - P0) %*% X

P0 %*% P %*% X
P %*% P0 %*% X
P0 %*% X

P %*% t(P)

(diag(4) - P) %*% t(P)

P %*% X
P0 %*% X
X0
P0 %*% P
(diag(4) - P0) %*% P %*% X
(diag(4) - P0) %*% X
X - P0 %*% X
(P - P0) %*% X

X <- cbind(c(1, 1, 1, 1),
           c(2, 1, 3, 2),
           c(1, 2, 3, 4),
           c(1, 4, 9, 14))
Rank(X)
(X0 <- X[, c(1, 3)])
(X0 <- cbind(X[, 1], X[, 1] + X[, 3]))
(X0 <- cbind(X[, 2], X[, 1] + X[, 3]))
(X0 <- X[, 1:2])
(X0 <- X[, 3:4])
(X0 <- cbind(3*X[, 2], -2*X[, 1] + 4*X[, 2]))
(X0 <- cbind(X[, 3], -2*X[, 3] + 4*X[, 4]))
(P <- proj_col(X))
(P0 <- proj_col(X0))
(Po <- P - P0)
(L_all <- zapsmall(Po %*% X))
Rank(L_all)
(nL <- ncol(X) - ncol(X0)) # assuming full rank of X and X0
(L <- L_save <- L_all[seq_len(nL), , drop=FALSE])

Rank(X)
Rank(X0)


get_Ldiffmat2(X0, X)
Rank(get_Ldiffmat2(X0, X))
(L <- zapsmall(get_Ldiffmat(X0, X)))

nicify_cols(L)
zapsmall(t(orthonorm(t(L))))
Rank(t(orthonorm(t(L))))

X
X0
zapsmall(ginv(X0) %*% X)
zapsmall(t(ginv(X) %*% X0))


X %*% t(L_save)
L_save %*% X

X <- L
zero_col <- apply(X, 2, function(x) all(abs(x) < 1e-8))
x <- X[, !zero_col, drop=FALSE]
t(gram_schmidt(t(x), norm=FALSE))

x[, 1] <- x[, 1] / x[1, 1]

orthonorm(t(L))
nicify_cols(L)

t(nicify_cols(t(L)))
X <- L




proj(x[, 1L], x[, 2L])
x2 <- x[, 2L] - proj(x[, 1L], x[, 2L])


#####################################################################
# step method:

data("TVbo", package="lmerTest")
fm <- lmer(Coloursaturation ~ TVset + Picture + (1|Assessor:TVset) + (1|Assessor),
           data=TVbo)
(step_fm <- step(fm))
names(step_fm)
str(step_fm, max.level = 1)
fm2 <- get_model(step_fm)
ranova(fm2)
drop1(fm2)

summary(fm)
ranova(fm)

drop1

data("ham", package="lmerTest")
fm <- lmer(Informed.liking ~ Product*Information*Gender+
             (1|Consumer) + (1|Product:Consumer) + (1|Product:Information), data=ham)
(step_fm <- step(fm))
names(step_fm)
str(step_fm, max.level = 1)
fm2 <- get_model(step_fm)
ranova(fm2)
drop1(fm2)


red <- reduce_random(fm)
redf <- reduce_fixed(attr(red, "model"))
ran_redTable(red)
fix_redTable(redf)
model <- fm
str(red_fixed)
names(attributes(red_fixed))


step(fm)
## Add headers to ran and fix anova tables.


slotNames(fm)
1/eigen(fm@vcov_varpar)$values
eigen(fm@optinfo$derivs$Hessian)$values
rand(fm)
drop1(fm)



#####################################################################

library(MASS)

A <- t(rbind(c(1, 1, 1, 1),
             c(1, 2, 3, 4),
             c(-3, -1, 1, 3)))
A <- t(rbind(c(1, 1, 1, 1),
             c(1, 2, 3, 4),
             c(2, 2, 1, 3),
             c(1, 4, 9, 16)))
A
qr(A)$rank

(A0 <- A[, 1:2])

zapsmall(A %*% ginv(crossprod(A)) %*% t(A))
(PA <- zapsmall(A %*% ginv(A)))
(PA0 <- zapsmall(A0 %*% ginv(A0)))

PA - PA0

svd_A <- svd(A)
# Basis for the column space of A:
(C_A <- svd_A$u[, 1:4])
# Basis for the column space of A0:
svd_A0 <- svd(A0, nu = nrow(A0))
(C_A0 <- svd_A0$u[, 1:2])

# Are column vectors in C_A0 in C_A? Can we express C_A0[, 1] as a linear combination
# of the columns in C_A?
C_A0[, 1]
yhat <- C_A %*% ginv(C_A) %*% C_A0[, 1]

resid(lm.fit(x = C_A, y=C_A0))

resid(lm.fit(x = C_A0, y=C_A))

A0A <- cbind(A0, A)
u <- svd(A0A)$u

tcrossprod(u[, 3:4]) %*% A

(L <- t(u[, 3:4]) %*% A)

library(pbkrtest)
?KRmodcomp

methods(model2restrictionMatrix)

model2restrictionMatrix.lm
pbkrtest:::.restrictionMatrixBA

#####################################################################

cake2 <- cake
cake2$temperature <- factor(cake2$temperature, ordered = FALSE)
fm <- lmer(angle ~ recipe + temperature + (1|recipe:replicate), cake2)
show_tests(anova(fm, type=1))
fm0 <- lmer(angle ~ temperature + (1|recipe:replicate), cake2)
X <- model.matrix(fm)
X0 <- model.matrix(fm0)
(L0 <- diag(ncol(X))[2:3, ])


head(X)
head(X0)
head(zapsmall(R <- resid(.lm.fit(x=X0, y=X))), 20)
R <- R[, colSums(abs(R)) > 1e-8]
Lt <- crossprod(X, R) # t(X) %*% resid(.lm.fit(x=X0, y=X))
head(zapsmall(X0e <- X %*% Lt), 20) # The 'diff' of 'X-X0'


head(R)
zapsmall(t(X) %*% R)

(L7 <- t(zapsmall(qr.Q(qr(t(X) %*% R)))[, 1:2]))
contestMD(fm, L7)

zapsmall(get_Ldiffmat(X0, X))
(L1 <- zapsmall(get_Ldiffmat(X0, X)))
contestMD(fm, L0)
contestMD(fm, L1)
#####

## Null-space:
M <- cbind(X0, X)
tmp <- qr(M)
set <- if (tmp$rank == 0L) seq_len(ncol(M)) else -seq_len(tmp$rank)
qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]

svd_M <- svd(M)
str(svd_M)
svd_M$d
svd_M$u[]

head(resid(lm.fit(x=X0, y=X)))


getAnywhere(zapsmall)
base::zapsmall

(L <- zapsmall(t(resid(lm.fit(x=X0, y=X))) %*% X))
contestMD(L, fm)
L2 <- L[rowSums(L) > 1e-8, , drop=FALSE]
L2 <- zapsmall(t(qr.Q(qr(t(L2))))) # Orthonormalize
L2
contestMD(L2, fm)

dim(X0)
dim(X)
dim(Xhat <- X0 %*% ginv(X0) %*% X)
R <- X - Xhat
head(R)
dim(R)
all(t(X0) %*% R < 1e-8) # = 0 so X0 is orthogonal to R

qr(zapsmall(R))[-1]

colSums(abs(R))
(L <- t(R[, 2:3]) %*% X)
t(qr.Q(qr(t(L))))
qr(R)[-1]

NN <- nullspace(cbind(X0, X))
dim(NN)
zapsmall(NN)
KK <- zapsmall(nullspace(cbind(X0, X), type="left"))
dim(KK)

dim(MASS::Null(cbind(X0, X)))
Q <- qr.Q(qr(cbind(X0, X)))
Q2 <- Q[, 7:8]

qr(cbind(X0, X))$rank
head(qr.Q(qr(cbind(X0, X)), complete = TRUE)[, 1:14])
head(MASS::Null(cbind(X0, X))[, 1:14])
head(Q)

head(QQ2 <- MASS::Null(cbind(X0, X))[, 7:8])
head(Q2)



t(QQ2) %*% X

L[, 3] * L[, 4]


u <- svd(cbind(X0, X))$u
rX0 <- qr(X0)$rank
rX <- qr(X)$rank
u2 <- u[, (rX0 + 1):rX]
(L3 <- t(u2) %*% X)
t(X) %*% u2
contestMD(L3, fm)

L4 <- t(qr.Q(qr(t(L3))))
contestMD(L4, fm)

L %*% fixef(fm)


#####################################################################

data("sleepstudy", package="lme4")
fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
           sleepstudy)
# F-test of third coeffcients - I(Days^2):
contest(c(0, 0, 1), fm)
# Equivalent t-test:
contest(L=c(0, 0, 1), fm, joint=FALSE)
# Test of 'Days + I(Days^2)':
contest(L=diag(3)[2:3, ], fm)
# Other options:
contest(L=diag(3)[2:3, ], fm, joint=FALSE)
contest(L=diag(3)[2:3, ], fm, joint=FALSE, collect=FALSE)

# Illustrate a list argument:
L <- list("First"=diag(3)[3, ], "Second"=diag(3)[-1, ])
contest(L, fm)
contest(L, fm, collect = FALSE)
contest(L, fm, joint=FALSE, confint = FALSE)
contest(L, fm, joint=FALSE, collect = FALSE, level=0.99)

# Illustrate testing of estimability:
# Consider the 'cake' dataset with a missing cell:
data("cake", package="lme4")
cake$temperature <- factor(cake$temperature, ordered=FALSE)
cake <- droplevels(subset(cake, temperature %in% levels(cake$temperature)[1:2] &
                            !(recipe == "C" & temperature == "185")))
with(cake, table(recipe, temperature))
fm <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
fixef(fm)
# The coefficient for recipeC:temperature185 is dropped:
attr(model.matrix(fm), "col.dropped")
# so any contrast involving this coefficient is not estimable:
Lmat <- diag(6)
contest(Lmat, fm, joint=FALSE, check_estimability = TRUE)

(lsm7 <- lsmeans(model))
show_contrasts(lsm7)
show_contrasts(lsm7, names=FALSE)
summary(model)

data("soup", package="ordinal")
soup$sure <- as.numeric(soup$SURENESS)
fm <- lmer(sure ~ PRODID * DAY + (1 | RESP), data=soup)

length(fixef(fm))
Lmat <- structure(diag(length(fixef(fm))),
                  dimnames=list(names(fixef(fm)), names(fixef(fm))))
contest(Lmat[2, ], fm)
contest(Lmat, fm, joint=FALSE, confint = FALSE)
coef(summary(fm))

Lmat <- diag(length(fixef(fm)) + 1L)
contest(Lmat, fm, joint=FALSE, check_est=TRUE)

load_all(r2path)
model <- fm
which <- NULL
ddf <- "Sat"
level <- 0.95
rhs=0
confint <- TRUE
check_estimability <- TRUE

L <- L_list <- lsmeans_contrasts(fm)
cc <- contest(L_list, fm, joint=FALSE, check_estimability = TRUE, collect = FALSE)
rbindall(unname(cc))

#####################################################################
# Marten data example:

library(lmerTest)

d <- readRDS("../misc/lmm_data.rds")
str(d)
with(d, table(A, B))
set.seed(1234)
d2 <- d[sample(nrow(d), 1000), ]
with(d2, table(A, B))

m1 <- lmer(value ~ A * B + (1 + A|subject) + (1|item) + (1|subject:item), data = d2)
anova(m1, type="3") # works fine
anova(m1, type="3", ddf="KR") # works fine

m1 <- lmer(value ~ A*B + (1 + A|subject) + (1|item) + (1|subject:item), data = d2)
m1
anova(m1) # works fine
anova(m1, type="3c") # works fine
anova(m1, type="3c", ddf="KR") # works fine

#####################################################################
## Marginal contrasts and new type 2b
load_all(r2path)

data("soup", package="ordinal")
soup$sure <- as.numeric(soup$SURENESS)
fm <- lmer(sure ~ PRODID * DAY + (1 | RESP), data=soup)
model <- fm
containment(fm)
need_yates(fm)

anova(fm, type="1")
anova(fm, type="2")
anova(fm, type="3")
anova(fm, type="marginal")
anova(fm, type="yates") # = 3b

anova(fm, type="2b")
anova(fm, type="3b")
anova(fm, type="3c")

show_tests(anova(fm, type="1"))
show_tests(anova(fm, type="2"))
show_tests(anova(fm, type="2b"))
show_tests(anova(fm, type="3"))
show_tests(anova(fm, type="3b"))
show_tests(anova(fm, type="marginal"))

###################
data("cake", package="lme4")
fm <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
containment(fm)
need_yates(fm)

anova(fm, type="1")
anova(fm, type="2")
anova(fm, type="2b")
anova(fm, type="3")
anova(fm, type="3b")
anova(fm, type="marginal")

show_tests(anova(fm, type="1"))
show_tests(anova(fm, type="2"))
show_tests(anova(fm, type="2b"))
show_tests(anova(fm, type="3"))
show_tests(anova(fm, type="3b"))
show_tests(anova(fm, type="marginal"))

########################
cake2 <- cake
cake2$temperature <- factor(cake2$temperature, ordered = FALSE)
fm <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake2)
containment(fm)
need_yates(fm)

anova(fm, type="1")
anova(fm, type="2")
anova(fm, type="2b")
anova(fm, type="3")
anova(fm, type="3b")
anova(fm, type="3c")
anova(fm, type="marginal")

show_tests(anova(fm, type="1"), names = FALSE, fractions = TRUE)
show_tests(anova(fm, type="2"), names = FALSE, fractions = TRUE)
show_tests(anova(fm, type="2b"), names = FALSE, fractions = TRUE)
show_tests(anova(fm, type="3"), names = FALSE, fractions = TRUE)
show_tests(anova(fm, type="3b"), names = FALSE, fractions = TRUE)
show_tests(anova(fm, type="3c"), names = FALSE, fractions = TRUE)
show_tests(anova(fm, type="marginal"), names = FALSE)

## How to get type III contrasts in general (suggestion):
## - Make a design matrix using contr.SAS
## - Get the type III contrasts for the contr.SAS coding
## - map the coefficients from original coding to contr.SAS
##    - and evaluate the test with the mapped coefs.
## - alt. map the contrasts to original coding
##    - and evaluate the test with the original coefs.
##
## NOTE: any type II or III contrast should be independent of the order
##   in which the terms are listed in the model formula!
## This has to hold true even when there are aliased coefficients!
## Is type 3b really order independent?


fm <- lmer(sure ~ PRODID * DAY + (1 | RESP), data=soup)

anova(fm, type="3b")
fm <- lmer(sure ~ DAY * PRODID + (1 | RESP), data=soup)
anova(fm, type="3b") # Order independet so far

#######################

#####################################################################
## Soup

data("soup", package="ordinal")
str(soup)
soup$sure <- as.numeric(soup$SURENESS)

fm <- lmer(sure ~ PRODID * DAY + (1 | RESP), data=soup)


fm0 <- lm(sure ~ PRODID * DAY, data=soup)
colnames(model.matrix(fm0))
names(coef(fm0))
names(fm0)
methods(coef)
fm0$assign
stats:::coef.default
methods(coefficients)
summary(fm0)
stats:::print.lm
fm0$coefficients
names(summary(fm0))
summary(fm0)$aliased
stats:::summary.lm


fm <- lmer(sure ~ PRODID * DAY + (1 | RESP), data=soup)
fm <- lmer(sure ~ PRODID * relevel(DAY, "2") + (1 | RESP), data=soup)
summary(fm)
fixef(fm)

lapply(get_contrasts_type3b(fm), ncol)
colnames(model.matrix(fm))
names(fixef(fm))
anova(fm)
anova(fm, type="1")
anova(fm, type="2")
anova(fm, type="3")
anova(fm, type="3b")
attr(anova(fm, type="1"), "hypotheses")
attr(anova(fm, type="2"), "hypotheses")
attr(anova(fm, type="3"), "hypotheses")
attr(anova(fm, type="3b"), "hypotheses")

# Fit basic model to the 'cake' data:
data("cake", package="lme4")
fm1 <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
model <- fm1
# Type 3 anova table:
(an <- anova(fm1, type="3"))
(an <- anova(fm1, type="3c"))
# Display tests/hypotheses for type 1, 2, and 3 ANOVA tables:
# (and illustrate effects of 'fractions' and 'names' arguments)
show_tests(anova(fm1, type="1"))
show_tests(anova(fm1, type="2"))
show_tests(anova(fm1, type="3"))
show_tests(anova(fm1, type="3c"))
show_tests(anova(fm1, type="3b"))

show_tests(anova(fm1, type="3"), fractions = TRUE)
show_tests(anova(fm1, type="3"), names=FALSE)


show_tests(anova(fm, type=1))
show_tests(anova(fm, type=1), fractions = TRUE)
show_tests(anova(fm, type=2))
show_tests(anova(fm, type=3), names=FALSE)
show_tests(anova(fm, type=3), fractions = FALSE)
show_tests(anova(fm, type=3), fractions = TRUE)


tmp <- show_tests(anova(fm, type=2), names=TRUE)[["DAY"]]
MASS::fractions(tmp)
MASS::fractions(pi)

model <- fm
term <- which <- "DAY"

library(lsmeans)
lsmeans(fm, ~DAY:PRODID)
search()


XX <- model.matrix(terms(model), model.frame(model))
head(XX)
nonest.basis(XX)
getAnywhere(nonest.basis.qr)
attr(anova(fm, type="3"), "hypotheses")$DAY

is.estble(c(attr(anova(fm, type="3b"), "hypotheses")$DAY), nonest.basis(XX))

(nullspaceX <- zapsmall(nullspace(XX)))

library(estimability)
?estimability

warp.lm1 <- lm(breaks ~ wool * tension, data = warpbreaks,
               subset = -(26:38),
               contrasts = list(wool = "contr.treatment", tension = "contr.treatment"))
zapsmall(nonest.basis(warp.lm1))

warp.lm2 <- update(warp.lm1,
                   contrasts = list(wool = "contr.sum", tension = "contr.helmert"))
zapsmall(nonest.basis(warp.lm2))

#####################################################################
## ls-means (and optionally type 3 tests):

data("cake", package="lme4")
model <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
containment(model)
numeric_terms(model)
need_yates(model)
no_yates(model)


anova(model)
anova(model, type="1")
anova(model, type="2")
anova(model, type="3")
anova(model, type="3b")
object <- model
ddf <- "Satterthwaite"

Llist <- get_contrasts_type3b(model)
rbindall(lapply(Llist, contestMD, model=model))
anova(model, type=3)


numeric_terms(model)

(type2_terms <- names(is_contained)[!is_contained])
(yates_terms <- names(is_contained)[is_contained])

X <- model.matrix(model)
data_classes <- attr(terms(m, fixed.only=FALSE), "dataClasses")
Terms <- terms(model)

Llist_2 <- get_contrasts_type2(X, Terms, data_classes, which=type2_terms)
rbindall(lapply(Llist_2, contestMD, model=model))

anova(model, type=3)
attr(anova(model, type=2), "hypotheses")

which <- need_yates(model)
L3list <- get_yates_contrast(model, which=need_yates(model))
rbindall(lapply(L3list, contestMD, model=model))


## FIXME: Add "Response: ..."-line to anova heading

var_list <- get_var_list(model)
grid <- get_min_data(model)
uX <- model.matrix(nobars(formula(model)[-2]), data=grid)
term <- yates_terms[1]
Llist_3 <- lapply(setNames(as.list(yates_terms), yates_terms), function(term) {
  Lt <- model.matrix(formula(paste0("~ 0 + ", term)), data=grid)
  t(sweep(Lt, 2, 1/colSums(Lt), "*")) %*% uX
})
rbindall(lapply(Llist_3, contestMD, model=model))


mf <- model.frame(model)
X <- model.matrix(model)
head(mf)
(data_classes <- attr(terms(m, fixed.only=FALSE), "dataClasses"))
Terms <- terms(model)
Terms
(term_names <- attr(Terms, "term.labels"))
# make list of all terms with levels of the factors and mean values of numeric
#   variables
var_list <- list(recipe=levels(mf[["recipe"]]),
                 temp = mean(mf[["temp"]]))
var_list
grid <- do.call(expand.grid, var_list)
grid
uX <- model.matrix(nobars(formula(model)[-2]), data=grid)
term <- term_names[3]
Llist <- lapply(setNames(as.list(term_names), term_names), function(term) {
  Lt <- model.matrix(formula(paste0("~ 0 + ", term)), data=grid)
  t(sweep(Lt, 2, 1/colSums(Lt), "*"))
})
Llist
lapply(Llist, function(L) L %*% uX)

Lt <- model.matrix(formula("~ 0 + recipe"), data=grid) # contrast vectors for each level of sex in cols
Lt <- sweep(Lt, 2, 1/colSums(Lt), "*") # Take flat average by cols to get 'LS-means'
(L <- t(Lt) %*% uX)
L2 <- rbind(L[2, ] - L[1, ],
            L[3, ] - L[1, ])
L %*% fixef(model)
rbindall(lapply(1:nrow(L), function(i) contest1D(L[i, ], model)))
L2 <- cbind(c(1, 1),
            c(1, -1),
            c(0, -1)) %*% L
contestMD(L2, model)
head(X)
## How do I get the type III test from this LS-means contrast matrix?

Lt <- model.matrix(formula("~ 0 + recipe:temp"), data=grid) # contrast vectors for each level of sex in cols
Lt <- sweep(Lt, 2, 1/colSums(Lt), "*") # Take flat average by cols to get 'LS-means'
L <- t(Lt) %*% uX
contestMD(L, model)

Lt <- model.matrix(formula("~ 0 + temp"), data=grid) # contrast vectors for each level of sex in cols
Lt <- sweep(Lt, 2, 1/colSums(Lt), "*") # Take flat average by cols to get 'LS-means'
(L <- t(Lt) %*% uX)
L %*% fixef(model)
contest1D(L, model)
contestMD(L, model)


(an <- anova(model, type=3))
attr(an, "hypotheses")
emmeans::emmeans(model, "recipe")
emmeans::emmeans(model, "temp")
emmeans::emmeans(model, ~temp:recipe)
library(lsmeans)
lsmeans(model, ~temp:recipe)

lsmeans(model, pairwise ~ temp:recipe)
test(pairs(lsmeans(model, ~ temp:recipe)), joint=TRUE)
test(pairs(lsmeans(model, ~ recipe)), joint=TRUE)
test(pairs(lsmeans(model, ~ temp)), joint=TRUE)

methods(contrast)
lsmeans:::pairs.lsm.list
lsmeans:::pairs.ref.grid
lsmeans:::contrast.ref.grid
lsmeans:::pairwise.lsmc
pairwise.lsmc(letters[1:3])


###################################
## Case with two factors:
cake2 <- cake
cake2$temperature <- factor(cake2$temperature, ordered=FALSE)
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake2)
containment(model)
numeric_terms(model)
need_yates(model)
no_yates(model)

Llist <- get_contrasts_type3b(model)
rbindall(lapply(Llist, contestMD, model=model))
anova(model, type=3)

mf <- model.frame(model)
X <- model.matrix(model)
(data_classes <- attr(terms(model, fixed.only=FALSE), "dataClasses"))
Terms <- terms(model)
(term_names <- attr(Terms, "term.labels"))
# make list of all terms with levels of the factors and mean values of numeric
#   variables
var_list <- list(recipe=levels(mf[["recipe"]]),
                 temperature = levels(mf[["temperature"]]))
var_list

get_var_list(model)
get_min_data(model)

grid <- do.call(expand.grid, var_list)
grid
uX <- model.matrix(nobars(formula(model)[-2]), data=grid)
term <- term_names[3]
Llist <- lapply(setNames(as.list(term_names), term_names), function(term) {
  Lt <- model.matrix(formula(paste0("~ 0 + ", term)), data=grid)
  t(sweep(Lt, 2, 1/colSums(Lt), "*"))
})
Llist
Llist2 <- lapply(Llist, function(L) L %*% uX)
lapply(Llist2, unname)
str(Llist2)

x <- Llist2[[3]]
LLL <- t(get_trts(rownames(x))) %*% x
qr(LLL)$rank
lapply(Llist2, function(L) {
  unname(t(get_trts(rownames(L))) %*% L)
})
contestMD(x, model)
contestMD(LLL, model)
anova(model, type=3)
anova(model, type=1)
unname(attr(anova(model, type=3), "hypotheses")[[3]])

lapply(attr(anova(model, type=2), "hypotheses"), unname)
lapply(attr(anova(model, type=3), "hypotheses"), unname)

lapply(Llist2, function(LL)
  rbindall(lapply(1:nrow(LL), function(i) contest1D(LL[i, ], model))) )

## Could a Type III algorithm be:
## - For all terms that ARE NOTE contained in other terms, se Type II
## - For terms that ARE contained in other terms use a full rank contrast of
##    the LS-means expression.

## Can continuous variables be contained in other terms? Yes.

test(pairs(lsmeans(model, ~ recipe)), joint=TRUE)
test(pairs(lsmeans(model, ~ temperature)), joint=TRUE)
test(pairs(lsmeans(model, ~ temperature:recipe)), joint=TRUE)
test(lsmeans(model, ~ temperature:recipe), joint=TRUE)


emmeans::emmeans(model, "recipe")
emmeans::emmeans(model, "temperature")
emmeans::emmeans(model, ~recipe*temperature)

###################################
## Case with ordered factor:
model <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
terms(model)
attr(terms(model, fixed.only=FALSE), "dataClasses")
head(model.matrix(model))
containment(model)
numeric_terms(model)
need_yates(model)
no_yates(model)
.getXlevels(terms(model), model.frame(model))

Llist <- get_contrasts_type3b(model)
lapply(Llist, unname)
rbindall(lapply(Llist, contestMD, model=model))
LTmodel <- lmerTest::lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
lmerTest::anova(LTmodel, type=3)
## FIXME make Type 3b work for other types of contrasts as well!

## FIXME: ensure Type 3b (numeric_terms) work for models with nested factors
##       It seems it should: terms(y ~ a/b)

lsmeans(model, ~recipe)
lsmeans(model, ~temperature)

###################################
## LS-means via

#####################################################################
## Type II anova:

## FIXME: Improve attr(anova(m, type=Z), "hypotheses") such that all terms
## - are matrices
## - have column names corresponding to the parameters.

## FIXME: Test what happens with aliased coefs (cf. soup example from ordinal)

##########
## Completely balanced dataset in which:
## - type II and III agree
## - Type I is (partially) different from Type II and III
## - car::Anova type II agree with lmerTest::anova Type I
## - Satterthwaite and Kenward-Roger are identical
## - lme4 computes the same F-values for the Type I table
## - changing the order of the terms in the model does not affect any tests
data("cake", package="lme4")
with(cake, table(recipe, temp)) # Balance
with(cake, table(recipe, replicate)) # Balance
with(cake, table(temp, replicate)) # Balance
m <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
anova(m, type=1)
anova(m, type=3)
anova(m, type=2)
car::Anova(m, test="F", type=2)
# Using Kenward-Roger:
anova(m, type=1, ddf="KR")
anova(m, type=2, ddf="KR")
anova(m, type=3, ddf="KR")
# lme4 Type I table:
anova(m, ddf="lme4")
# Changing the order of terms in the model:
m2 <- lmer(angle ~  temp * recipe + (1|recipe:replicate), cake)
anova(m2, type=1)
anova(m2, type=2)
anova(m2, type=3)
car::Anova(m2, test="F", type=2)

## Type I differ from II & III since
## - Type I computes the contrast for recipe as if it is contained in recipe:temp
## - Type II assumes temp is contained in recipe:temp since they are both
##    continuous but recipe:temp have more factors than temp
## - Type II assumes recipe is NOT contained in recipe:temp since these terms
##    do not have the same continuous variables AND recipe:temp does not have
##    more factors than recipe alone.
containment(m)
attr(anova(m, type=1), "hypotheses")
attr(anova(m, type=2), "hypotheses")

m3 <- lm(angle ~  recipe * temp, cake)
containment(m3)


##########


#########
# Balanced example in which
data("cake", package="lme4")
with(cake, table(recipe, temp))
with(cake, table(recipe, replicate))
with(cake, table(temp, replicate))
m <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
anova(m, type=1)

anova(m, type=3)
lmerTest::anova(m, type=3)

anova(m, type=2)
lmerTest::anova(m, type=2)
car::Anova(m, test="F", type=2)


m <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
m <- lmer(angle ~ temp * recipe + (1|recipe:replicate), cake)
X <- model.matrix(m)
terms <- terms(m)
data_classes <- attr(terms(m, fixed.only=FALSE), "dataClasses")
get_contrasts_type1(X, terms)
lapply(get_contrasts_type2(X, terms, data_classes), function(x) {
  colnames(x) <- NULL
  x
})

m1 <- lmer(angle ~ recipe * temp + (1|recipe:replicate), cake)
m2 <- lmer(angle ~ temp * recipe + (1|recipe:replicate), cake)
attr(anova(m1, type=2), "hypotheses")
attr(anova(m2, type=2), "hypotheses")
attr(anova(m1, type=1), "hypotheses")
attr(anova(m2, type=1), "hypotheses")
attr(anova(m1, type=3), "hypotheses")
attr(anova(m2, type=3), "hypotheses")

anova(m1, type=2)
anova(m2, type=2)
anova(m1, type=1)
anova(m2, type=1)
car::Anova(m1, test="F", type=2)
car::Anova(m2, test="F", type=2)
anova(m1, type=3)
anova(m2, type=3)

t(doolittle(crossprod(X))$L)

lme4::terms.merMod()

#########

data("TVbo", package="lmerTest")
fm <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor:TVset) + (1|Assessor),
           data=TVbo)
anova(fm, type=2)
anova(fm, type=3)
lmerTest:::anova(fm, type=2)

data("ham", package="lmerTest")
m <- lmer(Informed.liking ~ Gender * Information + Product +(1|Consumer), data=ham)
anova(m, type=2)
anova(m, type=3)
car::Anova(m, type=2, test="F")
anova(m, type=2, ddf="KR")
?car::Anova
lmerTest:::anova(m, type=2)

data("carrots", package="lmerTest")
carrots$gender <- factor(carrots$Gender)
carrots$income <- factor(carrots$Income)
str(carrots)

fm <- lmer(Preference ~ sens2 + Homesize*gender + (1+sens2 | Consumer),
           data=carrots)
anova(fm, type=2)
anova(fm, type=3)
car::Anova(fm, type=2, test="F")
anova(fm, type=2, ddf="KR")
car::Anova(fm, type=3, test="F") # NOT actually type III tests!
anova(fm, type=3, ddf="KR")

# Here is a difference between car::Anova(..., type=2) and anova2
fm <- lm(Preference ~ sens2 + Homesize*gender + Income*gender, data=carrots)
anova(fm)
car::Anova(fm, type=2)
anova2(fm)


fm <- lm(Informed.liking ~ Gender * Information * Product, data=ham)
anova2(fm)
anova(fm)
car::Anova(fm, type=2)
object <- fm

## A + B + C + A:B
##
## A:B - The conventional test of A:B (last) - R(A:B | A, B, C)
## C - Adjusted for A, B and A:B - R(C | A, B, AB)
## A - Adjusted for B and C, but not A:B - R(A | B, C)
## B - Adjusted for A and C, but not A:B - R(B | A, C)

## If a term ('T'; e.g. C or A:B) IS NOT CONTAINED in any other term, move the
## columns in X corresponding to T last, and obtain the Type I contrast for T
## using doolittle.
##
## If a term ('T'; e.g. A or B) IS CONTAINED in another term (such as A:B),
## remove the columns in X corresponding to the containing terms, move the
## columns in X corresponding to T last, and obtain the Type I contrast for T
## using doolittle.
##
## A Type II test is a test for a term T that adjusts for all terms that do not
## contain T.

## Containment:
## F1 (A) is contained in F2 (A:B) [F2 (A:B) contains F1 (A)] if
## - F1 and F2 involve the same X's (if any)
## - F2 involve more factors (say B) than F1
## - All factors in F1 (if any) are part of F2

## A + X + A:X
## - A:X is not contained in any term
## - A is not contained in A:X since they do not involve the same X's
## - X is contained in A:X
##
## R(A:X | A, X)
## R(A | X, A:X)
## R(X | A)

## A + B + C + X + A:B + A:X
##
## R(A | B, C, X)
## R(B | A, C)
## R(C | A, B, AB)
## R(X | A, B, C, A:B) - not A:X
## R(A:B | A, B, C, X, A:X)
## R(A:X | A, B, C, X, A:X)
##


## For each term
## -

#####################################################################
?ranova

## Add BIC to model list?
## Only add AIC and BIC when REML=FALSE?
## AIC: smaller is better!
## - lmer.anova seems to include REML-AIC when models are fitted with REML,
##   so I think we could do that as well. For lm lmer.anova only reports AIC
##   for ML fits.

################################

## Model fails to converge with 1 negative eigenvalue:
m <- lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
            (1|Consumer) + (1 |Consumer:Income), data=carrots)
anova(m, type=3)
summary(m)
eigen(m@vcov_varpar)$values # exactly singular
ranova(m)

## Nearly unidentifiable fit:
m <- lmer(Preference ~ sens2 + Homesize + Income + (1 + sens2 | Consumer) +
            (1 + Income |Consumer), data=carrots)
eigen(m@optinfo$derivs$Hessian)$values
ranova(m)
devfun <- update(m, devFunOnly=TRUE)
H <- numDeriv::hessian(func=devfun, x=m@theta)
eigen(H)$values

#####################################################################
# Looking at case where model fails to converge with a (slightly) negative
# eigenvalue

m <- lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
            (1|Consumer) + (1 |Consumer:Income), data=carrots)

# Compute Hessian:
m <- lme4::lmer(Preference ~ sens2 + Homesize + Income + #(1 + sens2 | Consumer) +
                  (1|Consumer) + (1 |Consumer:Income), data=carrots)
devfun <- update(m, devFunOnly=TRUE)
## Hessian for theta:
H <- numDeriv::hessian(func=devfun, x=m@theta)
eigen(H)
eigen(m@optinfo$derivs$Hessian)
m@optinfo$derivs$gradient
## Hessian for varpar:
is_reml <- getME(m, "is_REML")
varpar_opt <- unname(c(m@theta, sigma(m)))
h <- numDeriv::hessian(func=devfun_vp, x=varpar_opt, devfun=devfun,
                       reml=is_reml)
eigen(h)

################################
## lme4::lmer fails to detect lack of convergence:
m <- lmer(Preference ~ sens2 + Homesize +
            (1 + sens1 + sens2 | Consumer) + (1 | Homesize:Consumer),
          data=carrots) ## Error message.
m <- lme4::lmer(Preference ~ sens2 + Homesize +
                  (1 + sens1 + sens2 | Consumer) + (1 | Homesize:Consumer),
                data=carrots) ## No error message

eigen(m@optinfo$de$Hessian)$values
devfun <- update(m, devFunOnly=TRUE)
H <- numDeriv::hessian(func=devfun, x=m@theta)
eigen(H)$values
eigen(m@optinfo$derivs$Hessian)$values

#####################################################################
## lmerTest update:

# - New implementation of 90% of the code
# - Less error prone code
# - better test coverage
# - higher accuracy
# - better model/convergence and error checking
# - Bug fixes
#   o
# - New implementation of LS-means (available in lsmeansLT)
# - lsmeans has been defunct
# - ?New algorithm for computing type III ANOVA tables
# - ?New handling of terms such as poly() and ordered factors
# - New facilities:
#   o ranova - random-effect ANOVA-like table (formerly rand(); now deprecated)
#   o ?drop1 using Satterthwaite F-tests for fixed terms
#   o ?ditto add1
#   o contest() aka test_contrast() for tests of linear functions of the fixed
#     effects parameters using Satterthwaites method for degrees of freedom
#   o anova for model comparison use F-tests for when models differ only
#     in the fixed effects.
#   o
#####################################################################


####################
# Test KR method for contestMD:
fm <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
           sleepstudy)
anova(fm)
anova(fm, ddf="KR")
# Define 2-df contrast - since L has 2 (linearly independent) rows
# the F-test is on 2 (numerator) df:
L <- rbind(c(0, 1, 0), # Note: ncol(L) == length(fixef(fm))
           c(0, 0, 1))
# Make the 2-df F-test of any effect of Days:
contestMD(L, fm)
# Make the 1-df F-test of the effect of Days^2:
contestMD(L[2, , drop=FALSE], fm)

contestMD(L, fm, method="KR")
contestMD(L, fm, method="Sat")

contestMD(L[2, , drop=FALSE], fm, method="KR")

contestMD(L[0, , drop=FALSE], fm, method="Sat")
contestMD(L[0, , drop=FALSE], fm, method="KR")

contestMD(rbind(c(1, 0, 0)), fm)
contestMD(rbind(c(1, 0, 0)), fm, method="KR")
####################


anova(m, ddf="KR")
anova(m, ddf="lme4")
assertError(anova(m, ddf="Ken"))
anova(m, type=1)

## Looking at singular fits
nrow(sleepstudy)
g <- factor(rep(1:2, c(170, 10)))
m <- lmer(Reaction ~ Days +  (Days | Subject) + (1|g), sleepstudy)
m
anova(m)
vc <- as.data.frame(VarCorr(m))
stopifnot(isTRUE( # check that fit has a zero variance
  all.equal(0, vc[vc$grp == "g", "sdcor"])
))
# The hessian/vcov is actually positive definite:
stopifnot(isTRUE(
  all(eigen(m@A, only.values = TRUE)$values > 0)
))

model <- lme4::lmer(Reaction ~ Days +  (Days | Subject) + (1|g), sleepstudy)
model


####################

#################
## DamonData.R

load(system.file("testdata", "ttDamon.RData", package="lmerTest"))

head(ttDamon[, c("condition", "v_matri_plus_n_indef_freq_sq", "Subject")])
summary(ttDamon[, c("condition", "v_matri_plus_n_indef_freq_sq", "Subject")])
str(ttDamon)
dim(ttDamon)

x <- ttDamon$v_matri_plus_n_indef_freq_sq
xm <- x - mean(x)
xs <- scale(x, center = TRUE, scale = TRUE)

fm1 <- lmer(RT2LogR ~ condition * v_matri_plus_n_indef_freq_sq+
              (1+condition|Subject), data=ttDamon, REML=TRUE)
fm2 <- lmer(RT2LogR ~ condition * xm +
              (1+condition|Subject), data=ttDamon, REML=TRUE)
fm3 <- lmer(RT2LogR ~ condition * xs +
              (1+condition|Subject), data=ttDamon, REML=TRUE)

(an.sat <- anova(fm1))
(an.sat2 <- anova(fm2))
(an.sat3 <- anova(fm3))

anova(fm1, ddf="lme4")
anova(fm2, ddf="lme4")
anova(fm3, ddf="lme4")

TOL <- 1e-2 # for the check

## with 4 decimals should agree with SAS output
## numbers before decimals should agree with SAS output
stopifnot(
  all.equal(an.sat[,"DenDF"], c(76.206, 6074.6, 6146.1), tol = TOL),
  all.equal(an.sat[, "F.value"], c(9.2182, 2.8360, 6.6544), tol = TOL)
  , TRUE)

#########################
## MAMSatterthKR.R

# require(lmerTest)
testType1 <- TRUE

data("TVbo", package="lmerTest")

lm.pred <- lm(as.formula(paste("Lightlevel", "~",
                               paste(c("TVset","Picture"), collapse="*"), sep="")),
              data=TVbo)
TVbo$x <- scale(predict(lm.pred), scale=FALSE)
lmerTVpic <- lmer(Lightlevel ~ TVset*Picture +   Assessor:x  + (1|Assessor) +
                    (1|TVset:Assessor) + (1|Picture:Assessor) +
                    (1|TVset:Picture:Assessor), data=TVbo)
anova(lmerTVpic, ddf="lme4")
anova(lmerTVpic, type=1)
anova(lmerTVpic, ddf="Kenward-Roger", type=1)
# Analysis of Variance Table of type I  with  Kenward-Roger
# approximation for degrees of freedom
#               Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)
# TVset         37.743 18.8715     2  8.400  35.631 7.886e-05 ***
# Picture       16.290  5.4301     3 18.007  10.252 0.0003663 ***
# TVset:Picture 15.454  2.5756     6 39.144   4.863 0.0008495 ***
# Assessor:x     2.279  0.3255     7 26.800   0.615 0.7389038
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## TODO: check with SAS dfs for Satterthwaite and KR to agree
if(testType1){
  tools::assertWarning(anova(lmerTVpic, type=1)) ## warning: ddf=0 for TVset
  anova(lmerTVpic, type=1, ddf="Kenward-Roger")
}

