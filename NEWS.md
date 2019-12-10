lmerTest 3.1-1
------------------

- Sofie P Jensen is taking over as maintainer replacing Per B Brockhoff.
- Fixing "noLD" CRAN issue (a check that ensures the package works on systems without long doubles). This was caused by an over sensitive test.

lmerTest 3.1-0
------------------

- Adding support for legacy model fits, i.e. `merModLmerTest` objects generated with lmerTest version `< 3.0-0`. This includes defining the `merModLmerTest` class and `anova`, `summary`, `drop1`, `ls_means`, `lsmeansLT` and `difflsmeans` methods. The usual `lme4` methods also work with objects of class `merModLmerTest`.

lmerTest 3.0-1
------------------

- over-sensitive tests (failing on Solaris) have reduced tolerance
- `sigma` and `sigma.merMod` defined and exported for `R <= 3.3.0`
- Warn if Kenward-Roger is used with `R <= 3.3.0` since it may give incorrect results
- Add `lme4 (>= 1.1-10)` and `R (>= 3.2.5)` to `Depends` (last available version where `lmerTest` checks out)
- `pbkrtest` package loaded conditional on availability in tests


lmerTest 3.0-0
------------------
 
  * The new and completely re-written lmerTest package. Details of changes are 
    available in [pdf](https://github.com/runehaubo/lmerTestR/blob/master/pkg_notes/new_lmerTest.pdf) or [html](http://htmlpreview.github.io/?https://github.com/runehaubo/lmerTestR/blob/master/pkg_notes/new_lmerTest.html)
  

lmerTest < 3.0-0
------------------

  * Signficant news and changes for the 2.0-xx release series is provided below.

2.0-34

- included citation info for JSS

2.0-33

- lsmeans and difflsmeans are now deprecated functions. Changed the names to lsmeansLT and dlsmeansLT
- changed the maintainer field


2.0-32

- changed the message of identifiability to the more appropriate one

2.0-31

- removed lmerTestFunctions.R and restructured the package. added calcSatterth(model, L) for  calculating Satterthwaite's approximation for a specified L matrix 

2.0-30

- envir.R failed with the newest version of lme4. Changed the code to pass the check. TODO: remove updating the model

2.0-28

- changes in general summary function.  callNextMethod changed to as(model, "lmerMod)

2.0-25 

- updated according to comments from CRAN

2.0-24 changes:

- cleaned the code


2.0-23 changes:

- hessian and grad changed to mygrad and myhess (deriv.R functions of Rune)
- plots use ggplot2
- look for  previous changes in R-Forge

2.0-11 changes:

- elimRandEffs deleted. now the rand table contains all the information

2.0-9 changes:

- fixed.calc option is added to step function
- elimRand effs changed: random effects that are 1 approx to 1e-6 are eliminated
- las=2 in barplots: verical axis names
contrast with the name "l" changed to "l.lmerTest.private.contrasts"

2.0-8 changes:

- throws error for lsmeans, difflsmeans, rand and step functions if the model does not inherit lmerMod class

2.0-7 changes:

- in utils calcSatterth changed: solve of 0 dim matrix now catches in tryCatch - example MAMex.R in tests is added to check the bug
- messages are printed if some computational errors occurr in anova or summary and the ones from lme4 are returned (bugSummary.R for testing)

2.0.6 changes:

- added a number of tests in the tests folder and inst/datasets for the testing data sets - will not be included in the R-forge nor CRAN (for a moment)
- model is not updated automatically to REML (tests for random effects are ML!)
- man functions updated

2.0.5 changes:

- fixed bug from Ben - summary(model, "lme4") changed to summary(model, ddf="lme4")
- fixed bug for summary from Cyrus
- added in manual notes regarding random coefficient models simplification
- Rune changed solve to chol2inv in lmerTestFunctions.R
- changed updateModel function so that the bugs with the environmentgs are solved

2.0.4 new:

- rewritten rand table elimination
- added elimrand.R


Modifications in lmerTest 2.0.1

- The elim.num column now has KEEP instead of 0
- X'X deficiancy was fixed by Rune, lmerTest was fixed accordingly 
