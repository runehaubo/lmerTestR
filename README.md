# lmerTest - Tests in Linear Mixed Effects Models

This is the repo for the _new_ **lmerTest** package, the old package is available [here](https://github.com/runehaubo/lmerTest).

[![Build Status](https://travis-ci.org/runehaubo/lmerTestR.svg?branch=master)](https://travis-ci.org/runehaubo/lmerTestR)
[![cran version](http://www.r-pkg.org/badges/version/lmerTest)](https://cran.r-project.org/package=lmerTest)
[![downloads](https://cranlogs.r-pkg.org/badges/lmerTest)](https://cran.r-project.org/package=lmerTest)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/lmerTest)](http://cranlogs.r-pkg.org/badges/grand-total/lmerTest)
[![Research software impact](http://depsy.org/api/package/cran/lmerTest/badge.svg)](http://depsy.org/package/r/lmerTest)

## Main features

The **lmerTest** package provides _p_-values in type I, II or III `anova` and `summary`
tables for linear mixed models (`lmer` model fits cf. **lme4**) via Satterthwaite's degrees of freedom method; a Kenward-Roger method is also available via the **pbkrtest**
package. Model selection and assessment methods include `step`, `drop1`, anova-like 
tables for random effects (`ranova`), least-square means (LS-means; `ls_means`) 
and tests of linear contrasts of fixed effects (`contest`).

## Citation

To cite **lmerTest** in publications use:

Kuznetsova A., Brockhoff P.B. and Christensen R.H.B. (2017). "lmerTest Package: Tests in Linear Mixed Effects Models." _Journal of Statistical Software_, 82(13), pp. 1–26. doi: 10.18637/jss.v082.i13.

Corresponding BibTeX entry:

    @Article{,
      title = {{lmerTest} Package: Tests in Linear Mixed Effects Models},
      author = {Alexandra Kuznetsova and Per B. Brockhoff and Rune H. B.
        Christensen},
      journal = {Journal of Statistical Software},
      year = {2017},
      volume = {82},
      number = {13},
      pages = {1--26},
      doi = {10.18637/jss.v082.i13},
    }

## Discovered a bug?

Please raise a new issue! Preferably add code that illustrates the problem using one of the datasets from the **lmerTest**.

## Installation

Basically there are two options for installing **lmerTest**:

1. Released (stable version) from CRAN: in **R** run `install.packages("lmerTest")`.
2. Development version from GitHub: First load the **devtools** package (and install it if you do not have it) and install the default (master) branch:
```
library("devtools")
install_github("runehaubo/lmerTestR")
```
If you haven't already installed a previous version of **lmerTest** you need to also install dependencies (other packages that **lmerTest** depends on and requires you to install to function properly). We recommend that you install **lmerTest** from CRAN (using `install.packages("lmerTest")`) before installing from GitHub as described above. 

An alternative is to use 
```
library("devtools")
install_github("runehaubo/lmerTestR", dependencies=TRUE)
```
but that requires you to install all dependent packages from source (which only works if you have the correct compilers installed and set up correctly); installing the pre-compiled packages from CRAN is usually easier.


  

