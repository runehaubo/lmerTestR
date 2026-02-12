---

# Question (2026-02-11)

I have an R package lmerTest (see https://github.com/runehaubo/lmerTestR) which Depends on lme4 (see https://github.com/lme4/lme4). An update to lme4 from version 1.1-38 (see https://cran.r-project.org/web/packages/lme4/index.html) to version 2.0.0 on Github caused errors in lmerTest. I have tried to accommodate the lme4-update in lmerTest by adding the following to the .onLoad function: 

"  ns <- parent.env(environment())
  if (packageVersion("lme4") < "2.0-0")
    assign("forceNewMerMod", envir = ns, inherits = FALSE,
           function(object, reference) object)" 

and the following to the NAMESPACE file: 

"if (utils::packageVersion("lme4") >= "2.0-0")
  importFrom(lme4,forceNewMerMod)" 

but it seems that `forceNewMerMod` is not imported from lme4 if lmerTest was build and installed when lme4 < 2.0-0 was installed and `forceNewMerMod` is not assign'ed when loading lmerTest if lme4>=2.0-0 is (later) installed. Therefore the following sequence fails: 

"remove.packages("lme4")
remove.packages("lmerTest")
install.packages("lme4")
install.packages("lmerTest")
devtools::install_github("lme4/lme4")
data("sleepstudy", package = "lme4")
lmerTest::lmer(Reaction ~ 1 + (1|Subject), sleepstudy)" 

with the error: 

"Error in forceNewMerMod(res, reference = model) : 
  could not find function "forceNewMerMod"
". 

How can I fix lmerTest such that the package works regardless of install sequence while I am not prepared to make lmerTest depend on the version of lme4 by adding something like `Depends: lme4 (>= 2.0-0)` in the `DESCRIPTION` file? 

---

# Answer from ChatGPT-5.2

You are running into a fundamental R packaging constraint:

* `NAMESPACE` is evaluated at *build time*.
* `.onLoad()` runs at *load time*.
* The set of imported symbols is fixed in the installed package.
* Changing `lme4` **after** installing `lmerTest` does not re-resolve imports.

So your conditional `importFrom()` in `NAMESPACE` cannot work reliably across install sequences. Once `lmerTest` is built without importing `forceNewMerMod`, it will never be available via import—even if a newer `lme4` is installed later.

This is expected behavior.

---

# Correct Pattern: Resolve at Runtime via `getFromNamespace()`

If you do **not** want to version-depend on `lme4 (>= 2.0-0)`, the correct approach is:

* Do **not** use conditional `importFrom()`
* Do **not** rely on static imports
* Resolve the symbol dynamically at runtime

In other words: treat `forceNewMerMod()` as an optional API.

---

## Recommended Fix

### 1️⃣ Remove this from NAMESPACE

Remove entirely:

```r
if (utils::packageVersion("lme4") >= "2.0-0")
  importFrom(lme4,forceNewMerMod)
```

NAMESPACE must not contain conditional logic like this.

---

### 2️⃣ Define a wrapper in your package

Inside `lmerTest` (e.g., in `zzz.R` or a small utility file), define:

```r
forceNewMerMod_safe <- function(object, reference) {
  if (utils::packageVersion("lme4") >= "2.0-0") {
    getFromNamespace("forceNewMerMod", "lme4")(object, reference)
  } else {
    object
  }
}
```

Then replace all calls to:

```r
forceNewMerMod(...)
```

with:

```r
forceNewMerMod_safe(...)
```

---

## Why This Works

* `getFromNamespace()` resolves the symbol dynamically at runtime.
* No static import required.
* Works regardless of installation order.
* No dependency version bump required.
* Avoids stale namespace issues.

---

## Even More Robust Version (Recommended)

Avoid relying purely on `packageVersion()`; instead test existence:

```r
forceNewMerMod_safe <- function(object, reference) {
  if (exists("forceNewMerMod", envir = asNamespace("lme4"), inherits = FALSE)) {
    getFromNamespace("forceNewMerMod", "lme4")(object, reference)
  } else {
    object
  }
}
```

This is better because:

* It guards against development versions.
* It guards against API removal.
* It avoids fragile version string comparison.

---

# Why Your `.onLoad()` Approach Fails

This code:

```r
ns <- parent.env(environment())
assign("forceNewMerMod", envir = ns, inherits = FALSE,
       function(object, reference) object)
```

fails because:

1. It assigns into the wrong environment.
2. It assumes the symbol will be visible in the calling scope.
3. Static imports (or lack thereof) dominate namespace resolution.
4. Reinstalling `lme4` does not rebind symbols inside `lmerTest`.

In short: namespace resolution is frozen at install time.

---

# Alternative Pattern (If You Want Minimal Refactoring)

You can bind the function at load time:

```r
.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  if (exists("forceNewMerMod", envir = asNamespace("lme4"), inherits = FALSE)) {
    assign("forceNewMerMod",
           getFromNamespace("forceNewMerMod", "lme4"),
           envir = ns)
  } else {
    assign("forceNewMerMod",
           function(object, reference) object,
           envir = ns)
  }
}
```

This works because you always define a symbol in your namespace.

However, the explicit wrapper function is cleaner and easier to maintain.

---

# What You Should NOT Do

* ❌ Conditional `importFrom()` in NAMESPACE
* ❌ Depend on install order
* ❌ Assume reinstalling `lme4` refreshes bindings
* ❌ Version-string logic alone (fragile)

---

# Best-Practice Summary

If you support multiple upstream API versions:

* Avoid static imports of version-dependent symbols
* Use `getFromNamespace()` dynamically
* Guard using `exists()` in namespace
* Keep `DESCRIPTION` dependency broad

This is the standard compatibility pattern used in packages that straddle API transitions (e.g., Matrix, ggplot2 internals, rlang evolution).

---

If you want, I can also suggest a small internal compatibility layer pattern that keeps all lme4-version branching isolated in a single file.

