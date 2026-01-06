# test_devfun_vp.R

# library(devtools)
# # has_devel()
# r2path <- "~/GitHub/lmerTestR/lmerTest"
# # document(pkg=r2path)
# load_all(r2path)
# 
# # ?`Covariance-class`

library(lmerTest)

has_pkgs <- requireNamespace("utils", quietly = TRUE) && 
  requireNamespace("numDeriv", quietly = TRUE) && 
  requireNamespace("lme4", quietly = TRUE)
is_lme4_2_0_0 <- utils::packageVersion("lme4") >= "2.0.0"

if(has_pkgs && is_lme4_2_0_0) {
  ## Functions:
  devfun_vp <- lmerTest:::devfun_vp
  getOptPar <- lmerTest:::getOptPar
  getVarPar <- lmerTest:::getVarPar
  
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
  
  for(model in lme4models) { # model <- lme4models[[1]]
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
  }
  
  ## Here we also want to check that devfun and and devfun_vp returns the same 
  ## value at non-optimum values of varpar.
  ## Because sigma is profiled out of devfun this cannot be done right away.
  ## We need to optimize over all parameters to get equivalence.
  ## We should be able to set one of the parameters in common to a 
  ## particular value and optimize the rest. This will build confidence that
  ## the likelihood in devfun_vp is the same as that returned by 
  ## lme4::lmer(., devFunOnly=TRUE).
  
  do_trace <- 0
  for(i in seq_along(lme4models)) { # i <- 4
    if(do_trace) print(i) 
    model <- lme4models[[i]]
    devfun <- update(model, devFunOnly=TRUE)
    (optpar <- getOptPar(model))
    (varpar <- getVarPar(model))
    (optpar2 <- optpar * 1.1)
    (varpar2 <- c(optpar2, varpar[length(varpar)]))
    
    ## Evaluate gradients to ensure that devfun and devfun_vp are both 
    ## functions with optima at optpar and varpar respectively:
    (g_devfun <- numDeriv::grad(devfun, optpar))
    (g_devfun_vp <- numDeriv::grad(devfun_vp, varpar, devfun=devfun, reml=TRUE))
    stopifnot(
      if(i != 4) all(abs(g_devfun) < 1e-3) else all(abs(g_devfun) < 1e-2),
      if(i != 4) all(abs(g_devfun_vp) < 1e-3) else all(abs(g_devfun_vp) < 1e-2)
    )
    ## These are not zero as expected:
    (g_devfun <- numDeriv::grad(devfun, optpar2))
    (g_devfun_vp <- numDeriv::grad(devfun_vp, varpar2, devfun=devfun, reml=TRUE)) 
    ## Try optimizing devfun_vp:
    x <- nlminb(start=varpar2, objective = devfun_vp, devfun=devfun, reml=TRUE,
                control=list(trace=do_trace))
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
    if(length(optpar) > 1) for(j in seq_along(optpar)) { # j <- 1
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
      }, control = list(trace=do_trace),
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
      }, control = list(trace=do_trace),
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
  }
  
}


