#############################################################################
#    Copyright (c) 2013-2020 Alexandra Kuznetsova, Per Bruun Brockhoff, and
#    Rune Haubo Bojesen Christensen
#
#    This file is part of the lmerTest package for R (*lmerTest*)
#
#    *lmerTest* is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    *lmerTest* is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at
#    <https://www.r-project.org/Licenses/> and/or
#    <http://www.gnu.org/licenses/>.
#############################################################################
#
# terms_utils.R - utilities for computing on terms objects and friends

# ------- Contents: --------
#
# --- utility functions: ---
#
# term2colX
# need_yates
# no_yates
# numeric_terms
# get_model_matrix
# get_contrast_coding
# get_min_data
# get_var_list
# get_fac_list
# get_num_list
# get_pairs
# get_trts
#

##############################################
######## term2colX()
##############################################
term2colX <- function(terms, X) {
  # Compute map from terms to columns in X using the assign attribute of X.
  # Returns a list with one element for each term containing indices of columns
  #   in X belonging to that term.
  if(is.null(asgn <- attr(X, "assign")))
    stop("Invalid design matrix:",
         "design matrix 'X' should have a non-null 'assign' attribute",
         call. = FALSE)
  term_names <- attr(terms, "term.labels")
  has_intercept <- attr(terms, "intercept") > 0
  col_terms <- if(has_intercept) c("(Intercept)", term_names)[asgn + 1] else
    term_names[asgn[asgn > 0]]
  if(!length(col_terms) == ncol(X)) # should never happen.
    stop("An error happended when mapping terms to columns of X")
  # get names of terms (including aliased terms)
  nm <- union(unique(col_terms), term_names)
  res <- lapply(setNames(as.list(nm), nm), function(x) numeric(0L))
  map <- split(seq_along(col_terms), col_terms)
  res[names(map)] <- map
  res[nm] # order appropriately
}

##############################################
######## need_yates()
##############################################
need_yates <- function(model) {
  ## Do not need yates for:
  ## - continuous variables
  ## - factors that are not contained in other factors
  ## Need yates for all other terms, i.e. terms which are:
  ##  - contained in other terms, AND
  ##  - which are not numeric/continuous
  term_names <- attr(terms(model), "term.labels")
  cont <- containment(model)
  is_contained <- names(cont[sapply(cont, function(x) length(x) > 0)])
  nmt <- numeric_terms(model)
  num_terms <- names(nmt[nmt])
  term_names[!term_names %in% num_terms &
               term_names %in% is_contained]
}

##############################################
######## no_yates()
##############################################
no_yates <- function(model) {
  setdiff(attr(terms(model), "term.labels"), need_yates(model))
}

##############################################
######## numeric_terms()
##############################################
#' @importFrom stats delete.response terms
numeric_terms <- function(model) {
  ## Determines for all terms (not just all variables) if the 'dataClass'
  ## is numeric
  ## (interactions involving one or more numerics variables are numeric).
  Terms <- delete.response(terms(model))
  all_vars <- all.vars(attr(Terms, "variables"))
  data_classes <- attr(terms(model, fixed.only=FALSE), "dataClasses")
  var_class <- data_classes[names(data_classes) %in% all_vars]
  factor_vars <- names(var_class[var_class %in% c("factor", "ordered")])
  num_vars <- setdiff(all_vars, factor_vars)

  term_names <- attr(terms(model), "term.labels")
  # term_names <- setNames(as.list(term_names), term_names)
  sapply(term_names, function(term) {
    vars <- unlist(strsplit(term, ":"))
    any(vars %in% num_vars)
  })
}

##############################################
######## get_model_matrix()
##############################################
#' Extract or remake model matrix from model
#'
#' Extract or remake model matrix from model and potentially change the
#' contrast coding
#'
#' @param model an \code{lm} or \code{lmerMod} model object.
#' @param type extract or remake model matrix?
#' @param contrasts contrasts settings. These may be restored to those in the
#' model or they may be changed. If a length one character vector (e.g.
#' \code{"contr.SAS"}) this is applied to all factors in the model, but it can
#' also be a list naming factors for which the contrasts should be set as specified.
#'
#' @return the model (or 'design') matrix.
#' @keywords internal
#' @author Rune Haubo B Christensen
get_model_matrix <- function(model, type=c("extract", "remake"),
                             contrasts="restore") {
  type <- match.arg(type)
  stopifnot(inherits(model, "lm") || inherits(model, "lmerMod"))

  if(type == "extract") return(model.matrix(model))
  # Set appropriate contrasts:
  Contrasts <- get_contrast_coding(model, contrasts=contrasts)
  model.matrix(terms(model), data=model.frame(model),
               contrasts.arg = Contrasts)
}

##############################################
######## get_contrast_coding()
##############################################
get_contrast_coding <- function(model, contrasts="restore") {
  # Compute a list of contrasts for all factors in model
  Contrasts <- contrasts
  if(length(contrasts) == 1 && is.character(contrasts) &&
     contrasts == "restore") {
    Contrasts <- attr(model.matrix(model), "contrasts")
  } else if(length(contrasts) == 1 && is.character(contrasts) &&
            contrasts != "restore") {
    Contrasts <- .getXlevels(terms(model), model.frame(model))
    Contrasts[] <- contrasts
    Contrasts
  }
  Contrasts
}


get_min_data <- function(model, FUN=mean)
  # Get a minimum complete model.frame based on the variables in the model
  do.call(expand.grid, get_var_list(model, FUN=FUN))

get_var_list <- function(model, FUN=mean)
  # Extract a named list of variables in the model containing the levels of
  # factors and the mean value of numeric variables
  c(get_fac_list(model), get_num_list(model, FUN=FUN))

#' @importFrom stats .getXlevels
get_fac_list <- function(model) {
  # Extract a named list of factor levels for each factor in the model
  res <- .getXlevels(Terms=terms(model), m=model.frame(model))
  if(is.null(res)) list() else res
}

get_num_list <- function(model, FUN=mean) { # FUN=function(x) mean(x, na.rm=TRUE)) {
  # Extract named list of mean/FUN values of numeric variables in model
  Terms <- terms(model)
  mf <- model.frame(model)
  xvars <- sapply(attr(Terms, "variables"), deparse2)[-1L]
  if((yvar <- attr(Terms, "response")) > 0)
    xvars <- xvars[-yvar]
  if(!length(xvars)) return(list())
  xlev <- lapply(mf[xvars], function(x) {
    if (is.numeric(x)) FUN(x) else NULL
  })
  res <- xlev[!vapply(xlev, is.null, NA)]
  if(is.null(res)) list() else res
}

#' @importFrom utils combn
get_pairs <- function(levs) {
  stopifnot(is.character(levs), length(levs) > 1)
  combs <- combn(seq_along(levs), 2)
  ind <- seq_len(ncombs <- ncol(combs))
  A <- as.data.frame(array(0, dim=c(length(levs), ncombs)))
  dimnames(A) <- list(levs, paste(levs[combs[1, ]], levs[combs[2, ]], sep=" - "))
  A[cbind(combs[1, ], ind)] <- 1
  A[cbind(combs[2, ], ind)] <- -1
  A
}

get_trts <- function(levs) {
  nlevs <- length(levs)
  ans <- t(cbind(-1, diag(nlevs - 1)))
  rownames(ans) <- levs
  colnames(ans) <- paste(levs[-1], levs[1], sep=" - ")
  ans
}

# get_trts(letters[1:5])
# get_pairs(letters[1:5])

