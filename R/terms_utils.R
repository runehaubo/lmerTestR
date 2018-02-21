


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
  split(seq_along(col_terms), col_terms)[unique(col_terms)]
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


get_min_data <- function(model, FUN=mean)
  # Get a minimum complete model.frame based on the variables in the model
  do.call(expand.grid, get_var_list(model, FUN=FUN))

get_var_list <- function(model, FUN=mean)
  # Extract a named list of variables in the model containing the levels of
  # factors and the mean value of numeric variables
  c(get_fac_list(model), get_num_list(model, FUN=FUN))

#' @importFrom stats .getXlevels
get_fac_list <- function(model)
  # Extract a named list of factor levels for each factor in the model
  .getXlevels(Terms=terms(model), m=model.frame(model))

get_num_list <- function(model, FUN=mean) { # FUN=function(x) mean(x, na.rm=TRUE)) {
  # Extract named list of mean/FUN values of numeric variables in model
  deparse2 <- function(x) paste(safeDeparse(x), collapse = " ")
  Terms <- terms(model)
  mf <- model.frame(model)
  xvars <- sapply(attr(Terms, "variables"), deparse2)[-1L]
  if((yvar <- attr(Terms, "response")) > 0)
    xvars <- xvars[-yvar]
  if(!length(xvars)) return(NULL)
  xlev <- lapply(mf[xvars], function(x) {
    if (is.numeric(x)) FUN(x) else NULL
  })
  xlev[!vapply(xlev, is.null, NA)]
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

