#' @importFrom mmpf uniformGrid cartesianExpand
#' @importFrom data.table data.table as.data.table rbindlist set .SD
#' @importFrom stats predict as.formula terms
#' @importFrom glmnet glmnet
#' @importFrom Matrix sparse.model.matrix
#' @importFrom utils combn
#'
#' @title generalized functional analysis of variance
#' @description computes a generalized (weighted) functional ANOVA decomposition of a prediction function, giving the best additive decomposition of the prediction function in terms of squared error.
#'
#' @param data a \code{data.frame} or \code{data.table} which contains the features/covariates on which \code{predict.fun} was learned/estimated.
#' @param vars a \code{character} vector which corresponds to a subset of the columns of \code{data}.
#' @param n a \code{numeric} vector of length 2, where the first element corresponds to the dimension of the grid constructed for each of the elements of \code{vars} and the second element to the number of rows to sample from \code{data}.
#' @param model the first argument to \code{predict.fun}, presumably a model object which can make predictions.
#' @param predict.fun a \code{function} whose first two arguments are "object" and "newdata" which returns a numeric vector the same length as the number of rows in newdata. the default value is to call the \code{predict} method on the \code{model}.
#' @param weight.fun a \code{function} with two arguments, \code{design} and \code{data}, both of which are \code{data.tables} which returns a \code{numeric} of the same length as the number of rows in \code{design}. this is intended for use to use the \code{data} to estimate the distribution of the input features, and then to use that estimate to the probability of points in the \code{design} grid.
#' @return a \code{data.table} with columns for a grid of points of the \code{vars}, a (set of) column(s) that correspond to the estimated effect of those features/covariates on the prediction function, and a column \code{effect} which indicates which subset of the covariates/features each estimate belongs to.
#'
#' @references Giles Hooker. Generalized Functional ANOVA Diagnostics for High Dimensional Functions of Dependent Variables, Journal of Computational and Graphical Statistics, Vol. 16, No. 3 (2007), pp. 709-732.
#'
#' @export
functionalANOVA = function(data, vars, n = 10, model,
  predict.fun = function(object, newdata) predict(object, newdata = newdata),
  weight.fun = NULL) {

  ## define effects to estimate
  ## u is the subset of x with indices in u
  ## estimate effects for
  ## f_u, f_{-u}, f_{sub(u)}, f_{i in u, sub(-u) U -i}
  not.vars = colnames(data)[!(colnames(data) %in% vars)]
  subset.vars = getSubsets(vars)
  effects = list(list(vars),
    if (length(unlist(subset.vars)) > 0) subset.vars else NULL,
    if (length(not.vars) > 0) list(not.vars) else NULL,
    if (length(unlist(subset.vars)) > 0) sapply(vars, function(x) c(x, not.vars), simplify = FALSE) else NULL
  )
  effects = unique(unname(unlist(effects, recursive = FALSE)))
  effects.names = sapply(effects, function(x) paste0(x, collapse = ":"))
  effects.variables = unique(unlist(effects))
  
  ## make a grid of points for the effects
  ## treat non effect variables as a unit, permute
  ## effect variable grids wrt each other and the non effect variables
  ## as a group
  points = sapply(effects.variables, function(x)
    uniformGrid(data[, x], n), simplify = FALSE)
  grid = lapply(effects, function(x)
    cartesianExpand(expand.grid(points[x]),
      as.data.table(points[!names(points) %in% x])))
  names(grid) = effects.names
  grid = rbindlist(grid, fill = TRUE, idcol = "effect")

  ## evaluate model on grid
  preds = predict.fun(model, grid[, !"effect", with = FALSE])
  ## preds = predict.fun(model, grid)
  
  ## evaluate weight function on grid
  if (!is.null(weight.fun))
    w = weight.fun(grid, data)
  else
    w = rep(1, nrow(grid))

  ## this is the simple method from the paper
  ## turn each effect into a matrix of indicators
  ## weighted least squares of matrix on F with w
  ## extract coefficients and assign to f
  
  ## create sparse indicator matrix
  formula = as.formula(paste0("~", paste0(effects.names, collapse = "+")))
  design = sparse.model.matrix(formula,
    grid[, lapply(.SD, as.factor), .SDcols = effects.variables])
  effect.idx = attributes(design)$assign
  
  ## anova fit
  fit = glmnet(design, preds, "gaussian", weights = w, lambda = 0)
  ## extract coefficients and associated variable values and
  ## form them into something useable
  betas = c(fit$a0, fit$beta[-1, 1])
  betas = lapply(unique(effect.idx), function(x)
    data.table("f" = betas[effect.idx == x]))
  names(betas) = c("intercept", attributes(terms(formula))$term.labels)
  betas = rbindlist(betas, fill = TRUE, idcol = "effect")
  
  id = strsplit(attributes(design)$Dimnames[[2]], ":")
  re = paste0("^", effects.variables, collapse = "|")
  values = sapply(id, function(x) {
    m = gregexpr(re, x)
    variables = regmatches(x, m)
    values = lapply(regmatches(x, m, TRUE), function(z) z[2])
    names(values) = variables
    as.data.table(values)
  }, simplify = FALSE)
  values = rbindlist(values, fill = TRUE)[, effects.variables, with = FALSE]
  
  ## cast all effect variables to the correct class (from the input data)
  data.types = sapply(data[, effects.variables], class)
  for (variable in names(data.types))
    set(values, j = variable,
      value = suppressWarnings(as.vector(values[[variable]], data.types[variable])))

  ## combine everything for return
  ret = cbind(betas, values)
  ret.effects = sapply(c(list(vars), getSubsets(vars)),
    function(x) paste0(x, collapse = ":"))
  ret[ret$effect %in% c(ret.effects, "intercept"),
    c("f", "effect", unique(getSubsets(vars, FALSE))), with = FALSE]
}
