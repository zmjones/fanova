#' @importFrom mmpf uniformGrid cartesianExpand
#' @importFrom data.table data.table as.data.table rbindlist set .SD setDT
#' @importFrom stats predict as.formula terms
#' @importFrom glmnet glmnet
#' @importFrom MatrixModels model.Matrix
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
functionalANOVA = function(data, vars, n = c(10, 2), model,
  predict.fun = function(object, newdata) predict(object, newdata = newdata),
  weight.fun = function(grid, data) rep(1, nrow(grid))) {

  setDT(data)
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
  ## effect variable grids wrt each other and the non effect variables as a group
  points = expand.grid(uniformGrid(data[, vars, with = FALSE], n[1]))
  if (length(not.vars) > 0) {
    not.points = uniformGrid(data[, not.vars, with = FALSE], n[2])
    grid = cartesianExpand(points, not.points)
  } else {
    grid = points
    setDT(grid)
  }

  ## for each of the vars assign a unique identifier for each value in the original grid
  ## so if x = 1, 2, 3 and z = a, b, c
  ## then grid = 1a, 2a, 3a, 1b, 2b, 3b, 1c, 2c, 3c and the ids should be
  ## 11, 21, 31, 12, 22, 32, 13, 23, 33
  ## if i made a matrix out of this where...

  ## evaluate model on grid
  preds = predict.fun(model, grid)
  
  ## evaluate weight function on grid
  w = weight.fun(grid, data)
  
  ## this is the simple method from the paper
  ## turn each effect into a matrix of indicators
  ## weighted least squares on sparse design matrix
  ## extract coefficients and assign to f

  ## create sparse indicator matrix
  formula = as.formula(paste0("~ ", paste0(effects.names, collapse = "+")))
  design = model.Matrix(formula,
    grid[, lapply(.SD, as.factor), .SDcols = effects.variables], sparse = TRUE)
  
  ## anova fit
  fit = glmnet(design, preds, "gaussian", weights = w, lambda = 0)
  ## extract coefficients and associated variable values and
  ## form them into something useable
  terms.to.extract = c(subset.vars, paste0(vars, collapse = ":"))
  terms.estimated = attributes(terms(formula))$term.labels
  max.idx = sum(sapply(terms.to.extract, function(x)
    grepl(paste0("^", x, "$"), terms.estimated)))

  effect.idx = attributes(design)$assign
  idx = effect.idx <= max.idx
  betas = fit$beta[idx, 1]
  betas[betas == 0] = NA
  betas = lapply(unique(effect.idx[idx]),
    function(x) data.table("f" = betas[effect.idx == x]))
  names(betas) = c("intercept", terms.to.extract)
  betas[["intercept"]] = data.table("f" = fit$a0)
  betas = rbindlist(betas, fill = TRUE, idcol = "effect")

  id = strsplit(attributes(design)$Dimnames[[2]][idx], ":")
  re = paste0("^", vars, collapse = "|")
  values = sapply(id, function(x) {
    m = gregexpr(re, x)
    variables = regmatches(x, m)
    values = lapply(regmatches(x, m, TRUE), function(z) z[2])
    names(values) = variables
    as.data.table(values)
  }, simplify = FALSE)
  values = rbindlist(values, fill = TRUE)[, vars, with = FALSE]
  
  ## cast all effect variables to the correct class (from the input data)
  data.types = sapply(data[, vars, with = FALSE], class)

  for (x in vars)
    set(values, j = x, value = cast(values[[x]], data.types[x]))

  ## combine everything for return
  cbind(betas, values)
}


## functionalANOVA = function(data, vars, n = c(10, 2), model,
##   predict.fun = function(object, newdata) predict(object, newdata = newdata),
##   weight.fun = function(grid, data) rep(1, nrow(grid))) {

##   setDT(data)
##   ## define effects to estimate
##   ## u is the subset of x with indices in u
##   ## estimate effects for
##   ## f_u, f_{-u}, f_{sub(u)}, f_{i in u, sub(-u) U -i}
##   not.vars = colnames(data)[!(colnames(data) %in% vars)]
##   subset.vars = getSubsets(vars)
##   effects = list(list(vars),
##     if (length(unlist(subset.vars)) > 0) subset.vars else NULL,
##     if (length(not.vars) > 0) list(not.vars) else NULL,
##     if (length(unlist(subset.vars)) > 0) sapply(vars, function(x) c(x, not.vars), simplify = FALSE) else NULL
##   )
##   effects = unique(unname(unlist(effects, recursive = FALSE)))
##   effects.names = sapply(effects, function(x) paste0(x, collapse = ":"))
##   effects.variables = unique(unlist(effects))

##   ## make a grid of points for the effects
##   ## treat non effect variables as a unit, permute
##   ## effect variable grids wrt each other and the non effect variables as a group
##   points = expand.grid(uniformGrid(data[, vars, with = FALSE], n[1]))
##   not.points = uniformGrid(data[, not.vars, with = FALSE], n[2])
##   grid = cartesianExpand(points, not.points)

##   ## evaluate model on grid
##   preds = predict.fun(model, grid)

##   ## evaluate weight function on grid
##   w = weight.fun(grid, data)
  
##   ## create sparse indicator matrix
##   main.effects = c(subset.vars, paste0(vars, collapse = ":"))
##   term.dimension = attributes(terms(formula))$order
##   N = n[1]^(length(vars)) * n[2]

##   ## design matrix
##   lapply(effects, function(x) {
##     if ()

##   })
  
  
##   x1 <- ones(n2 * n3, 1) %x% diag(n1)
##   x2 <- ones(n3, 1) %x% (diag(n2) %x% ones(n1, 1))
##   x3 <- diag(n3) %x% ones(n1 * n2, 1)
##   x12 <- ones(n3, 1) %x% diag(n1 * n2)
##   x13 <- diag(n3) %x% (ones(n2, 1) %x% diag(n1))
##   x23 <- diag(n3 * n2) %x% ones(n1, 1)
  
  

##   ## create constraint matrix
## }




## gfanova <- function(X, vars, n1, n2, n3, model) {
##   ## grid
##   idx <- 1:ncol(X) # column index
##   nidx <- idx[!idx %in% vars] # not vars column index
##   xx1 <- seq(min(X[, 1]), max(X[, 1]), length.out = n1) # points for x
##   xx2 <- seq(min(X[, 2]), max(X[, 2]), length.out = n2) # points for y
##   allpts <- cbind(ones(n2, 1) %x% xx1, xx2 %x% ones(n1, 1)) ## cartesian product of xy poitns
##   Xu <- make_unif(X, vars, n1 * n2, n3, TRUE)
##   ## make uniform design with data, vars, n[1] = n^length(vars), n[2] = n for not.vars
##   ## insert into vars slot of the result the points selected
  
##   Xu[, vars] <- ones(n3, 1) %x% allpts

##   ## predictions
##   n <- c(n1, n2, n3)
##   F <- predict(model, newdata = Xu)
##   G <- array(F, dim = n)
##   gg <- cell(3)
##   gg[[1]] <- arrayMean(arrayMean(G, 2), 3)
##   gg[[2]] <- arrayMean(arrayMean(G, 1), 3)
##   gg[[3]] <- arrayMean(G, 3)

##   ## design matrix
##   x1 <- ones(n2 * n3, 1) %x% diag(n1)
##   test <- ones(n1 * n3, 1) %x% diag(n2)
##   x2 <- ones(n3, 1) %x% (diag(n2) %x% ones(n1, 1))
##   x3 <- diag(n3) %x% ones(n1 * n2, 1)
##   x12 <- ones(n3, 1) %x% diag(n1 * n2)
##   x13 <- diag(n3) %x% (ones(n2, 1) %x% diag(n1))
##   x23 <- diag(n3 * n2) %x% ones(n1, 1)
  
##   ## intercept, univariate, bivariate
##   X <- cbind(ones(prod(n1, n2, n3), 1), x1, x2, x3, x12, x13, x23)
  
##   ## constraint matrix
##   y <- cell(6)
##   y[[1]] <- ones(1, n1)
##   y[[2]] <- ones(1, n2)
##   y[[3]] <- ones(1, n3)
##   y[[4]] <- list(diag(n2) %x% ones(1, n1), ones(1, n2) %x% diag(n1))
##   y[[5]] <- list(diag(n3) %x% ones(1, n1), ones(1, n3) %x% diag(n1))
##   y[[6]] <- list(diag(n3) %x% ones(1, n2), ones(1, n3) %x% diag(n2))

##   Y <- cbind(zeroes(3 + 2 * sum(n), 1), mattdiag_cell(y,))

  

##   Y <- zeroes(3 + 2 * sum(n), 1)
  
## }
