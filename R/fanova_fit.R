#' Fit a Generalized Functional ANOVA model
#'
#' Decomposes a function f(x) into a hierarchically orthogonal, ANOVA-like
#' representation of main effects and interactions.
#'
#' @param X A numeric matrix or data frame of features (n x p).
#' @param f A numeric vector of responses, or a function that takes X and returns a numeric vector.
#' @param weights Optional numeric vector of weights for each observation.
#' @param max_order The maximum interaction order to include (e.g., 2 for main effects and pairwise interactions).
#' @param df1 The number of spline basis functions for each univariate main effect.
#' @param degree The degree of the splines.
#' @param include An optional list of integer vectors specifying the exact interaction sets to include.
#' @return An object of class 'fanova'.
#' @export
#' @importFrom utils combn
fanova_fit <- function(X, f, weights = NULL, max_order = 2L, df1 = 6L, degree = 3L, include = NULL) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (is.null(weights)) weights <- rep(1, n)
  w <- as.numeric(weights) / sum(weights)

  # Evaluate f if function
  y <- if (is.function(f)) as.numeric(f(X)) else as.numeric(f)
  stopifnot(length(y) == n)

  # Intercept
  intercept <- rep(1, n)
  mu <- sum(w * y)
  y_center <- y - mu  # remove f_∅

  # Build univariate bases (centered => mean zero under w)
  uni <- vector("list", p)
  for (j in seq_len(p)) uni[[j]] <- build_basis_1d(X[, j], df = df1, degree = degree, w = w)

  # Assemble list of sets to include
  if (is.null(include)) {
    include <- list()
    # Main effects are built separately. This loop generates interaction terms.
    if (max_order > 1L) {
      for (m in 2:max_order) {
        include <- c(include, combn(p, m, simplify = FALSE))
      }
    }
  }

  # Lower-order spaces by index
  lower_span <- list(`0` = matrix(intercept, n, 1))
  for (j in seq_len(p)) lower_span[[paste0("1_", j)]] <- uni[[j]]

  # Build design matrices respecting hierarchical orthogonality:
  # - Univariates: already mean-zero wrt intercept
  # - Pairs: build tensor and residualize vs intercept + both relevant univariates
  blocks <- list()
  group_names <- character()

  # order 1
  for (j in seq_len(p)) {
    blocks <- c(blocks, list(uni[[j]]))
    group_names <- c(group_names, paste0("1_", j))
  }

  # order 2
  if (max_order >= 2L) {
    for (S in include) if (length(S) == 2L) {
      i <- S[1]; j <- S[2]
      Bij_raw <- tensor2(uni[[i]], uni[[j]])
      # residualize vs intercept, Bi, Bj
      lower <- cbind(intercept, uni[[i]], uni[[j]])
      Bij <- w_residualize(Bij_raw, lower, w)
      # Drop any all-zero columns that result from projection
      keep <- which(colSums(abs(Bij)) > 1e-10)
      if (length(keep) > 0) {
        blocks <- c(blocks, list(Bij[, keep, drop = FALSE]))
        group_names <- c(group_names, paste0("2_", i, "_", j))
      }
    }
  }

  # (Higher orders can be added similarly by residualizing against all proper-subset spans.)

  # Final orthogonal design (across orders) = cbind of blocks
  Xdesign <- do.call(cbind, blocks)

  # Solve WLS
  beta <- wls_solve(Xdesign, y_center, w)

  # Extract group fits
  fitted_groups <- list()
  offset <- 0L
  for (g in seq_along(blocks)) {
    k <- ncol(blocks[[g]])
    b <- beta[(offset + 1):(offset + k)]
    fitted_groups[[group_names[g]]] <- drop(blocks[[g]] %*% b)
    offset <- offset + k
  }

  structure(list(
    call = match.call(),
    X = X,
    w = w,
    y = y,
    mu = mu,
    basis = uni,
    blocks = blocks,
    beta = beta,
    groups = fitted_groups,
    group_names = group_names,
    df1 = df1,
    degree = degree,
    max_order = max_order
  ), class = "fanova")
}

#' (Internal) Extract component functions
#'
#' @param object A fanova object.
#' @return A matrix where each column is an evaluated component function.
fanova_components <- function(object) {
  comps <- cbind(`f0` = rep(object$mu, length(object$y)))
  for (nm in object$group_names) {
    comps <- cbind(comps, object$groups[[nm]])
    colnames(comps)[ncol(comps)] <- nm
  }
  comps
}
