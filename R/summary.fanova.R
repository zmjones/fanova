#' Summarize a fanova fit
#'
#' Provides a summary of the fitted functional ANOVA model, including the
#' variance decomposition.
#'
#' @param object A 'fanova' object.
#' @param ... Other arguments (unused).
#'
#' @return An object of class 'summary.fanova', which is a list containing
#'   the model call, dimensions, and the importance data frame.
#' @export
#' @method summary fanova
summary.fanova <- function(object, ...) {
  imp <- fanova_importance(object)

  info <- list(
    call = object$call,
    n_obs = nrow(object$X),
    n_vars = ncol(object$X),
    max_order = object$max_order,
    importance = imp
  )

  class(info) <- "summary.fanova"
  return(info)
}

#' @export
#' @method print summary.fanova
#' @rdname summary.fanova
print.summary.fanova <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(sprintf("fANOVA model on %d observations and %d variables up to order %d.\n\n", x$n_obs, x$n_vars, x$max_order))
  cat("Variance Decomposition:\n")
  print(x$importance)
  invisible(x)
}