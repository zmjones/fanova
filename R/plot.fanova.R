#' Plot fanova importance
#'
#' Creates a bar plot of the variance shares for each component in the
#' fANOVA decomposition.
#'
#' @param x A 'fanova' object.
#' @param ... Other arguments passed to \code{barplot}.
#'
#' @return Invisibly returns the result of \code{barplot}.
#' @importFrom graphics barplot par
#' @export
#' @method plot fanova
plot.fanova <- function(x, ...) {
  imp <- fanova_importance(x)
  imp <- imp[order(imp$variance), ]

  par(mar = c(5, 6, 4, 2) + 0.1)
  b <- barplot(imp$variance, names.arg = imp$component,
               main = "fANOVA Variance Shares", horiz = TRUE, las = 1, ...)
  invisible(b)
}