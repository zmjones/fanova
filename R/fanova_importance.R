#' Calculate variable importance from a fanova fit
#'
#' Computes the fraction of total variance explained by each main effect
#' and interaction term in the model.
#'
#' @param object A 'fanova' object from \code{fanova_fit}.
#' @return A data frame with the component names and their corresponding
#'   variance shares.
#' @export
fanova_importance <- function(object) {
  comps <- fanova_components(object)
  # Weighted variance of each component
  vars <- apply(comps, 2, function(c) sum(object$w * (c - sum(object$w * c))^2))
  total_var <- sum(vars[-1]) # Total variance is sum of all non-intercept components
  shares <- vars[-1] / total_var

  data.frame(component = names(shares), variance = shares)
}