#' Plot fanova model components
#'
#' Creates plots for a fitted fanova object. By default, it shows a bar plot
#' of variance shares. If the `which` argument is specified, it plots the
#' shape of an individual component function.
#'
#' @param x A 'fanova' object.
#' @param which The component to plot. If `NULL` (the default), a bar plot of
#'   variance shares is shown. Otherwise, `which` should be the name of a
#'   component (e.g., "1_1" or "2_1_2").
#' @param ... Other arguments passed to the underlying plot function.
#'
#' @return Invisibly returns the plotted data or the result of `barplot`.
#' @importFrom graphics barplot par rug
#' @importFrom grDevices colorRampPalette
#' @export
#' @method plot fanova
plot.fanova <- function(x, which = NULL, ...) {
  if (is.null(which)) {
    # Original importance bar plot
    imp <- fanova_importance(x)
    imp <- imp[order(imp$variance), ]

    par(mar = c(5, 6, 4, 2) + 0.1)
    b <- barplot(imp$variance, names.arg = imp$component,
                 main = "fANOVA Variance Shares", horiz = TRUE, las = 1, ...)
    return(invisible(b))
  }

  # Plotting a specific component
  if (!is.character(which) || length(which) != 1 || !(which %in% x$group_names)) {
    stop("'which' must be a valid component name string: ", paste(x$group_names, collapse = ", "))
  }

  parts <- as.integer(strsplit(which, "_")[[1]])
  order <- parts[1]
  vars <- parts[-1]

  if (order == 1) {
    # 1D main effect plot
    j <- vars[1]
    plot_data <- data.frame(x = x$X[, j], y = x$groups[[which]])
    plot_data <- plot_data[order(plot_data$x), ]

    plot(plot_data$x, plot_data$y, type = 'l',
         xlab = paste0("x", j), ylab = paste0("f(", which, ")"),
         main = paste("Main Effect:", which), ...)
    rug(x$X[, j])
    return(invisible(plot_data))

  } else if (order == 2) {
    # 2D interaction plot
    i <- vars[1]; j <- vars[2]
    f_ij <- x$groups[[which]]

    # Use color to represent the value of the interaction component
    pal <- grDevices::colorRampPalette(c("blue", "white", "red"))
    colors <- pal(100)[cut(f_ij, 100, labels = FALSE)]

    plot(x$X[, i], x$X[, j], col = colors, pch = 16,
         xlab = paste0("x", i), ylab = paste0("x", j),
         main = paste("Interaction:", which), ...)
    return(invisible(data.frame(x1 = x$X[, i], x2 = x$X[, j], f_ij = f_ij)))

  } else {
    stop("Plotting for order > 2 is not implemented.")
  }
}