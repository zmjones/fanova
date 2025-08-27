test_that("fanova recovers simple structure", {
  set.seed(1); n <- 300
  # Use independent inputs to make effects separable for this unit test
  x1 <- rnorm(n); x2 <- rnorm(n)
  # Function is dominated by main effect of x1 and interaction
  ftrue <- function(X) 2*X[,1]^2 + 0.1*X[,2] + 1.2*X[,1]*X[,2]
  X <- cbind(x1, x2); y <- ftrue(X) + rnorm(n, .01)
  fit <- fanova_fit(X, y, max_order=2)
  imps <- fanova_importance(fit)

  # Check that the class is correct
  expect_s3_class(fit, "fanova")

  # Check that importance is a data frame
  expect_s3_class(imps, "data.frame")

  # Check that the most important components are x1 and the interaction
  imps_ordered <- imps[order(imps$variance, decreasing = TRUE), ]
  top_two <- imps_ordered$component[1:2]
  expect_true(all(c("1_1", "2_1_2") %in% top_two))
})

test_that("importance scores are valid", {
  set.seed(2); n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  y <- X[,1] + X[,2]^2 + rnorm(n, 0.1)
  fit <- fanova_fit(X, y, max_order = 2)
  imps <- fanova_importance(fit)

  expect_s3_class(imps, "data.frame")
  expect_named(imps, c("component", "variance"))
  expect_true(all(imps$variance >= 0))
  expect_equal(sum(imps$variance), 1.0, tolerance = 1e-9)
})

test_that("additive models (max_order = 1) work", {
  set.seed(3); n <- 100
  X <- matrix(rnorm(n * 3), ncol = 3)
  y <- X[,1] + X[,2] + rnorm(n, 0.1)
  fit <- fanova_fit(X, y, max_order = 1)
  imps <- fanova_importance(fit)

  # Should only have main effect components
  expect_equal(nrow(imps), 3)
  expect_true(all(grepl("^1_", imps$component)))
})

test_that("fitting with a function for 'f' works", {
  set.seed(4); n <- 50
  X <- matrix(rnorm(n * 2), ncol = 2)
  my_func <- function(x_mat) x_mat[,1] * x_mat[,2]

  fit <- fanova_fit(X, f = my_func, max_order = 2)
  expect_s3_class(fit, "fanova")
  expect_equal(fit$y, my_func(X))
})
