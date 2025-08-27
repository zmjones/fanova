# functionalANOVA

Empirical implementation of Hooker (2007) generalized functional ANOVA for dependent inputs.

```r
devtools::load_all(".")
set.seed(1); n <- 400; x1 <- rnorm(n); x2 <- 0.7*x1 + rnorm(n, sd=.7)
ftrue <- function(x) x[,1]^2 + x[,2] + 1.5*x[,1]*x[,2]
X <- cbind(x1, x2); y <- ftrue(X)
fit <- fanova_fit(X, y, max_order=2)
fanova_importance(fit)
plot(fit)
