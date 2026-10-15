simulated_data <- function(seed = 1, n = 200) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.05) # correlated with x1 -> collinearity
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  # large intercept, small slopes so eta = 1/mu stays well above 0
  eta <- 3 + 0.15 * x1 + 0.15 * x2 - 0.1 * x3 + 0.05 * x4
  mu <- 1 / eta
  y <- rgamma(n, shape = 5, rate = 5 / mu)
  data.frame(y, x1, x2, x3, x4)
}

# Condition index (sqrt scale) of a matrix, from its singular values
sv_ci <- function(M) {
  d <- svd(M)$d
  max(d) / min(d)
}

unit_length <- function(M) sweep(M, 2L, sqrt(colSums(M^2)), "/")
