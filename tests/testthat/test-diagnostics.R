simulated_data <- function(seed = 1, n = 200) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.05) # correlated with x1 -> collinearity
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  # large intercept, small slopes so eta = 1/mu stays well above 0
  # (mu > 0 always) despite the x1/x2 correlation
  eta <- 3 + 0.15 * x1 + 0.15 * x2 - 0.1 * x3 + 0.05 * x4
  mu <- 1 / eta
  y <- rgamma(n, shape = 5, rate = 5 / mu)
  data.frame(y, x1, x2, x3, x4)
}

test_that("condition_number.glm reproduces the original manual calculation", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  # --- manual calculation, exactly as in the original code ---
  X <- model.matrix(mod)
  W <- diag(mod$weights)
  Wsqrt <- diag(sqrt(diag(W)))
  Xw <- Wsqrt %*% X
  X_unit <- scale(Xw, center = FALSE, scale = sqrt(colSums(Xw^2)))
  M_star <- t(X_unit) %*% X_unit
  eig_star <- sort(eigen(M_star, symmetric = TRUE)$values, decreasing = TRUE)
  CN_manual <- max(eig_star) / min(eig_star)

  res <- condition_number(mod)

  expect_s3_class(res, "multicollglm_cn")
  expect_equal(unname(res$eigenvalues), eig_star, tolerance = 1e-8)
  expect_equal(res$condition_number, CN_manual, tolerance = 1e-8)
  expect_equal(res$condition_index, sqrt(CN_manual), tolerance = 1e-8)
})

test_that("bkw_diagnostics.glm produces proportions that sum to 1 per row", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  res <- bkw_diagnostics(mod)

  expect_s3_class(res, "multicollglm_bkw")
  expect_equal(unname(rowSums(res$proportions)), rep(1, nrow(res$proportions)), tolerance = 1e-8)
  expect_equal(res$condition_index[1], 1, tolerance = 1e-8)
  expect_equal(max(res$condition_index), sqrt(max(res$eigenvalues) / min(res$eigenvalues)), tolerance = 1e-8)
  # x1, strongly collinear with x2, should show up in a flagged dimension
  flagged_vars <- unlist(lapply(res$flagged, `[[`, "variables"))
  expect_true(any(c("x1", "x2") %in% flagged_vars))
})

test_that("rvif_diagnostics.glm matches a direct call to rvif::rvifs() on X_unit", {
  skip_if_not_installed("rvif")

  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  res <- rvif_diagnostics(mod)

  X <- model.matrix(mod)
  Xw <- X * sqrt(mod$weights)
  X_unit <- scale(Xw, center = FALSE, scale = sqrt(colSums(Xw^2)))
  expected <- rvif::rvifs(X_unit, ul = TRUE, intercept = TRUE)

  expect_s3_class(res, "multicollglm_rvif")
  expect_equal(unname(res$table$RVIF), unname(expected$RVIF), tolerance = 1e-8)
  expect_equal(unname(res$table$`%`), unname(expected$`%`), tolerance = 1e-8)
  expect_equal(rownames(res$table), colnames(X))
  # x1/x2 should carry most of the near-collinearity, given how they were simulated
  expect_true(all(res$table[c("x1", "x2"), "%"] > res$table["x3", "%"]))
})

test_that("scale = TRUE (root-mean-square) preserves the condition number ratio but rescales eigenvalues", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  res_unit <- condition_number(mod) # default: center = FALSE, scale = "unit"
  res_rms <- condition_number(mod, center = FALSE, scale = TRUE)

  n <- nrow(model.matrix(mod))
  # scale = TRUE divides by sqrt(colSums(x^2) / (n - 1)) instead of
  # sqrt(colSums(x^2)), so eigenvalues differ by a factor of (n - 1)
  expect_equal(res_rms$eigenvalues, res_unit$eigenvalues * (n - 1), tolerance = 1e-8)
  # the ratio max/min eigenvalue, and therefore the condition number, is unaffected
  expect_equal(res_rms$condition_number, res_unit$condition_number, tolerance = 1e-8)
})

test_that("scale = FALSE and a numeric 'scale' vector are delegated to base::scale()", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  X <- model.matrix(mod)
  Xw <- X * sqrt(mod$weights)

  res_false <- condition_number(mod, center = FALSE, scale = FALSE)
  M_false <- crossprod(Xw)
  eig_false <- sort(eigen(M_false, symmetric = TRUE)$values, decreasing = TRUE)
  expect_equal(unname(res_false$eigenvalues), eig_false, tolerance = 1e-6)

  divisors <- c(2, 3, 4, 5, 6)
  res_num <- condition_number(mod, center = FALSE, scale = divisors)
  X_num <- scale(Xw, center = FALSE, scale = divisors)
  M_num <- crossprod(X_num)
  eig_num <- sort(eigen(M_num, symmetric = TRUE)$values, decreasing = TRUE)
  expect_equal(unname(res_num$eigenvalues), eig_num, tolerance = 1e-6)
})

test_that("center = TRUE on an intercept model raises an informative error instead of returning NaN", {
  # For a canonical-link GLM (e.g. Gamma with the inverse link), the
  # IRLS-weighted design matrix satisfies Xw %*% coef(mod) == 1 exactly
  # (since sqrt(w) = mu = 1/eta, so sqrt(w) * eta = 1), so centering an
  # intercept column always drops the rank by exactly one. This is the
  # classical reason Belsley-Kuh-Welsch recommend against centering.
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  X <- model.matrix(mod)
  Xw <- X * sqrt(mod$weights)
  Xc <- sweep(Xw, 2L, colMeans(Xw), "-")
  expect_equal(qr(Xc)$rank, ncol(Xc) - 1L)

  expect_error(bkw_diagnostics(mod, center = TRUE, scale = "unit"), "rank-deficient")
  expect_error(rvif_diagnostics(mod, center = TRUE, scale = "unit"), "rvifs")
})

test_that("condition_number.glmnet refits on the active set and matches the equivalent glm", {
  skip_if_not_installed("glmnet")
  library(glmnet)

  d <- simulated_data()
  x <- as.matrix(d[, c("x1", "x2", "x3", "x4")])
  y <- d$y

  fit <- glmnet(x, y, family = Gamma(link = "inverse"), lambda = c(0.05, 0.01, 0.001))
  chosen_s <- fit$lambda[2]

  res <- condition_number(fit, x = x, y = y, family = Gamma(link = "inverse"), s = chosen_s)

  # reproduce the refit by hand to compare
  cf <- coef(fit, s = chosen_s)
  active <- rownames(cf)[as.numeric(cf) != 0]
  active <- setdiff(active, "(Intercept)")
  form <- reformulate(active, response = "y")
  manual_mod <- glm(form, data = d, family = Gamma(link = "inverse"))
  expected <- condition_number(manual_mod)

  expect_s3_class(res, "multicollglm_cn")
  expect_equal(res$eigenvalues, expected$eigenvalues, tolerance = 1e-6)
  expect_equal(res$condition_number, expected$condition_number, tolerance = 1e-6)
  expect_equal(sort(res$active_vars), sort(active))
})

test_that("condition_number.glmnet requires 's' and a valid 'family'", {
  skip_if_not_installed("glmnet")
  library(glmnet)

  d <- simulated_data()
  x <- as.matrix(d[, c("x1", "x2", "x3", "x4")])
  y <- d$y
  fit <- glmnet(x, y, family = Gamma(link = "inverse"), lambda = c(0.05, 0.01))

  expect_error(condition_number(fit, x = x, y = y, family = Gamma(link = "inverse")), "lambda")
  expect_error(condition_number(fit, x = x, y = y, s = fit$lambda[1]), "family")
})
