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

test_that("method = \"RAW\" computes the condition number directly on X'WX, with no centering or scaling", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  res_raw <- condition_number(mod, method = "RAW")

  X <- model.matrix(mod)
  Xw <- X * sqrt(mod$weights) # X'WX = t(Xw) %*% Xw, no further transformation
  ev_manual <- sort(eigen(crossprod(Xw), symmetric = TRUE)$values, decreasing = TRUE)
  cn_manual <- max(ev_manual) / min(ev_manual)

  expect_identical(res_raw$nc_label, "NC_RAW")
  expect_identical(res_raw$method, "RAW")
  expect_equal(unname(res_raw$eigenvalues), ev_manual, tolerance = 1e-8)
  expect_equal(res_raw$condition_number, cn_manual, tolerance = 1e-8)

  # Same as calling condition_number(mod, center = FALSE, scale = FALSE) directly
  res_manual_call <- condition_number(mod, center = FALSE, scale = FALSE)
  expect_equal(res_raw$condition_number, res_manual_call$condition_number, tolerance = 1e-8)

  # Unlike NC_WS, NC_RAW is not unit-length scaled, so it is (usually) a
  # genuinely different, typically much larger number.
  res_ws <- condition_number(mod, method = "WS")
  expect_false(isTRUE(all.equal(res_raw$condition_number, res_ws$condition_number)))
})

test_that("method = \"WS\" matches center = FALSE, scale = \"unit\" (the package default) and labels the result", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  res_default <- condition_number(mod)
  res_ws <- condition_number(mod, method = "WS")

  expect_null(res_default$nc_label) # unchanged default: no method -> no label
  expect_identical(res_ws$nc_label, "NC_WS")
  expect_identical(res_ws$method, "WS")
  expect_equal(res_ws$eigenvalues, res_default$eigenvalues)
  expect_equal(res_ws$condition_number, res_default$condition_number)
})

test_that("method = \"OZ\" drops the intercept before centering/scaling, following every worked example in Ozkale (2019)", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  # Reproduce Ozkale's own procedure by hand: drop the intercept column
  # from the design matrix (keeping the IRLS weights of the original,
  # with-intercept fit), then center and scale to unit length -- this is
  # what method = "OZ" should do automatically, without requiring the user
  # to refit the model without an intercept themselves.
  X <- model.matrix(mod)
  X_no_intercept <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  Xw <- X_no_intercept * sqrt(mod$weights)
  Xc <- sweep(Xw, 2, colMeans(Xw), "-")
  X_unit <- sweep(Xc, 2, sqrt(colSums(Xc^2)), "/")
  ev_manual <- sort(eigen(crossprod(X_unit), symmetric = TRUE)$values, decreasing = TRUE)
  cn_manual <- max(ev_manual) / min(ev_manual)

  res_oz <- condition_number(mod, method = "OZ")

  expect_identical(res_oz$nc_label, "NC_OZ")
  # Because the intercept was dropped, this ordinary (non-degenerate) case
  # should NOT hit the exact-rank-deficiency warning/Inf anymore.
  expect_true(is.finite(res_oz$condition_number))
  expect_equal(unname(res_oz$condition_number), cn_manual, tolerance = 1e-6)

  # NC_OZ and NC_WS are still genuinely different definitions
  res_ws <- condition_number(mod, method = "WS")
  expect_false(isTRUE(all.equal(res_oz$condition_number, res_ws$condition_number)))
})

test_that("manual center = TRUE (bypassing method = \"OZ\") still keeps the intercept and can hit exact rank-deficiency", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  # Without going through method = "OZ", center = TRUE keeps the intercept
  # column, which for a canonical-link GLM with an intercept is
  # mathematically guaranteed to drop the rank by exactly one (see
  # bkw_diagnostics()): the smallest eigenvalue is exactly zero up to
  # floating-point noise, which can round to a tiny *negative* value.
  # condition_number() must clamp that to a proper Inf rather than silently
  # returning a nonsensical negative "condition number".
  expect_warning(
    res_manual <- condition_number(mod, center = TRUE, scale = "unit"),
    "rank-deficient"
  )
  expect_true(is.infinite(res_manual$condition_number) && res_manual$condition_number > 0)
})

test_that("method = \"MS\" refits on centered-and-unit-scaled variables and matches a manual refit", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  res_ms <- condition_number(mod, method = "MS")
  expect_identical(res_ms$nc_label, "NC_MS")
  expect_identical(res_ms$method, "MS")

  # Reproduce the Marx and Smith (1990) procedure by hand: center and
  # unit-length-scale the non-intercept columns, refit with an ordinary
  # intercept, then read the condition number directly off the refit's own
  # (unweighted-transformation) IRLS-weighted information matrix.
  X <- model.matrix(mod)
  Xvars <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  Xc <- sweep(Xvars, 2, colMeans(Xvars), "-")
  Xcu <- sweep(Xc, 2, sqrt(colSums(Xc^2)), "/")
  refit_manual <- glm(d$y ~ Xcu, family = Gamma(link = "inverse"))
  Xw_manual <- model.matrix(refit_manual) * sqrt(refit_manual$weights)
  ev_manual <- sort(eigen(crossprod(Xw_manual), symmetric = TRUE)$values, decreasing = TRUE)
  cn_manual <- max(ev_manual) / min(ev_manual)

  expect_equal(unname(res_ms$eigenvalues), ev_manual, tolerance = 1e-6)
  expect_equal(res_ms$condition_number, cn_manual, tolerance = 1e-6)

  # MS, MP, WS and OZ are all genuinely different definitions
  cn_ws <- condition_number(mod, method = "WS")$condition_number
  expect_false(isTRUE(all.equal(res_ms$condition_number, cn_ws)))
})

test_that("method = \"MS\" gives an informative error for a no-intercept Gamma/inverse-link model that fails to reconverge", {
  # Centering removes the constant baseline a no-intercept model relies on;
  # for the Gamma family with the inverse link (which needs positive fitted
  # values), the refit's IRLS iteration can diverge entirely. This is a
  # genuine limitation of the method for this family/link combination (Marx
  # and Smith 1990 developed it for logistic regression with an intercept),
  # not a bug -- the error should be informative rather than a bare glm()
  # failure.
  d <- data.frame(
    y  = c(6, 5, 5, 3, 7, 9, 6, 2, 10, 7, 3, 4, 13, 10, 7, 3, 6, 5, 4, 9, 11, 8, 9, 6, 2),
    x1 = c(11.1, 12.1, 12, 17.8, 9.5, 7.2, 11.5, 13.4, 10.8, 13.8, 14.6, 12.1, 8, 8.8, 12.9,
           12.7, 12.1, 11.1, 11.3, 9, 9.2, 8.4, 8, 13.8, 17.8),
    x2 = c(90, 86, 80, 70, 90, 100, 92, 74, 87, 78, 73, 85, 94, 91, 84, 68, 81, 78, 74, 78,
           84, 90, 90, 80, 68),
    x3 = c(382, 380, 372, 352, 358, 362, 302, 316, 339, 328, 278, 339, 241, 193, 268, 113,
           313, 317, 324, 312, 349, 290, 295, 283, 259),
    x4 = c(12, 20, 19, 16, 10, 12, 15, 15, 14, 14, 5, 17, 16, 13, 8, -9, 6, 10, 1, 5, 4, 14,
           9, 5, -10)
  )
  mod <- glm(y ~ x1 + x2 + x3 + x4 - 1, family = Gamma(link = "inverse"), data = d)
  expect_error(condition_number(mod, method = "MS"), "method = \"MS\"")
})

test_that("method = \"MS\" works for a no-intercept binomial/logit model (Marx and Smith's own setting)", {
  # Lee (1974) cancer remission data, fit without an intercept as in Ozkale
  # (2019, Sec. 5.2) -- unlike the Gamma/inverse case above, binomial/logit
  # has no structural obstruction to centering a no-intercept model, so
  # method = "MS" should converge to a finite, well-defined result here.
  d <- data.frame(
    y  = c(1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0),
    x1 = c(0.8, 0.9, 0.8, 1, 0.9, 1, 0.95, 0.95, 1, 0.95, 0.85, 0.7, 0.8, 0.2, 1, 1, 0.65, 1,
           0.5, 1, 1, 0.9, 1, 0.95, 1, 1, 1),
    x2 = c(0.83, 0.36, 0.88, 0.87, 0.75, 0.65, 0.97, 0.87, 0.45, 0.36, 0.39, 0.76, 0.46, 0.39,
           0.9, 0.84, 0.42, 0.75, 0.44, 0.63, 0.33, 0.93, 0.58, 0.32, 0.6, 0.69, 0.73),
    x3 = c(0.66, 0.32, 0.7, 0.87, 0.68, 0.65, 0.92, 0.83, 0.45, 0.34, 0.33, 0.53, 0.37, 0.08,
           0.9, 0.84, 0.27, 0.75, 0.22, 0.63, 0.33, 0.84, 0.58, 0.3, 0.6, 0.69, 0.73),
    x5 = c(1.1, 0.74, 0.176, 1.053, 0.519, 0.519, 1.23, 1.354, 0.322, 0, 0.279, 0.146, 0.38,
           0.114, 1.037, 2.064, 0.114, 1.322, 0.114, 1.072, 0.176, 1.591, 0.531, 0.886, 0.964,
           0.398, 0.398),
    x6 = c(0.996, 0.992, 0.982, 0.986, 0.98, 0.982, 0.992, 1.02, 0.999, 1.038, 0.988, 0.982,
           1.006, 0.99, 0.99, 1.02, 1.014, 1.004, 0.99, 0.986, 1.01, 1.02, 1.002, 0.988, 0.99,
           0.986, 0.986)
  )
  mod <- glm(y ~ x1 + x2 + x3 + x5 + x6 - 1, family = binomial(), data = d)
  res_ms <- condition_number(mod, method = "MS")
  expect_true(is.finite(res_ms$condition_number))
  expect_identical(res_ms$nc_label, "NC_MS")
})

test_that("method = \"MP\" refits on unit-length-rescaled variables and matches the direct (no-refit) formula", {
  skip_if_not_installed("multiColl")

  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)

  res_mp <- condition_number(mod, method = "MP")
  expect_identical(res_mp$nc_label, "NC_MP")

  # Direct formula, no refit needed: a GLM's fitted values -- and hence its
  # IRLS weights at convergence -- are invariant to any linear rescaling of
  # the columns of X, so MP can equivalently be computed by scaling the
  # ORIGINAL X to unit length first and weighting with mod's own IRLS
  # weights afterwards, with no further transformation.
  X <- model.matrix(mod)
  X_unit_first <- multiColl::lu(X)
  Xw <- X_unit_first * sqrt(mod$weights)
  M <- crossprod(Xw)
  ev <- sort(eigen(M, symmetric = TRUE)$values, decreasing = TRUE)
  cn_direct <- max(ev) / min(ev)

  expect_equal(unname(res_mp$condition_number), cn_direct, tolerance = 1e-6)

  # MP, WS and OZ are three genuinely different definitions
  cn_ws <- condition_number(mod, method = "WS")$condition_number
  expect_false(isTRUE(all.equal(res_mp$condition_number, cn_ws)))
})

test_that("method = \"MP\" gives an informative error when 'multiColl' is not installed", {
  skip_if(requireNamespace("multiColl", quietly = TRUE), "multiColl is installed; cannot test the missing-package branch")

  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  expect_error(condition_number(mod, method = "MP"), "multiColl")
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

test_that("method = \"OZ\" can still be exactly rank-deficient for a Gamma/inverse-link model with NO intercept (Ozkale 2019, Sec. 5.3)", {
  # Nitrogen dioxide data (Chatterjee & Hadi 1998), used by Ozkale (2019,
  # Sec. 5.3) with a Gamma/inverse-link GLM fit WITHOUT an intercept. This
  # is the case that method = "OZ"'s intercept-drop cannot rescue: the
  # structural identity sqrt(w) * linear.predictors == 1 holds for the
  # Gamma-inverse link regardless of whether an intercept was ever present,
  # so the (centered, unit-scaled) design matrix is exactly rank-deficient
  # here too.
  d <- data.frame(
    y  = c(6, 5, 5, 3, 7, 9, 6, 2, 10, 7, 3, 4, 13, 10, 7, 3, 6, 5, 4, 9, 11, 8, 9, 6, 2),
    x1 = c(11.1, 12.1, 12, 17.8, 9.5, 7.2, 11.5, 13.4, 10.8, 13.8, 14.6, 12.1, 8, 8.8, 12.9,
           12.7, 12.1, 11.1, 11.3, 9, 9.2, 8.4, 8, 13.8, 17.8),
    x2 = c(90, 86, 80, 70, 90, 100, 92, 74, 87, 78, 73, 85, 94, 91, 84, 68, 81, 78, 74, 78,
           84, 90, 90, 80, 68),
    x3 = c(382, 380, 372, 352, 358, 362, 302, 316, 339, 328, 278, 339, 241, 193, 268, 113,
           313, 317, 324, 312, 349, 290, 295, 283, 259),
    x4 = c(12, 20, 19, 16, 10, 12, 15, 15, 14, 14, 5, 17, 16, 13, 8, -9, 6, 10, 1, 5, 4, 14,
           9, 5, -10)
  )
  mod <- glm(y ~ x1 + x2 + x3 + x4 - 1, family = Gamma(link = "inverse"), data = d)

  eta_hat <- as.numeric(model.matrix(mod) %*% coef(mod))
  # unname(): mod$weights carries the data's row names, which would
  # otherwise make expect_equal() report a spurious mismatch (differing
  # names attribute) even though every value matches within tolerance.
  expect_equal(unname(sqrt(mod$weights) * eta_hat), rep(1, nrow(d)),
               tolerance = 1e-6)

  expect_warning(res_oz <- condition_number(mod, method = "OZ"), "rank-deficient")
  expect_true(is.infinite(res_oz$condition_number) && res_oz$condition_number > 0)

  # method = "WS" (the package default) has no centering step, so it is
  # unaffected by this degeneracy and remains finite.
  res_ws <- condition_number(mod, method = "WS")
  expect_true(is.finite(res_ws$condition_number))
})

test_that("method = \"RAW\" reproduces the mine fracture data condition indices (Myers 1990; Marx 1992; Kurtoglu and Ozkale 2017)", {
  # Mine fracture data (Myers 1990), 44 coal mine areas in the Appalachian
  # region of western Virginia, used by Marx (1992) and Kurtoglu and Ozkale
  # (2017) to illustrate severe multicollinearity in a Poisson GLM. y is the
  # number of upper seam injuries/fractures; x1 = inner burden thickness,
  # x2 = percent extraction of the lower previously mined seam, x3 = lower
  # seam height, x4 = time (years) the mine has been open.
  d <- data.frame(
    y  = c(2, 1, 0, 4, 1, 2, 0, 0, 4, 4, 1, 4, 1, 5, 2, 5, 5, 5, 0, 5, 1, 1, 3, 3, 2, 2, 0, 1,
           5, 2, 3, 3, 3, 0, 0, 2, 0, 0, 3, 2, 3, 5, 0, 3),
    x1 = c(50, 230, 125, 75, 70, 65, 65, 350, 350, 160, 145, 145, 180, 43, 42, 42, 45, 83, 300,
           190, 145, 510, 65, 470, 300, 275, 420, 65, 40, 900, 95, 40, 140, 150, 80, 80, 145,
           100, 150, 150, 210, 11, 100, 50),
    x2 = c(70, 65, 70, 65, 65, 70, 60, 60, 90, 80, 65, 85, 70, 80, 85, 85, 85, 85, 65, 90, 90,
           80, 75, 90, 80, 90, 50, 80, 75, 90, 88, 85, 90, 50, 60, 85, 65, 65, 80, 80, 75, 75,
           65, 88),
    x3 = c(52, 42, 45, 68, 53, 46, 62, 54, 54, 38, 38, 38, 42, 40, 51, 51, 42, 48, 68, 84, 54,
           57, 68, 90, 165, 40, 44, 48, 51, 48, 36, 57, 38, 44, 96, 96, 72, 72, 48, 48, 42, 42,
           60, 60),
    x4 = c(1, 6, 1, 0.5, 0.5, 3, 1, 0.5, 0.5, 0, 10, 0, 2, 0, 12, 0, 0, 10, 10, 6, 12, 10, 5, 9,
           9, 4, 17, 15, 15, 35, 20, 10, 7, 5, 5, 5, 9, 9, 3, 0, 2, 0, 25, 20)
  )

  # The published OLS estimator (used as the IRLS starting value in the
  # papers above) reproduces exactly, which confirms the data set is
  # transcribed correctly (two of its 44 rows are easy to mistype: row 30
  # has y = 2, not 22, and row 36 has x1 = 80, not 0).
  ols <- lm(y ~ x1 + x2 + x3 + x4, data = d)
  expect_equal(
    unname(coef(ols)),
    c(-4.579885, -0.001896, 0.104574, -0.007993, -0.050250),
    tolerance = 1e-5
  )

  mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(link = "log"), data = d)
  res_raw <- condition_number(mod, method = "RAW")

  expect_identical(res_raw$nc_label, "NC_RAW")
  # Published eigenvalues of X'WX (Marx 1992): 4164686, 338966.8, 30607.81,
  # 3824.866, 0.9504078 -- matched to within IRLS convergence tolerance.
  expect_equal(
    unname(res_raw$eigenvalues),
    c(4164686, 338966.8, 30607.81, 3824.866, 0.9504078),
    tolerance = 1e-4
  )
  # Published condition indices: 1.0000, 3.5052, 11.6647, 32.9976, 2093.3225
  condition_indices <- sqrt(max(res_raw$eigenvalues) / res_raw$eigenvalues)
  expect_equal(
    condition_indices,
    c(1.0000, 3.5052, 11.6647, 32.9976, 2093.3225),
    tolerance = 1e-3
  )
  expect_equal(res_raw$condition_index, 2093.3225, tolerance = 1e-3)

  # Ridge parameter k = p / (b'b) with p = number of coefficients INCLUDING
  # the intercept (5, not 4) reproduces the published k = 0.3871427.
  k <- length(coef(mod)) / sum(coef(mod)^2)
  expect_equal(k, 0.3871427, tolerance = 1e-4)
})
