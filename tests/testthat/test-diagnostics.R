test_that("each method matches a manual calculation of its matrix", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  X <- model.matrix(mod)
  Xw <- X * sqrt(mod$weights)

  # WS: unit-length columns of W^(1/2) X (= S X'WX S with S = diag(X'WX)^(-1/2))
  res_ws <- condition_number(mod) # the default
  expect_identical(res_ws$method, "WS")
  expect_identical(res_ws$nc_label, "NC_WS")
  S <- diag(1 / sqrt(diag(crossprod(Xw))))
  ev_ws <- eigen(S %*% crossprod(Xw) %*% S, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(unname(res_ws$eigenvalues), ev_ws, tolerance = 1e-8)
  expect_equal(res_ws$condition_index, sv_ci(unit_length(Xw)), tolerance = 1e-8)
  expect_equal(res_ws$condition_number, res_ws$condition_index^2, tolerance = 1e-8)

  # X: unit-length columns of X, no weights
  res_x <- condition_number(mod, method = "X")
  expect_equal(res_x$condition_index, sv_ci(unit_length(X)), tolerance = 1e-8)

  # RAW: X'WX with no transformation
  res_raw <- condition_number(mod, method = "RAW")
  expect_equal(res_raw$condition_index, sv_ci(Xw), tolerance = 1e-8)

  # MP: scale X first, refit, then X_ul' W X_ul
  Xul <- unit_length(X)
  refit <- glm.fit(Xul, mod$y, family = mod$family, intercept = FALSE)
  res_mp <- condition_number(mod, method = "MP")
  expect_equal(res_mp$condition_index, sv_ci(Xul * sqrt(refit$weights)), tolerance = 1e-5)

  # the four definitions are genuinely different
  cis <- sapply(c("WS", "MP", "X", "RAW"), function(m) condition_number(mod, method = m)$condition_index)
  expect_equal(length(unique(round(cis, 6))), 4L)
})

test_that("the default is WS and 'method' is validated", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  expect_equal(condition_number(mod)$condition_number, condition_number(mod, method = "WS")$condition_number)
  expect_equal(bkw_diagnostics(mod)$condition_index, bkw_diagnostics(mod, method = "WS")$condition_index)
  expect_error(condition_number(mod, method = "nope"))
})

test_that("removed methods and arguments fail with an informative message", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  expect_error(condition_number(mod, method = "OZ"), "removed in multiCollglm 0.4.0")
  expect_error(condition_number(mod, method = "MS"), "removed in multiCollglm 0.4.0")
  expect_error(bkw_diagnostics(mod, method = "OZ"), "removed")
  expect_error(rvif_diagnostics(mod, method = "MS"), "removed")
  expect_error(condition_number(mod, center = TRUE), "center")
  expect_error(condition_number(mod, scale = TRUE), "scale")
  expect_error(bkw_diagnostics(mod, center = TRUE), "center")
  expect_error(rvif_diagnostics(mod, scale = "unit"), "scale")
})

test_that("with equal weights (Gaussian, identity) WS, MP, X coincide and match a linear-model diagnostic", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = gaussian(), data = d)
  cis <- sapply(c("WS", "MP", "X"), function(m) condition_number(mod, method = m)$condition_index)
  expect_equal(unname(cis[1]), unname(cis[3]), tolerance = 1e-8)
  expect_equal(unname(cis[2]), unname(cis[3]), tolerance = 1e-6)
  # Belsley-Kuh-Welsch on scale(X, center = FALSE)
  Xs <- scale(model.matrix(mod), center = FALSE, scale = sqrt(colSums(model.matrix(mod)^2)))
  expect_equal(unname(cis[3]), sv_ci(Xs), tolerance = 1e-8)
})

test_that("under orthogonality in the GLM metric WS attains its minimum, 1, for unequal weights", {
  # two disjoint indicator columns (no intercept): X'WX is diagonal for ANY
  # weights, so WS = 1 and X = 1, while MP (weights after scaling) is not 1
  set.seed(10)
  g <- rep(c("a", "b"), c(12, 28))
  y <- c(rpois(12, 2), rpois(28, 9))
  mod <- glm(y ~ 0 + g, family = poisson())
  expect_gt(diff(range(mod$weights)), 1)
  expect_equal(condition_number(mod, method = "WS")$condition_index, 1, tolerance = 1e-8)
  expect_equal(condition_number(mod, method = "X")$condition_index, 1, tolerance = 1e-8)
  expect_gt(condition_number(mod, method = "MP")$condition_index, 1.5)
})

test_that("a perfectly collinear design gives Inf (with a warning) in condition_number and an error in bkw_diagnostics", {
  d <- simulated_data()
  d$x5 <- d$x1 + d$x3
  mod <- suppressWarnings(glm(y ~ x1 + x3 + x5, family = Gamma(link = "inverse"), data = d))
  expect_true(anyNA(coef(mod))) # x5 is aliased
  expect_warning(res <- condition_number(mod, method = "X"), "rank-deficient")
  expect_true(is.infinite(res$condition_number) || res$condition_index > 1e7)
  expect_error(bkw_diagnostics(mod, method = "X"), "rank-deficient")
})

test_that("method = \"RAW\" reproduces the published eigenvalues of the mine fracture data", {
  # Marx (1992); Kurtoglu and Ozkale (2017): eigenvalues of X'WX
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(link = "log"), data = Mine)
  res <- condition_number(mod, method = "RAW")
  expect_identical(res$nc_label, "NC_RAW")
  expect_equal(
    unname(res$eigenvalues),
    c(4164686, 338966.8, 30607.81, 3824.866, 0.9504078),
    tolerance = 1e-4
  )
  expect_equal(res$condition_index, 2093.3225, tolerance = 1e-3)
  k <- length(coef(mod)) / sum(coef(mod)^2)
  expect_equal(k, 0.3871427, tolerance = 1e-4)
})

# ---- Table 1 of the manuscript ------------------------------------------------
table1 <- list(
  lee_yes     = c(X = 190.7765, RAW = 10808.1330, MP = 329.9542, WS = 315.1366),
  lee_no      = c(X = 11.7696,  RAW = 1347.1513,  MP = 33.0269,  WS = 32.3722),
  nitrogen_y  = c(X = 67.7053,  RAW = 8790.6504,  MP = 75.9081,  WS = 77.9772),
  nitrogen_n  = c(X = 14.6751,  RAW = 119.2248,   MP = 13.7742,  WS = 14.3986),
  mine_yes    = c(X = 19.8950,  RAW = 2093.3158,  MP = 26.4538,  WS = 28.5822),
  mine_no     = c(X = 6.4023,   RAW = 31.3666,    MP = 6.9798,   WS = 7.4399),
  plastic_yes = c(X = 595.7429, RAW = 73561.6243, MP = 899.0447, WS = 876.8393),
  plastic_no  = c(X = 342.7521, RAW = 17495.3548, MP = 325.5947, WS = 323.4382)
)

cn_all <- function(f, family, data) {
  mod <- glm(f, family = family, data = data)
  sapply(c("X", "RAW", "MP", "WS"), function(m) condition_number(mod, method = m)$condition_index)
}

test_that("the four methods reproduce Table 1 of the manuscript (packaged datasets)", {
  lee <- LeeCancer
  lee$x6100 <- lee$x6 * 100
  chk <- function(obs, key) expect_equal(round(obs, 4), table1[[key]], tolerance = 1e-6)

  chk(cn_all(y ~ x1 + x4 + x6100, binomial(), lee), "lee_yes")
  chk(cn_all(y ~ x1 + x4 + x6100 - 1, binomial(), lee), "lee_no")
  chk(cn_all(y ~ x1 + x2 + x3 + x4, Gamma("inverse"), Nitrogen), "nitrogen_y")
  chk(cn_all(y ~ x1 + x2 + x3 + x4 - 1, Gamma("inverse"), Nitrogen), "nitrogen_n")
  chk(cn_all(y ~ x1 + x2 + x3 + x4, poisson(), Mine), "mine_yes")
  chk(cn_all(y ~ x1 + x2 + x3 + x4 - 1, poisson(), Mine), "mine_no")
  chk(cn_all(Y ~ X1 + X2 + X3 + X4, poisson(), Plastic), "plastic_yes")
  chk(cn_all(Y ~ X1 + X2 + X3 + X4 - 1, poisson(), Plastic), "plastic_no")
})

test_that("Bodyfat (Gamma, log): weights are all 1, so X = MP = WS", {
  skip_if_not_installed("TH.data")
  data("bodyfat", package = "TH.data", envir = environment())
  mod <- glm(DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
               anthro3a + anthro3b + anthro3c + anthro4,
             family = Gamma(link = "log"), data = bodyfat)
  expect_equal(range(mod$weights), c(1, 1), tolerance = 1e-8)
  cis <- sapply(c("X", "MP", "WS"), function(m) condition_number(mod, method = m)$condition_index)
  expect_equal(round(unname(cis), 4), rep(352.0413, 3), tolerance = 1e-6)
  expect_equal(round(condition_number(mod, method = "RAW")$condition_index, 4), 3410.6317, tolerance = 1e-6)
})

test_that("centering would make the Gamma-inverse Nitrogen design singular (why the package does not center)", {
  mod <- glm(y ~ x1 + x2 + x3 + x4 - 1, family = Gamma(link = "inverse"), data = Nitrogen)
  # sqrt(w) * eta == 1 for every observation
  expect_equal(unname(sqrt(mod$weights) * mod$linear.predictors), rep(1, nrow(Nitrogen)), tolerance = 1e-6)
  Xw <- model.matrix(mod) * sqrt(mod$weights)
  Xc <- scale(Xw, center = TRUE, scale = FALSE)
  expect_lt(min(svd(Xc)$d) / max(svd(Xc)$d), 1e-6)
  # whereas the uncentered WS diagnostic is finite and well defined
  expect_true(is.finite(condition_number(mod)$condition_index))
})

test_that("bkw_diagnostics: proportions sum to one per coefficient and flag the collinear pair", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  for (m in c("WS", "MP", "X", "RAW")) {
    res <- bkw_diagnostics(mod, method = m)
    expect_s3_class(res, "multicollglm_bkw")
    expect_equal(unname(rowSums(res$proportions)), rep(1, 5), tolerance = 1e-8)
    expect_equal(max(res$condition_index), condition_number(mod, method = m)$condition_index, tolerance = 1e-8)
    expect_identical(res$method, m)
  }
  res <- bkw_diagnostics(mod)
  expect_true(length(res$flagged) >= 1)
  expect_true(all(c("x1", "x2") %in% res$flagged[[1]]$variables))
})

test_that("bkw_diagnostics(method = 'X') equals the textbook computation on scale(X, center = FALSE)", {
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(), data = Mine)
  X <- model.matrix(mod)
  s <- svd(scale(X, center = FALSE, scale = sqrt(colSums(X^2))))
  ci <- max(s$d) / s$d
  phi <- sweep(s$v^2, 2L, s$d^2, "/")
  pi_ <- phi / rowSums(phi)
  res <- bkw_diagnostics(mod, method = "X")
  expect_equal(unname(res$condition_index), ci, tolerance = 1e-8)
  expect_equal(unname(res$proportions), unname(pi_), tolerance = 1e-8)
})

test_that("rvif_diagnostics agrees with rvif::rvifs on the method's matrix", {
  skip_if_not_installed("rvif")
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  X <- model.matrix(mod)
  Xw <- X * sqrt(mod$weights)

  res <- rvif_diagnostics(mod)
  expect_s3_class(res, "multicollglm_rvif")
  expect_identical(res$method, "WS")
  expect_equal(rownames(res$table), colnames(X))
  ref <- rvif::rvifs(Xw, ul = TRUE, intercept = TRUE)
  expect_equal(unname(as.matrix(res$table)), unname(as.matrix(ref)), tolerance = 1e-8)

  # rvifs() re-normalizes to unit length, so WS, MP and RAW coincide; X differs
  r_mp  <- rvif_diagnostics(mod, method = "MP")$table
  r_raw <- rvif_diagnostics(mod, method = "RAW")$table
  r_x   <- rvif_diagnostics(mod, method = "X")$table
  expect_equal(unname(as.matrix(r_raw)), unname(as.matrix(res$table)), tolerance = 1e-8)
  expect_equal(unname(as.matrix(r_mp)), unname(as.matrix(res$table)), tolerance = 1e-4)
  expect_false(isTRUE(all.equal(unname(as.matrix(r_x)), unname(as.matrix(res$table)), tolerance = 1e-6)))
})

test_that("print methods run and mention the method", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  expect_output(print(condition_number(mod, method = "X")), "X_ul'X_ul")
  expect_output(print(condition_number(mod)), "NC_WS")
  expect_output(print(bkw_diagnostics(mod, method = "MP")), "method = \"MP\"")
})

test_that("VisCollin::colldiag(add.intercept = TRUE) equals method = 'X'", {
  skip_if_not_installed("VisCollin")
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(), data = Mine)
  cd <- VisCollin::colldiag(mod, add.intercept = TRUE)
  res <- bkw_diagnostics(mod, method = "X")
  expect_equal(as.numeric(cd$condindx), unname(res$condition_index), tolerance = 1e-8)
  # colldiag stores components in rows and variables in columns
  expect_equal(unname(t(unclass(cd$pi))), unname(res$proportions), tolerance = 1e-8)
})

test_that("glmnet methods refit on the active set and take 'method'", {
  skip_if_not_installed("glmnet")
  d <- simulated_data()
  x <- as.matrix(d[, c("x1", "x2", "x3", "x4")])
  y <- d$y
  fit <- glmnet::glmnet(x, y, family = Gamma(link = "inverse"), lambda = c(0.05, 0.01, 0.001))
  s <- fit$lambda[2]
  cf <- coef(fit, s = s)
  active <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
  manual <- glm(reformulate(active, response = "y"), data = d, family = Gamma(link = "inverse"))

  for (m in c("WS", "MP", "X", "RAW")) {
    res <- condition_number(fit, x = x, y = y, family = Gamma(link = "inverse"), s = s, method = m)
    expect_s3_class(res, "multicollglm_cn")
    expect_equal(res$eigenvalues, condition_number(manual, method = m)$eigenvalues, tolerance = 1e-5)
    expect_equal(sort(res$active_vars), sort(active))
  }
  expect_error(condition_number(fit, x = x, y = y, family = Gamma(link = "inverse")), "lambda")
  expect_error(condition_number(fit, x = x, y = y, s = s), "family")
  expect_error(condition_number(fit, x = x, y = y, family = Gamma(link = "inverse"), s = s, method = "OZ"), "removed")

  bk <- bkw_diagnostics(fit, x = x, y = y, family = Gamma(link = "inverse"), s = s, method = "X")
  expect_equal(unname(bk$condition_index), unname(bkw_diagnostics(manual, method = "X")$condition_index), tolerance = 1e-5)

  cvfit <- glmnet::cv.glmnet(x, y, family = Gamma(link = "inverse"), nfolds = 5)
  expect_s3_class(condition_number(cvfit, x = x, y = y, family = Gamma(link = "inverse")), "multicollglm_cn")
})
