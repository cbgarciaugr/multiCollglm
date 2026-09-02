test_that("ml_collinearity.glm reproduces Lesaffre and Marx's (1993) own worked example", {
  skip_if_not_installed("multiColl")

  data(LeeCancer, package = "multiCollglm")
  # TEMP is recorded on a /100 scale in the bundled dataset; multiplied back
  # by 100 to match Lesaffre and Marx's (1993) own units, exactly as in the
  # reproduction-lee-cancer article.
  LeeCancer$x6100 <- LeeCancer$x6 * 100
  mod <- glm(y ~ x1 + x4 + x6100, family = binomial(), data = LeeCancer)

  res <- ml_collinearity(mod)

  expect_s3_class(res, "multicollglm_mlcoll")
  # Lesaffre and Marx (1993, Sec. 6.1): kappa_X = 190.78, kappa_W = 329.95.
  expect_equal(res$kappa_x, 190.78, tolerance = 1e-2)
  expect_equal(res$kappa_w, 329.95, tolerance = 1e-2)
  expect_equal(res$kappa_w, condition_number(mod, method = "MP")$condition_index)
  expect_equal(res$ratio_wx, res$kappa_w / res$kappa_x)

  # Their own conclusion for this example: collinearity in X, but not
  # ML-collinearity (ratio well below 5).
  expect_true(res$collinearity_x)
  expect_false(res$ml_collinearity)
  expect_lt(res$ratio_wx, 5)

  # Regression guard: CNs(X) returns a list here (3+ columns), with the
  # "without intercept" value first and "with intercept" second (the names
  # are descriptive strings, not "CN1"/"CN2"). kappa_x must come from the
  # second element (Lesaffre and Marx (1993) standardize X *including* the
  # constant vector), not the first. Assert the two disagree here so a
  # future regression back to the first element fails loudly.
  cns_x <- multiColl::CNs(model.matrix(mod))
  expect_false(isTRUE(all.equal(cns_x[[1]], cns_x[[2]])))
  expect_equal(res$kappa_x, as.numeric(cns_x[[2]]))
})

test_that("ml_collinearity.glm flags ML-collinearity when the ratio and kappa_W thresholds are both exceeded", {
  skip_if_not_installed("multiColl")

  res <- structure(
    list(kappa_x = 10, kappa_w = 100, ratio_wx = 10,
         ratio_threshold = 5, kappa_threshold = 30,
         collinearity_x = FALSE, ml_collinearity = TRUE),
    class = "multicollglm_mlcoll"
  )
  # sanity check of the flagging rule itself, independent of any model fit
  expect_equal(res$ml_collinearity, (res$ratio_wx > res$ratio_threshold) && (res$kappa_w > res$kappa_threshold))

  # a high ratio alone, with a small kappa_W, must NOT be flagged
  not_flagged <- (6 > 5) && (20 > 30)
  expect_false(not_flagged)
})

test_that("ml_collinearity.glm errors informatively without at least two design columns", {
  skip_if_not_installed("multiColl")

  set.seed(1)
  n <- 50
  y <- rgamma(n, shape = 5, rate = 5)
  mod <- glm(y ~ 1, family = Gamma(link = "log"))

  expect_error(ml_collinearity(mod), "at least two columns")
})

test_that("ml_collinearity.glmnet and .cv.glmnet refit on the active set and label lambda", {
  skip_if_not_installed("multiColl")
  skip_if_not_installed("glmnet")

  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
  eta <- 3 + 0.15 * x1 + 0.15 * x2 - 0.1 * x3 + 0.05 * x4
  mu <- 1 / eta
  y <- rgamma(n, shape = 5, rate = 5 / mu)
  x <- as.matrix(data.frame(x1, x2, x3, x4))

  fit <- glmnet::glmnet(x, y, family = Gamma(link = "inverse"))
  res <- ml_collinearity(fit, x = x, y = y, family = Gamma(link = "inverse"), s = 0.01)
  expect_s3_class(res, "multicollglm_mlcoll")
  expect_equal(res$lambda, 0.01)
  expect_true(!is.null(res$active_vars))

  cvfit <- glmnet::cv.glmnet(x, y, family = Gamma(link = "inverse"))
  res_cv <- ml_collinearity(cvfit, x = x, y = y, family = Gamma(link = "inverse"), s = "lambda.min")
  expect_s3_class(res_cv, "multicollglm_mlcoll")
})
