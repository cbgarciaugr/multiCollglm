test_that("ml_collinearity reproduces Lesaffre and Marx (1993) for the Lee cancer data", {
  lee <- LeeCancer
  lee$x6100 <- lee$x6 * 100
  mod <- glm(y ~ x1 + x4 + x6100, family = binomial(), data = lee)
  res <- ml_collinearity(mod)

  expect_s3_class(res, "multicollglm_mlcoll")
  expect_equal(res$kappa_x, 190.78, tolerance = 1e-4)
  expect_equal(res$kappa_w, 329.95, tolerance = 1e-4)
  expect_equal(res$ratio_wx, 329.9542 / 190.7765, tolerance = 1e-5)
  expect_true(res$collinearity_x)
})

test_that("ml_collinearity uses the native 'X' and 'MP' condition numbers", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  res <- ml_collinearity(mod)
  expect_equal(res$kappa_x, condition_number(mod, method = "X")$condition_index)
  expect_equal(res$kappa_w, condition_number(mod, method = "MP")$condition_index)
  expect_equal(res$ratio_wx, res$kappa_w / res$kappa_x)
})

test_that("ML-collinearity needs BOTH a high ratio and a high kappa_w", {
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(), data = Mine)
  none <- ml_collinearity(mod, ratio_threshold = 1000, kappa_threshold = 1000)
  expect_false(none$ml_collinearity)
  expect_false(none$collinearity_x)
  # ratio passes but kappa_w does not
  ratio_only <- ml_collinearity(mod, ratio_threshold = 0, kappa_threshold = 1e6)
  expect_false(ratio_only$ml_collinearity)
  both <- ml_collinearity(mod, ratio_threshold = 0, kappa_threshold = 0)
  expect_true(both$ml_collinearity)
})

test_that("ml_collinearity does not require the 'multiColl' package", {
  d <- simulated_data()
  mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = d)
  expect_no_error(ml_collinearity(mod))
  expect_output(print(ml_collinearity(mod)), "ML-collinearity")
})

test_that("ml_collinearity.glmnet / cv.glmnet refit on the active set", {
  skip_if_not_installed("glmnet")
  d <- simulated_data()
  x <- as.matrix(d[, c("x1", "x2", "x3", "x4")])
  y <- d$y
  fit <- glmnet::glmnet(x, y, family = Gamma(link = "inverse"), lambda = c(0.05, 0.01, 0.001))
  s <- fit$lambda[2]
  res <- ml_collinearity(fit, x = x, y = y, family = Gamma(link = "inverse"), s = s)
  expect_s3_class(res, "multicollglm_mlcoll")
  expect_error(ml_collinearity(fit, x = x, y = y, family = Gamma(link = "inverse")), "lambda")
  cvfit <- glmnet::cv.glmnet(x, y, family = Gamma(link = "inverse"), nfolds = 5)
  expect_s3_class(ml_collinearity(cvfit, x = x, y = y, family = Gamma(link = "inverse")), "multicollglm_mlcoll")
})
