#' ML-collinearity diagnostic (Lesaffre and Marx, 1993)
#'
#' Lesaffre and Marx (1993) distinguish two, structurally different causes
#' of an ill-conditioned (estimated) information matrix in a GLM: ordinary
#' collinearity among the explanatory variables themselves, and what they
#' call *ML-collinearity* -- ill-conditioning that is not attributable to
#' the explanatory variables, but instead arises from the combination of
#' the response, the fitted coefficients and the link function. To tell
#' the two apart they compare `sqrt(CN)` of the original (unweighted)
#' design matrix `X` against `sqrt(CN)` of the IRLS-weighted information
#' matrix under Mackinnon and Puterman's (1989) definition (`method =
#' "MP"`; both quantities standardize the columns of `X`, including the
#' intercept, to unit length first).
#'
#' `ml_collinearity()` computes `sqrt(CN)` of `X` with
#' [multiColl::CNs()], following Lesaffre and Marx (1993, Sec. 4.2)'s own
#' recommendation. `CNs(X)` (with `X` including the intercept column)
#' returns a list of two elements, the condition number without and with
#' the intercept (in that order; when `X` has only one regressor besides
#' the intercept, it instead returns a bare numeric scalar, the
#' with-intercept value, since the without-intercept quantity is undefined
#' for a single column). Since Lesaffre and Marx's `kappa_X` is explicitly
#' computed on `X` **with** the constant vector standardized to unit
#' length along with the rest, this uses the with-intercept value, which
#' reproduces Lesaffre and Marx's (1993) own published `kappa_X = 190.78`
#' exactly for their Lee cancer-remission example, while the
#' without-intercept value does not. `sqrt(CN)` of the MacKinnon-Puterman
#' information matrix is computed with `condition_number(model, method =
#' "MP")`, and the two are combined into their ratio, `ratio_wx =
#' sqrt(CN)_W / sqrt(CN)_X`. Following
#' Lesaffre and Marx (1993, Sec. 4.2), a model is flagged with
#' ML-collinearity only when **both** `ratio_wx > ratio_threshold`
#' (default 5) **and** `sqrt(CN)_W > kappa_threshold` (default 30): the
#' ratio alone is not sufficient, since a high ratio with a small
#' `sqrt(CN)_W` does not indicate an ill-conditioned information matrix to
#' begin with. Separately, `sqrt(CN)_X > kappa_threshold` flags ordinary
#' collinearity among the explanatory variables. If both flags are set and
#' the ratio is close to 1, Lesaffre and Marx (1993) note that both types
#' of collinearity are simultaneously present.
#'
#' @param model A `glm`, `glmnet` or `cv.glmnet` object.
#' @param ratio_threshold Threshold for `ratio_wx` above which the ratio is
#'   considered high (default 5, as suggested by Lesaffre and Marx, 1993).
#' @param kappa_threshold Threshold for `sqrt(CN)_W`/`sqrt(CN)_X` above
#'   which the corresponding condition number is considered to indicate
#'   ill-conditioning (default 30, as suggested by Lesaffre and Marx, 1993,
#'   following Mackinnon and Puterman, 1989 and the classical Belsley, Kuh
#'   and Welsch, 1980 rule of thumb on the `sqrt(CN)` scale).
#' @param ... Additional arguments passed to the methods (see `x`, `y`,
#'   `family`, `s`, `weights`, `offset` for `glmnet`/`cv.glmnet`).
#' @param x Predictor matrix (the same one used to fit the `glmnet`
#'   model), with `colnames`.
#' @param y Response variable used to fit the `glmnet` model.
#' @param family GLM family (a `family` object, e.g.
#'   `Gamma(link = "inverse")`) used to refit a `glm` on the active set of
#'   variables at the given lambda.
#' @param s Lambda value (numeric) for a `glmnet` object. For `cv.glmnet`
#'   it can also be `"lambda.min"` or `"lambda.1se"`.
#' @param weights Optional prior weights for the refit.
#' @param offset Optional offset for the refit.
#'
#' @return An object of class `multicollglm_mlcoll` with `kappa_x`
#'   (`sqrt(CN)` of the original design matrix `X`, via
#'   `multiColl::CNs(X)`, the "with intercept" condition number),
#'   `kappa_w` (`sqrt(CN)` of the
#'   MacKinnon-Puterman information matrix, i.e.
#'   `condition_number(model, method = "MP")$condition_index`), `ratio_wx`
#'   (`kappa_w / kappa_x`), the two thresholds used (`ratio_threshold`,
#'   `kappa_threshold`), and the logical flags `collinearity_x` (`kappa_x >
#'   kappa_threshold`) and `ml_collinearity` (`ratio_wx > ratio_threshold`
#'   and `kappa_w > kappa_threshold`).
#'
#' @references
#' Lesaffre, E. and Marx, B.D. (1993). Collinearity in generalized linear
#' regression. Communications in Statistics - Theory and Methods, 22(7),
#' 1933-1952. \doi{10.1080/03610929308831126}
#'
#' Mackinnon, M.J. and Puterman, M.L. (1989). Collinearity in generalized
#' linear models. Communications in Statistics - Theory and Methods, 18(9),
#' 3463-3472. \doi{10.1080/03610928908830102}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
#' mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
#' y <- rgamma(n, shape = 5, rate = 5 / mu)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
#' \dontrun{
#' ml_collinearity(mod) # requires the 'multiColl' package
#' }
#'
#' @export
ml_collinearity <- function(model, ...) {
  UseMethod("ml_collinearity")
}

#' @export
ml_collinearity.glm <- function(model, ratio_threshold = 5, kappa_threshold = 30, ...) {
  if (!requireNamespace("multiColl", quietly = TRUE)) {
    stop(
      "ml_collinearity() requires the 'multiColl' package, which supplies CNs(), used to ",
      "compute the classical condition number of the original design matrix 'X'. Please ",
      "install it.",
      call. = FALSE
    )
  }

  X <- stats::model.matrix(model)
  if (ncol(X) < 2L) {
    stop("ml_collinearity() needs at least two columns in the design matrix.", call. = FALSE)
  }

  # CNs(X) (X including the intercept column) returns a list of two named
  # elements -- "Condition Number without intercept" and "Condition Number
  # with intercept" (in that order; the names are descriptive strings, not
  # "CN1"/"CN2") -- when X has 3+ columns, but a bare numeric scalar when X
  # has only 2 columns (intercept + a single regressor), since dropping the
  # intercept then leaves a single column for which collinearity is
  # undefined. Lesaffre and Marx's (1993) kappa_X standardizes X
  # *including* the constant vector, i.e. the "with intercept" value -- the
  # second list element, or the lone scalar in the 2-column case. Using the
  # second element reproduces their own published kappa_X = 190.78 for the
  # Lee cancer-remission example; the "without intercept" value does not.
  cns_x <- multiColl::CNs(X)
  kappa_x <- if (is.list(cns_x)) as.numeric(cns_x[[2]]) else as.numeric(cns_x[1])

  cn_w <- condition_number(model, method = "MP")
  kappa_w <- cn_w$condition_index

  ratio_wx <- kappa_w / kappa_x

  structure(
    list(
      kappa_x = kappa_x,
      kappa_w = kappa_w,
      ratio_wx = ratio_wx,
      ratio_threshold = ratio_threshold,
      kappa_threshold = kappa_threshold,
      collinearity_x = kappa_x > kappa_threshold,
      ml_collinearity = (ratio_wx > ratio_threshold) && (kappa_w > kappa_threshold),
      var_names = colnames(X)
    ),
    class = "multicollglm_mlcoll"
  )
}

#' @rdname ml_collinearity
#' @export
ml_collinearity.glmnet <- function(model, x, y, family, s, weights = NULL, offset = NULL,
                                    ratio_threshold = 5, kappa_threshold = 30, ...) {
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- ml_collinearity.glm(refit$model, ratio_threshold = ratio_threshold, kappa_threshold = kappa_threshold)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname ml_collinearity
#' @export
ml_collinearity.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                       ratio_threshold = 5, kappa_threshold = 30, ...) {
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- ml_collinearity.glm(refit$model, ratio_threshold = ratio_threshold, kappa_threshold = kappa_threshold)
  out$lambda <- if (is.character(s)) model[[s]] else s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @export
print.multicollglm_mlcoll <- function(x, digits = 4, ...) {
  cat("ML-collinearity diagnostic (Lesaffre and Marx, 1993)\n")
  if (!is.null(x$lambda)) {
    cat(sprintf(
      "lambda = %s | active variables: %s\n",
      format(x$lambda), paste(x$active_vars, collapse = ", ")
    ))
    if (length(x$dropped_vars) > 0) {
      cat(sprintf("dropped variables (coef = 0): %s\n", paste(x$dropped_vars, collapse = ", ")))
    }
  }

  cat(sprintf("\nsqrt(CN) of X (original design matrix, multiColl::CNs(X), with intercept): %.*f\n", digits, x$kappa_x))
  cat(sprintf("sqrt(CN) of W (MacKinnon-Puterman information matrix, method = \"MP\"): %.*f\n", digits, x$kappa_w))
  cat(sprintf("ratio (sqrt(CN)_W / sqrt(CN)_X): %.*f\n", digits, x$ratio_wx))

  cat(sprintf(
    "\nCollinearity in X (sqrt(CN)_X > %g): %s\n",
    x$kappa_threshold, if (x$collinearity_x) "YES" else "no"
  ))
  cat(sprintf(
    "ML-collinearity (ratio_wx > %g and sqrt(CN)_W > %g): %s\n",
    x$ratio_threshold, x$kappa_threshold, if (x$ml_collinearity) "YES" else "no"
  ))

  if (x$collinearity_x && x$ml_collinearity) {
    cat("\n>> Both collinearity in X and ML-collinearity appear to be present.\n")
  }
  invisible(x)
}
