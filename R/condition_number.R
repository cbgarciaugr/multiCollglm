#' Condition number of a model (collinearity diagnostic)
#'
#' Computes the condition number over the design matrix weighted by the
#' IRLS weights at convergence and, by default, scaled to unit length by
#' column (without centering). This generalizes, to any GLM family and
#' link function, the classical condition number of linear regression.
#'
#' @param model A `glm`, `glmnet` or `cv.glmnet` object.
#' @param center Either `FALSE` (default; no centering, as originally
#'   proposed and required to detect non-essential collinearity involving
#'   the intercept) or any value accepted by [scale()]'s own `center`
#'   argument: `TRUE` to center each column on its mean, or a numeric
#'   vector of per-column values to subtract. For a canonical-link GLM
#'   with an intercept, centering the IRLS-weighted design matrix is
#'   mathematically guaranteed to drop its rank by exactly one (since
#'   `sqrt(w) * eta` is then constant), which [bkw_diagnostics()] and
#'   [rvif_diagnostics()] report as an error rather than silent `NaN`s;
#'   this is the classical reason Belsley, Kuh and Welsch advise against
#'   centering for collinearity diagnostics.
#' @param scale Either `"unit"` (default; each column is scaled to
#'   Euclidean/L2 unit length, as originally proposed and required for the
#'   classical condition-index thresholds and the RVIF formula) or any
#'   value accepted by [scale()]'s own `scale` argument: `TRUE` for
#'   root-mean-square scaling, `FALSE` for no scaling, or a numeric vector
#'   of per-column divisors.
#' @param ... Additional arguments passed to the methods (see `x`, `y`,
#'   `family`, `s`, `weights`, `offset` for `glmnet`/`cv.glmnet`).
#'
#' @return An object of class `multicollglm_cn` with, among others, the
#'   eigenvalues of \eqn{X_{unit}'X_{unit}}, the condition number on the
#'   eigenvalue scale (`condition_number`), and the classical condition
#'   index on the singular-value scale (`condition_index`, i.e.
#'   `sqrt(condition_number)`).
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
#' mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
#' y <- rgamma(n, shape = 5, rate = 5 / mu)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
#' condition_number(mod)
#'
#' @export
condition_number <- function(model, ...) {
  UseMethod("condition_number")
}

#' @export
condition_number.glm <- function(model, center = FALSE, scale = "unit", ...) {
  core <- .build_design(model, center = center, scale = scale)
  ev <- core$eigenvalues
  cn <- max(ev) / min(ev)

  structure(
    list(
      eigenvalues = ev,
      condition_number = cn,
      condition_index = sqrt(cn),
      var_names = core$var_names
    ),
    class = "multicollglm_cn"
  )
}

#' @rdname condition_number
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
#' @export
condition_number.glmnet <- function(model, x, y, family, s, weights = NULL, offset = NULL,
                                     center = FALSE, scale = "unit", ...) {
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- condition_number.glm(refit$model, center = center, scale = scale)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname condition_number
#' @export
condition_number.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                        center = FALSE, scale = "unit", ...) {
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- condition_number.glm(refit$model, center = center, scale = scale)
  out$lambda <- if (is.character(s)) model[[s]] else s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @export
print.multicollglm_cn <- function(x, digits = 4, ...) {
  cat("Eigenvalues of X_unit'X_unit (decreasing order):\n")
  print(round(x$eigenvalues, digits))

  if (!is.null(x$lambda)) {
    cat(sprintf(
      "\nlambda = %s | active variables: %s\n",
      format(x$lambda), paste(x$active_vars, collapse = ", ")
    ))
    if (length(x$dropped_vars) > 0) {
      cat(sprintf("dropped variables (coef = 0): %s\n", paste(x$dropped_vars, collapse = ", ")))
    }
  }

  cat(sprintf("\nCondition number (eigenvalue scale): %.*f\n", digits, x$condition_number))
  cat(sprintf("Condition index (classical, sqrt): %.*f\n", digits, x$condition_index))
  invisible(x)
}
