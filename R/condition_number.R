#' Condition number of a model (collinearity diagnostic)
#'
#' Computes the condition number of a generalized linear model, following
#' one of four definitions chosen with `method`. By default
#' (`method = "WS"`) it follows Weissfeld and Sereika (1991): the model is
#' fitted on the original design matrix and the information matrix
#' \eqn{X'\hat{W}X} at convergence is rescaled so that its diagonal is all
#' ones, which is equivalent to scaling the IRLS-weighted design matrix
#' \eqn{\hat{W}^{1/2}X} to unit length by column, *without* centering. This
#' generalizes to any GLM family and link function the classical condition
#' number of linear regression. The other definitions are the one of
#' MacKinnon and Puterman (1989), the one of the original (unweighted)
#' design matrix \eqn{X'X}, and the one of the untransformed information
#' matrix \eqn{X'\hat{W}X}.
#'
#' @details
#' ## Methods
#'
#' Let \eqn{X} be the design matrix (including the intercept column, if
#' any), \eqn{\hat{W}} the diagonal matrix of IRLS weights at convergence,
#' and let the subscript \eqn{ul} denote scaling every column to unit
#' Euclidean length. None of the four methods centers the data, so the
#' intercept always takes part in the diagnostic and non-essential
#' collinearity (between the intercept and the explanatory variables) can
#' be detected.
#'
#' | `method` | Matrix whose eigenvalues are used | Reference |
#' |:---|:---|:---|
#' | `"WS"` (default) | \eqn{(X'\hat{W}X)_{ul} = S\,X'\hat{W}X\,S}, with \eqn{S = \mathrm{diag}(X'\hat{W}X)^{-1/2}} | Weissfeld and Sereika (1991) |
#' | `"MP"` | \eqn{X_{ul}'\hat{W}X_{ul}}, with \eqn{X} scaled to unit length *before* refitting the model | MacKinnon and Puterman (1989) |
#' | `"X"` | \eqn{X_{ul}'X_{ul}}: the original design matrix, no weights | Belsley, Kuh and Welsch (1980); \eqn{\kappa_X} of Lesaffre and Marx (1993) |
#' | `"RAW"` | \eqn{X'\hat{W}X}: no transformation at all | -- |
#'
#' * `"WS"` fits the GLM on the raw `X` and rescales afterwards.
#'   Under orthogonality in the GLM metric (\eqn{X'\hat{W}X} diagonal) it
#'   attains its minimum value, 1, for any weights, exactly as the
#'   unit-length condition number does in the linear model. This is the
#'   recommended definition.
#' * `"MP"` scales `X` to unit length first, refits the GLM on that matrix
#'   (this is done internally) and reads the information matrix off the
#'   refit. Because the weights enter after the scaling, the result is 1 under
#'   orthogonality only if all the weights are equal.
#' * `"X"` ignores the model's weights, so it measures the collinearity of
#'   the explanatory variables themselves, not of the information matrix.
#'   It is what a linear-model diagnostic computes on the same predictors,
#'   e.g. `VisCollin::colldiag(mod, add.intercept = TRUE)`.
#' * `"RAW"` is the baseline without any transformation. It depends on the
#'   units of the explanatory variables and is usually far larger than the
#'   others, even when the underlying collinearity is the same.
#'
#' The condition number is returned on two scales, as in the literature:
#' `condition_number`, the ratio \eqn{\lambda_{max}/\lambda_{min}} of the
#' eigenvalues, and `condition_index`, its square root, which is the
#' classical condition number computed from the singular values of the
#' matrix. Thresholds such as 10 (weak), 30 (moderate to strong) or 100
#' (Belsley, Kuh and Welsch, 1980) are stated on the square-root scale.
#'
#' ## Removed in version 0.4.0
#'
#' `method = "OZ"` (Ozkale, 2021) and `method = "MS"` (Marx and Smith, 1990),
#' as well as the arguments `center` and `scale`, no longer exist. All of
#' them centered the data, which removes the intercept from the analysis
#' and so rules out the detection of non-essential collinearity.
#'
#' @param model A `glm`, `glmnet` or `cv.glmnet` object.
#' @param method One of `"WS"` (default), `"MP"`, `"X"` or `"RAW"`. See
#'   Details.
#' @param ... Additional arguments passed to the methods (see `x`, `y`,
#'   `family`, `s`, `weights`, `offset` for `glmnet`/`cv.glmnet`).
#'
#' @return An object of class `multicollglm_cn` with the eigenvalues of the
#'   matrix defined by `method` (`eigenvalues`), the condition number as
#'   the ratio between the largest and the smallest eigenvalue
#'   (`condition_number`), its square root on the classical singular-value
#'   scale (`condition_index`), the variable names (`var_names`), the
#'   method used (`method`) and its label (`nc_label`, e.g. `"NC_WS"`). If
#'   the matrix is exactly rank-deficient (perfect collinearity), the
#'   condition number is `Inf` and a warning is issued.
#'
#' @references
#' Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). *Regression
#' Diagnostics: Identifying Influential Data and Sources of Collinearity*.
#' Wiley.
#'
#' Lesaffre, E. and Marx, B.D. (1993). Collinearity in generalized linear
#' regression. *Communications in Statistics - Theory and Methods*, 22(7),
#' 1933-1952. \doi{10.1080/03610929308831126}
#'
#' MacKinnon, M.J. and Puterman, M.L. (1989). Collinearity in generalized
#' linear models. *Communications in Statistics - Theory and Methods*,
#' 18(9), 3463-3472. \doi{10.1080/03610928908830102}
#'
#' Weissfeld, L.A. and Sereika, S.M. (1991). A multicollinearity diagnostic
#' for generalized linear models. *Communications in Statistics - Theory and
#' Methods*, 20(4), 1183-1198. \doi{10.1080/03610929108830558}
#'
#' @seealso [bkw_diagnostics()], [rvif_diagnostics()], [ml_collinearity()].
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
#' mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
#' y <- rgamma(n, shape = 5, rate = 5 / mu)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
#'
#' condition_number(mod) # Weissfeld and Sereika (1991), the default
#'
#' # The four definitions side by side (sqrt scale)
#' sapply(c("WS", "MP", "X", "RAW"),
#'        function(m) condition_number(mod, method = m)$condition_index)
#'
#' @export
condition_number <- function(model, ...) {
  UseMethod("condition_number")
}

#' @export
condition_number.glm <- function(model, method = c("WS", "MP", "X", "RAW"), ...) {
  .reject_removed_args(...)
  method <- .check_method(method)

  core <- .build_design(model, method = method)
  ev_raw <- core$eigenvalues

  # crossprod() of a real matrix is symmetric positive semi-definite, so its
  # true eigenvalues are never negative. An eigenvalue that is <= 0, or below
  # the numerical-rank tolerance max(ev) * p * .Machine$double.eps (the usual
  # rule of qr()/rankMatrix), is floating-point noise around an exactly zero
  # eigenvalue, i.e. an exactly singular matrix (perfect collinearity).
  # Report the mathematically correct Inf instead of a nonsensical negative
  # or platform-dependent CN that is really just noise.
  tol <- max(ev_raw) * length(ev_raw) * .Machine$double.eps
  singular <- min(ev_raw) <= tol
  if (singular) {
    warning(
      "The (weighted, transformed) design matrix is numerically rank-deficient (an ",
      "eigenvalue is zero up to floating-point error), i.e. there is perfect collinearity; ",
      "CN is mathematically infinite and is reported as Inf.",
      call. = FALSE
    )
  }
  ev <- pmax(ev_raw, 0)
  cn <- if (singular) Inf else max(ev) / min(ev)

  structure(
    list(
      eigenvalues = ev,
      condition_number = cn,
      condition_index = sqrt(cn),
      var_names = core$var_names,
      method = method,
      nc_label = paste0("NC_", method)
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
                                     method = c("WS", "MP", "X", "RAW"), ...) {
  .reject_removed_args(...)
  method <- .check_method(method)
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- condition_number.glm(refit$model, method = method)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname condition_number
#' @export
condition_number.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                        method = c("WS", "MP", "X", "RAW"), ...) {
  .reject_removed_args(...)
  method <- .check_method(method)
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- condition_number.glm(refit$model, method = method)
  out$lambda <- if (is.character(s)) model[[s]] else s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

# Label of the matrix whose eigenvalues each method uses (shared by the
# print methods).
.matrix_label <- function(method) {
  switch(
    method,
    "WS"  = "(X'WX)_ul",
    "MP"  = "X_ul'WX_ul",
    "X"   = "X_ul'X_ul",
    "RAW" = "X'WX"
  )
}

#' @export
print.multicollglm_cn <- function(x, digits = 4, ...) {
  cat(sprintf("Eigenvalues of %s (decreasing order):\n", .matrix_label(x$method)))
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

  cat(sprintf("\nCN (%s, eigenvalue scale): %.*f\n", x$nc_label, digits, x$condition_number))
  cat(sprintf("sqrt(CN) (%s, classical, singular-value scale): %.*f\n", x$nc_label, digits, x$condition_index))
  invisible(x)
}
