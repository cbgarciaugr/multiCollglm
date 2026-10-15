#' Redefined Variance Inflation Factor (RVIF) diagnostics
#'
#' Computes the Redefined Variance Inflation Factor (RVIF) of Salmeron,
#' Garcia and Garcia (2025), via [rvif::rvifs()], on the same matrix used
#' by [condition_number()] and [bkw_diagnostics()], chosen with `method`.
#' By default this is the IRLS-weighted design matrix scaled to unit length
#' by column, without centering (Weissfeld and Sereika, 1991). Unlike the
#' classical VIF, the RVIF is able to detect both essential collinearity
#' (near-linear relationships among the regressors) and non-essential
#' collinearity (near-linear relationships between the intercept and one or
#' more regressors), and also decomposes, for each variable, the percentage
#' of near collinearity it is responsible for.
#'
#' [rvif::rvifs()] is always called with `ul = TRUE`, so it re-scales the
#' matrix to Euclidean unit length right before computing the RVIF. As a
#' consequence `"WS"`, `"MP"` and `"RAW"`, which only differ in how the
#' weighted matrix is scaled, give the same RVIF; what changes the result
#' is whether the weights take part at all, i.e. `"X"` (the original
#' design matrix, ignoring the weights) against the other three.
#'
#' @inheritParams condition_number
#' @param tol Tolerance used by [rvif::rvifs()] to decide whether the
#'   system is computationally singular (default `1e-30`).
#'
#' @return An object of class `multicollglm_rvif` wrapping the data frame
#'   returned by [rvif::rvifs()] (`table`, with columns `RVIF` and `%`,
#'   one row per column of the design matrix, named after the model's
#'   variables) and the `method` used.
#'
#' @references
#' Salmeron, R., Garcia, C.B. and Garcia, J. (2025). A redefined Variance
#' Inflation Factor: overcoming the limitations of the Variance Inflation
#' Factor. Computational Economics, 65, 337-363.
#' \doi{10.1007/s10614-024-10575-8}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
#' mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
#' y <- rgamma(n, shape = 5, rate = 5 / mu)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
#' rvif_diagnostics(mod)
#'
#' @export
rvif_diagnostics <- function(model, ...) {
  UseMethod("rvif_diagnostics")
}

#' @export
rvif_diagnostics.glm <- function(model, method = c("WS", "MP", "X", "RAW"), tol = 1e-30, ...) {
  .reject_removed_args(...)
  method <- .check_method(method)

  core <- .build_design(model, method = method)

  # ul = TRUE: rvifs() re-scales core$M to Euclidean unit length right before
  # computing the RVIF, which is a no-op for the matrices that are already
  # unit length ("WS" is, "X" is) and equalizes the scaling of "MP" and
  # "RAW" (see the Details above).
  tab <- rvif::rvifs(core$M, ul = TRUE, intercept = core$has_intercept, tol = tol)
  if (is.null(tab)) {
    stop(
      "rvif::rvifs() could not compute the RVIF for this (weighted, transformed) design matrix ",
      "(see the message above from rvifs() itself, e.g. exact or near-exact multicollinearity).",
      call. = FALSE
    )
  }
  rownames(tab) <- core$var_names

  structure(
    list(table = tab, var_names = core$var_names, method = method),
    class = "multicollglm_rvif"
  )
}

#' @rdname rvif_diagnostics
#' @export
rvif_diagnostics.glmnet <- function(model, x, y, family, s, weights = NULL, offset = NULL,
                                     method = c("WS", "MP", "X", "RAW"), tol = 1e-30, ...) {
  .reject_removed_args(...)
  method <- .check_method(method)
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- rvif_diagnostics.glm(refit$model, method = method, tol = tol)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname rvif_diagnostics
#' @export
rvif_diagnostics.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                        method = c("WS", "MP", "X", "RAW"), tol = 1e-30, ...) {
  .reject_removed_args(...)
  method <- .check_method(method)
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- rvif_diagnostics.glm(refit$model, method = method, tol = tol)
  out$lambda <- if (is.character(s)) model[[s]] else s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @export
print.multicollglm_rvif <- function(x, digits = 3, ...) {
  cat(sprintf("Redefined Variance Inflation Factor (RVIF) on %s [method = \"%s\"]\n",
              .matrix_label(x$method), x$method))
  if (!is.null(x$lambda)) {
    cat(sprintf(
      "lambda = %s | active variables: %s\n",
      format(x$lambda), paste(x$active_vars, collapse = ", ")
    ))
    if (length(x$dropped_vars) > 0) {
      cat(sprintf("dropped variables (coef = 0): %s\n", paste(x$dropped_vars, collapse = ", ")))
    }
  }
  cat("\n")
  print(round(x$table, digits))
  invisible(x)
}
