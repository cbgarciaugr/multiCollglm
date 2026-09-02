#' Redefined Variance Inflation Factor (RVIF) diagnostics
#'
#' Computes the Redefined Variance Inflation Factor (RVIF) of Salmeron,
#' Garcia and Garcia (2025), via [rvif::rvifs()], on the same design
#' matrix used by [condition_number()] and [bkw_diagnostics()]: weighted
#' by the IRLS weights at convergence and, by default, scaled to unit
#' length by column (without centering); see `center`/`scale` below.
#' [rvif::rvifs()] is always called with `ul = TRUE`, so it re-scales
#' that matrix to Euclidean unit length right before computing the RVIF,
#' whatever `center`/`scale` were used to build it (a matrix that is
#' already unit length, the default, is left unchanged by this step).
#' As a consequence, **`scale` has no effect on the RVIF result**: any
#' value other than `"unit"` is silently overridden by that final
#' re-scaling, so a warning is issued whenever `scale != "unit"` is
#' requested (`center` is unaffected and still changes the result, since
#' `rvifs()` does not undo it). Unlike the classical VIF, the RVIF is able to
#' detect both essential collinearity (near-linear relationships among
#' the regressors) and non-essential collinearity (near-linear
#' relationships between the intercept and one or more regressors), and
#' also decomposes, for each variable, the percentage of near
#' collinearity it is responsible for.
#'
#' @inheritParams condition_number
#' @param scale Either `"unit"` (default) or any value accepted by
#'   [scale()]'s own `scale` argument (`TRUE`, `FALSE`, or a numeric
#'   vector of per-column divisors) -- but note that, unlike in
#'   [condition_number()] and [bkw_diagnostics()], **this has no effect
#'   on the RVIF result**: [rvif::rvifs()] is always called with `ul =
#'   TRUE`, which re-scales the design matrix to Euclidean unit length
#'   right before computing the RVIF regardless of `scale`. Requesting
#'   anything other than `"unit"` issues a warning for this reason.
#' @param tol Tolerance used by [rvif::rvifs()] to decide whether the
#'   system is computationally singular (default `1e-30`).
#'
#' @return An object of class `multicollglm_rvif` wrapping the data frame
#'   returned by [rvif::rvifs()] (`table`, with columns `RVIF` and `%`,
#'   one row per column of the design matrix, named after the model's
#'   variables).
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
rvif_diagnostics.glm <- function(model, center = FALSE, scale = "unit", tol = 1e-30, ...) {
  if (!identical(scale, "unit")) {
    warning(
      "'scale' has no effect on rvif_diagnostics(): rvif::rvifs() is always called with ",
      "ul = TRUE, which re-scales the design matrix to Euclidean unit length right before ",
      "computing the RVIF, overriding whatever 'scale' was requested here. Only 'center' ",
      "changes the result. Use scale = \"unit\" (the default) to avoid this warning.",
      call. = FALSE
    )
  }

  core <- .build_design(model, center = center, scale = scale)

  # ul = TRUE: rvifs() re-scales core$X_unit to Euclidean unit length right
  # before computing the RVIF, which is a no-op when scale = "unit" (the
  # default) and overrides any other 'scale' choice (see the warning above).
  tab <- rvif::rvifs(core$X_unit, ul = TRUE, intercept = core$has_intercept, tol = tol)
  if (is.null(tab)) {
    stop(
      "rvif::rvifs() could not compute the RVIF for this (weighted, transformed) design matrix ",
      "(see the message above from rvifs() itself, e.g. exact or near-exact multicollinearity). ",
      "Note that centering (center = TRUE) an IRLS-weighted design matrix that includes an ",
      "intercept is mathematically guaranteed to drop its rank by one for canonical-link GLMs; ",
      "use the default center = FALSE, or remove the exact collinearity from the data.",
      call. = FALSE
    )
  }
  rownames(tab) <- core$var_names

  structure(
    list(table = tab, var_names = core$var_names),
    class = "multicollglm_rvif"
  )
}

#' @rdname rvif_diagnostics
#' @export
rvif_diagnostics.glmnet <- function(model, x, y, family, s, weights = NULL, offset = NULL,
                                     center = FALSE, scale = "unit", tol = 1e-30, ...) {
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- rvif_diagnostics.glm(refit$model, center = center, scale = scale, tol = tol)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname rvif_diagnostics
#' @export
rvif_diagnostics.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                        center = FALSE, scale = "unit", tol = 1e-30, ...) {
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- rvif_diagnostics.glm(refit$model, center = center, scale = scale, tol = tol)
  out$lambda <- if (is.character(s)) model[[s]] else s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @export
print.multicollglm_rvif <- function(x, digits = 3, ...) {
  cat("Redefined Variance Inflation Factor (RVIF)\n")
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
