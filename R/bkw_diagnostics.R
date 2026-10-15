#' Full Belsley-Kuh-Welsch collinearity diagnostics
#'
#' Computes the condition index of every dimension and the Belsley, Kuh and
#' Welsch (1980) variance-decomposition proportion matrix, which identifies
#' which coefficients are involved in each collinearity relationship. The
#' matrix being decomposed is chosen with `method`, with the same four
#' definitions as [condition_number()]; by default it is the IRLS-weighted
#' design matrix scaled to unit length by column, without centering
#' (Weissfeld and Sereika, 1991).
#'
#' A dimension is flagged as problematic when its `sqrt(CN)` is
#' greater than or equal to `index_threshold` and at least two
#' coefficients have a variance proportion greater than or equal to
#' `proportion_threshold` on that dimension. The default `index_threshold
#' = 10` corresponds to a `CN` of 100 on the eigenvalue
#' (non-square-rooted) scale, since each component's `sqrt(CN)` is simply
#' the square root of that component's `CN`.
#'
#' The variance-decomposition proportions depend on how the columns are
#' scaled, so they are only comparable with those of other software when
#' the same matrix is decomposed. For instance, with `method = "X"` and a
#' model with an intercept the proportions are those of
#' `VisCollin::colldiag(mod, add.intercept = TRUE)`, since both decompose
#' the unit-length (unweighted) design matrix. The thresholds above are
#' meant for the unit-length definitions (`"WS"`, `"MP"` and `"X"`); with
#' `method = "RAW"` the condition indices are dominated by the units of the
#' variables.
#'
#' @inheritParams condition_number
#' @param index_threshold `sqrt(CN)` threshold used to flag a
#'   dimension as suspicious (default 10, i.e. a `CN` of 100
#'   before taking the square root).
#' @param proportion_threshold Variance-proportion threshold used to
#'   consider a coefficient involved in a suspicious dimension (default
#'   0.5).
#'
#' @return An object of class `multicollglm_bkw` with the eigenvalues,
#'   singular values, `sqrt(CN)` per component (`condition_index`), the
#'   proportions matrix (`proportions`, one row per coefficient and one
#'   column per component), the list of flagged dimensions (`flagged`) and
#'   the `method` used.
#'
#' @references
#' Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). *Regression
#' Diagnostics: Identifying Influential Data and Sources of Collinearity*.
#' Wiley.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
#' mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
#' y <- rgamma(n, shape = 5, rate = 5 / mu)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
#' bkw_diagnostics(mod)
#'
#' # Same decomposition on the original (unweighted) design matrix
#' bkw_diagnostics(mod, method = "X")
#'
#' @export
bkw_diagnostics <- function(model, ...) {
  UseMethod("bkw_diagnostics")
}

#' @export
bkw_diagnostics.glm <- function(model, index_threshold = 10, proportion_threshold = 0.5,
                                 method = c("WS", "MP", "X", "RAW"), ...) {
  .reject_removed_args(...)
  method <- .check_method(method)

  core <- .build_design(model, method = method)
  ev <- core$eigenvalues
  d <- core$singular_values
  V <- core$eigenvectors

  if (min(d) <= sqrt(.Machine$double.eps) * max(d)) {
    stop(
      "The (weighted, transformed) design matrix is exactly rank-deficient, so the ",
      "variance-decomposition proportions are undefined. This means there is perfect ",
      "collinearity among the columns of the design matrix: remove the redundant ",
      "variable(s) and fit the model again.",
      call. = FALSE
    )
  }

  condition_index <- max(d) / d

  # phi[j,k] = V[j,k]^2 / d[k]^2 ; pi[j,k] = phi[j,k] / sum_k phi[j,k]
  phi <- sweep(V^2, 2L, d^2, "/")
  row_sum <- rowSums(phi)
  proportions <- sweep(phi, 1L, row_sum, "/")
  dimnames(proportions) <- list(core$var_names, paste0("dim", seq_along(ev)))

  suspicious_dims <- which(condition_index >= index_threshold)
  flagged <- lapply(suspicious_dims, function(k) {
    vars <- core$var_names[proportions[, k] >= proportion_threshold]
    if (length(vars) >= 2) {
      list(dim = k, condition_index = condition_index[k], variables = vars)
    } else {
      NULL
    }
  })
  flagged <- Filter(Negate(is.null), flagged)

  structure(
    list(
      eigenvalues = ev,
      singular_values = d,
      condition_index = condition_index,
      proportions = proportions,
      index_threshold = index_threshold,
      proportion_threshold = proportion_threshold,
      flagged = flagged,
      var_names = core$var_names,
      method = method
    ),
    class = "multicollglm_bkw"
  )
}

#' @rdname bkw_diagnostics
#' @export
bkw_diagnostics.glmnet <- function(model, x, y, family, s, weights = NULL, offset = NULL,
                                    index_threshold = 10, proportion_threshold = 0.5,
                                    method = c("WS", "MP", "X", "RAW"), ...) {
  .reject_removed_args(...)
  method <- .check_method(method)
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- bkw_diagnostics.glm(refit$model, index_threshold = index_threshold,
                              proportion_threshold = proportion_threshold, method = method)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname bkw_diagnostics
#' @export
bkw_diagnostics.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                       index_threshold = 10, proportion_threshold = 0.5,
                                       method = c("WS", "MP", "X", "RAW"), ...) {
  .reject_removed_args(...)
  method <- .check_method(method)
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- bkw_diagnostics.glm(refit$model, index_threshold = index_threshold,
                              proportion_threshold = proportion_threshold, method = method)
  out$lambda <- if (is.character(s)) model[[s]] else s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @export
print.multicollglm_bkw <- function(x, digits = 3, ...) {
  tab <- cbind(
    "Eigenvalue" = x$eigenvalues,
    "Sing.value" = x$singular_values,
    "sqrt(CN)" = x$condition_index
  )
  rownames(tab) <- paste0("dim", seq_len(nrow(tab)))

  cat(sprintf("Collinearity diagnostics (Belsley-Kuh-Welsch) on %s [method = \"%s\"]\n",
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
  print(round(tab, digits))

  cat("\nVariance-decomposition proportions (row = coefficient, column = component):\n")
  print(round(x$proportions, digits))

  if (length(x$flagged) > 0) {
    cat(sprintf(
      "\n>> Possible collinearity problems (sqrt(CN) >= %g and proportion >= %g on >= 2 variables):\n",
      x$index_threshold, x$proportion_threshold
    ))
    for (p in x$flagged) {
      cat(sprintf("  - dim%d (sqrt(CN) = %.2f): %s\n", p$dim, p$condition_index, paste(p$variables, collapse = ", ")))
    }
  } else {
    cat(sprintf(
      "\nNo problems detected under the current thresholds (sqrt(CN) >= %g, proportion >= %g).\n",
      x$index_threshold, x$proportion_threshold
    ))
  }
  invisible(x)
}
