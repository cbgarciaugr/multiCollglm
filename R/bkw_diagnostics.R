#' Full Belsley-Kuh-Welsch collinearity diagnostics
#'
#' Computes, over the design matrix weighted by the IRLS weights at
#' convergence and, by default, scaled to unit length (without
#' centering), `sqrt(CN)` per component and the Belsley, Kuh and Welsch (1980)
#' variance-decomposition proportion matrix, which identifies which
#' coefficients are involved in each collinearity relationship.
#'
#' A dimension is flagged as problematic when its `sqrt(CN)` is
#' greater than or equal to `index_threshold` and at least two
#' coefficients have a variance proportion greater than or equal to
#' `proportion_threshold` on that dimension. The default `index_threshold
#' = 10` corresponds to a `CN` of 100 on the eigenvalue
#' (non-square-rooted) scale, since each component's `sqrt(CN)` is simply
#' the square root of that component's `CN`.
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
#'   column per component), and the list of flagged dimensions (`flagged`).
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
#' @export
bkw_diagnostics <- function(model, ...) {
  UseMethod("bkw_diagnostics")
}

#' @export
bkw_diagnostics.glm <- function(model, index_threshold = 10, proportion_threshold = 0.5,
                                 center = FALSE, scale = "unit", ...) {
  core <- .build_design(model, center = center, scale = scale)
  ev <- core$eigenvalues
  d <- core$singular_values
  V <- core$eigenvectors

  if (min(d) <= sqrt(.Machine$double.eps) * max(d)) {
    stop(
      "The (weighted, transformed) design matrix is exactly rank-deficient, so the ",
      "variance-decomposition proportions are undefined. Two known causes: (1) centering ",
      "(center = TRUE) an IRLS-weighted design matrix that includes an intercept, for some ",
      "family/link combinations (this is not a general canonical-link fact -- e.g. binomial ",
      "with the logit link is unaffected); and (2) the Gamma family with its canonical ",
      "(inverse) link, for which sqrt(w) * linear.predictors is exactly 1 for every ",
      "observation regardless of whether an intercept is present, so removing the intercept ",
      "does not avoid it there. This is exactly why Belsley, Kuh and Welsch recommend against ",
      "centering for collinearity diagnostics; use the default center = FALSE, or remove the ",
      "exact collinearity from the data.",
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
      var_names = core$var_names
    ),
    class = "multicollglm_bkw"
  )
}

#' @rdname bkw_diagnostics
#' @export
bkw_diagnostics.glmnet <- function(model, x, y, family, s, weights = NULL, offset = NULL,
                                    index_threshold = 10, proportion_threshold = 0.5,
                                    center = FALSE, scale = "unit", ...) {
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- bkw_diagnostics.glm(refit$model, index_threshold = index_threshold, proportion_threshold = proportion_threshold,
                              center = center, scale = scale)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname bkw_diagnostics
#' @export
bkw_diagnostics.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                       index_threshold = 10, proportion_threshold = 0.5,
                                       center = FALSE, scale = "unit", ...) {
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- bkw_diagnostics.glm(refit$model, index_threshold = index_threshold, proportion_threshold = proportion_threshold,
                              center = center, scale = scale)
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

  cat("Collinearity diagnostics (Belsley-Kuh-Welsch)\n")
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
