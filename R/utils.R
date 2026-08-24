# Internal functions, not exported.

# Given an already-fitted glm, builds the design matrix weighted by the
# IRLS weights at convergence and transformed via 'center'/'scale', and
# computes the spectral decomposition of X_unit' X_unit needed for the
# condition number and the Belsley-Kuh-Welsch diagnostics.
#
# 'center' and 'scale' accept the same values as base::scale()'s own
# 'center'/'scale' arguments (a logical or a numeric vector of per-column
# values), with one addition: scale = "unit" (the default) scales each
# (possibly centered) column to Euclidean/L2 unit length -- the
# transformation originally proposed, and the one required for the
# classical condition-index thresholds and the RVIF formula. Any other
# 'scale' value is delegated to base::scale() as-is (e.g. scale = TRUE
# for root-mean-square scaling).
.build_design <- function(mod, center = FALSE, scale = "unit") {
  if (!inherits(mod, "glm")) {
    stop("A 'glm' object was expected.", call. = FALSE)
  }

  X <- stats::model.matrix(mod)
  w <- mod$weights

  if (is.null(w) || length(w) != nrow(X)) {
    stop(
      "The model does not contain valid IRLS weights ('mod$weights'). ",
      "Fit the model with glm() and do not reuse a hand-modified object.",
      call. = FALSE
    )
  }
  if (any(w < 0) || anyNA(w)) {
    stop("The model's IRLS weights contain negative or NA values.", call. = FALSE)
  }

  # X %*% W^(1/2): equivalent to diag(sqrt(w)) %*% X without building the diagonal matrix
  Xw <- X * sqrt(w)

  if (identical(scale, "unit")) {
    Xc <- if (identical(center, FALSE)) {
      Xw
    } else if (isTRUE(center)) {
      sweep(Xw, 2L, colMeans(Xw), "-")
    } else {
      sweep(Xw, 2L, center, "-")
    }

    col_norms <- sqrt(colSums(Xc^2))
    if (any(col_norms == 0)) {
      zero_cols <- colnames(X)[col_norms == 0]
      stop(
        sprintf(
          "Columns [%s] have zero norm after the requested transformation; they cannot be scaled to unit length.",
          paste(zero_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    X_unit <- sweep(Xc, 2L, col_norms, "/")
  } else {
    X_unit <- base::scale(Xw, center = center, scale = scale)
    if (anyNA(X_unit)) {
      stop(
        "scale() produced NA/NaN values for the requested 'center'/'scale'; ",
        "check for zero-variance or zero-scale columns.",
        call. = FALSE
      )
    }
    attributes(X_unit)[c("scaled:center", "scaled:scale")] <- NULL
  }

  M <- crossprod(X_unit)
  eig <- eigen(M, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)

  eigenvalues <- eig$values[ord]
  eigenvectors <- eig$vectors[, ord, drop = FALSE]
  colnames(eigenvectors) <- paste0("dim", seq_len(ncol(eigenvectors)))
  rownames(eigenvectors) <- colnames(X)

  list(
    model = mod,
    X = X,
    X_unit = X_unit,
    eigenvalues = eigenvalues,
    eigenvectors = eigenvectors,
    singular_values = sqrt(pmax(eigenvalues, 0)),
    var_names = colnames(X),
    has_intercept = identical(colnames(X)[1], "(Intercept)")
  )
}

# Given a coefficient vector/matrix from glmnet (coef(fit, s = ...)),
# identifies the active set (coefficients != 0, excluding the intercept)
# and refits a standard glm() (IRLS) on those variables so the diagnostics
# from .build_design() can be applied.
.refit_active_glm <- function(coefs, x, y, family, weights = NULL, offset = NULL) {
  coefs <- as.matrix(coefs)
  cf <- stats::setNames(as.numeric(coefs), rownames(coefs))

  active <- names(cf)[cf != 0]
  active <- setdiff(active, "(Intercept)")

  if (length(active) == 0) {
    stop(
      "No coefficient in the active set is different from zero at the given lambda; ",
      "the model cannot be refit.",
      call. = FALSE
    )
  }
  if (is.null(colnames(x))) {
    stop("'x' must have column names (colnames) matching those of the glmnet model.", call. = FALSE)
  }
  if (!all(active %in% colnames(x))) {
    missing_cols <- setdiff(active, colnames(x))
    stop(
      sprintf(
        "The columns of 'x' do not include the active variables [%s]. Pass the same 'x' matrix used to fit the glmnet model.",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  frame <- data.frame(.y = y, x[, active, drop = FALSE], check.names = FALSE)
  active_formula <- stats::reformulate(
    termlabels = paste0("`", active, "`"),
    response = ".y"
  )

  args <- list(formula = active_formula, data = frame, family = family)
  if (!is.null(weights)) args$weights <- weights
  if (!is.null(offset)) args$offset <- offset

  model <- do.call(stats::glm, args)

  list(
    model = model,
    active = active,
    dropped = setdiff(colnames(x), active)
  )
}

.check_family <- function(family) {
  if (missing(family) || is.null(family)) {
    stop(
      "You must supply 'family' (e.g. Gamma(link = \"inverse\")) to refit the glm on the glmnet active set.",
      call. = FALSE
    )
  }
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) {
    stop("'family' must be an object of class 'family' (see ?family), e.g. Gamma(link = \"inverse\").", call. = FALSE)
  }
  family
}
