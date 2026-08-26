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
.build_design <- function(mod, center = FALSE, scale = "unit", drop_intercept = FALSE) {
  if (!inherits(mod, "glm")) {
    stop("A 'glm' object was expected.", call. = FALSE)
  }

  X <- stats::model.matrix(mod)

  if (isTRUE(drop_intercept) && "(Intercept)" %in% colnames(X)) {
    # Ozkale (2019) centers and scales the IRLS-weighted design matrix, but
    # every one of her worked examples fits the GLM *without* an intercept
    # term to begin with. Centering a weighted design matrix that includes
    # an intercept is mathematically guaranteed to drop its rank by exactly
    # one for a canonical-link GLM (see the warning below in
    # condition_number.glm()), which her examples never exhibit. To match
    # her procedure without requiring the user to refit the model, the
    # intercept column is simply excluded here before centering/scaling;
    # the IRLS weights used are still those of the original (with-intercept)
    # fit passed in by the caller.
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    if (ncol(X) == 0L) {
      stop(
        "Dropping the intercept (method = \"OZ\") left no explanatory variables to diagnose; ",
        "this model has no other terms.",
        call. = FALSE
      )
    }
  }

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

# Refits a glm with the same response, family, prior weights and offset as
# 'mod', but on its explanatory variables (including the intercept column)
# each rescaled to Euclidean unit length via multiColl::lu(), with no
# further intercept added (the rescaled intercept column already plays
# that role). This is the MacKinnon and Puterman (1989) transformation:
# scale the columns of X to unit length *before* fitting, then read
# collinearity directly off X'W(beta)X with no further transformation --
# used by condition_number(..., method = "MP").
#
# Note this refit is mathematically guaranteed to reproduce the same
# fitted values (and hence the same IRLS weights) as 'mod': a GLM's fitted
# values are invariant under any invertible linear reparametrisation of the
# columns of X (here, dividing each column by a positive constant), since
# the linear predictor X %*% beta is unchanged by rescaling a column and
# inversely rescaling its coefficient. The refit is nonetheless carried out
# explicitly, rather than reusing mod$weights directly, to mirror the
# original published definition exactly and to stay correct even if that
# invariance were ever broken in practice (e.g. non-convergence).
.refit_lu_glm <- function(mod) {
  if (!inherits(mod, "glm")) {
    stop("A 'glm' object was expected.", call. = FALSE)
  }
  if (!requireNamespace("multiColl", quietly = TRUE)) {
    stop(
      "Computing the MacKinnon-Puterman (1989) condition number (method = \"MP\") requires the ",
      "'multiColl' package, which supplies lu(), the unit-length transformation used to rescale ",
      "the explanatory variables before refitting the model. Please install it.",
      call. = FALSE
    )
  }

  X <- stats::model.matrix(mod)
  if (ncol(X) < 2L) {
    stop("method = \"MP\" needs at least two columns in the design matrix.", call. = FALSE)
  }

  X_lu <- multiColl::lu(X)
  colnames(X_lu) <- colnames(X)

  y <- mod$y
  if (is.null(y)) {
    stop(
      "The 'glm' object does not retain the response ('model$y' is NULL); refit it with y = TRUE ",
      "(the default in glm()) before calling condition_number(model, method = \"MP\").",
      call. = FALSE
    )
  }

  frame <- data.frame(.y = y, X_lu, check.names = FALSE)
  formula_lu <- stats::reformulate(
    termlabels = paste0("`", colnames(X_lu), "`"),
    response = ".y",
    intercept = FALSE
  )

  args <- list(formula = formula_lu, data = frame, family = mod$family)

  w0 <- mod$prior.weights
  if (!is.null(w0) && !isTRUE(all.equal(unname(w0), rep(1, length(w0))))) {
    args$weights <- w0
  }
  off <- mod$offset
  if (!is.null(off) && !isTRUE(all.equal(unname(off), rep(0, length(off))))) {
    args$offset <- off
  }

  do.call(stats::glm, args)
}
