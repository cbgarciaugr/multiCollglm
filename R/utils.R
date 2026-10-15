# Internal functions, not exported.

# The four condition-number definitions supported by condition_number(),
# bkw_diagnostics() and rvif_diagnostics(). The first one is the default.
.methods <- c("WS", "MP", "X", "RAW")

# Validates and normalizes the 'method' argument. 'OZ' (Ozkale 2021) and
# 'MS' (Marx and Smith 1990) existed up to version 0.3.0 and were removed in
# 0.4.0; they get a specific message instead of a generic match.arg() error.
.check_method <- function(method) {
  if (is.character(method) && length(method) == 1L && method %in% c("OZ", "MS")) {
    stop(
      sprintf("method = \"%s\" was removed in multiCollglm 0.4.0. ", method),
      "Available methods: ", paste0("\"", .methods, "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }
  match.arg(method, .methods)
}

# 'center' and 'scale' were removed in 0.4.0 (the transformation is now
# fixed by 'method'). Because the S3 methods take '...', passing them would
# otherwise be silently swallowed, which would silently change results for
# code written against 0.3.0; fail loudly instead.
.reject_removed_args <- function(...) {
  removed <- intersect(names(list(...)), c("center", "scale"))
  if (length(removed) > 0L) {
    stop(
      sprintf("Argument(s) %s were removed in multiCollglm 0.4.0. ",
              paste0("'", removed, "'", collapse = ", ")),
      "Choose the transformation with 'method' instead: ",
      paste0("\"", .methods, "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Divides each column of M by its Euclidean norm (unit-length scaling, no
# centering).
.unit_length <- function(M) {
  norms <- sqrt(colSums(M^2))
  if (any(norms == 0)) {
    stop(
      sprintf(
        "Columns [%s] have zero norm; they cannot be scaled to unit length.",
        paste(colnames(M)[norms == 0], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sweep(M, 2L, norms, "/")
}

# MacKinnon and Puterman (1989): the columns of X (intercept included) are
# scaled to unit length BEFORE fitting, the GLM is refitted on the scaled
# matrix (with no further intercept: the scaled column of ones already plays
# that role) and the information matrix is read off that refit. The refit
# keeps the response, family, prior weights and offset of 'mod'. Returns the
# refit's IRLS weights together with the scaled matrix.
.refit_unit_glm <- function(mod) {
  X <- stats::model.matrix(mod)
  if (ncol(X) < 2L) {
    stop("method = \"MP\" needs at least two columns in the design matrix.", call. = FALSE)
  }
  if (is.null(mod$y)) {
    stop(
      "The 'glm' object does not retain the response ('model$y' is NULL); refit it with ",
      "y = TRUE (the default in glm()) before using method = \"MP\".",
      call. = FALSE
    )
  }
  X_ul <- .unit_length(X)
  fit <- stats::glm.fit(
    x = X_ul, y = mod$y, weights = mod$prior.weights,
    offset = mod$offset, family = mod$family, intercept = FALSE
  )
  list(X_ul = X_ul, weights = fit$weights)
}

# Given an already-fitted glm, builds the matrix M whose cross-product M'M is
# the one the diagnostic is computed on, according to 'method':
#
#   "WS"  unit-length scaling of W^(1/2) X          -> (X'WX)_ul = S X'WX S
#   "MP"  W*^(1/2) X_ul (X scaled before refitting) -> X_ul' W X_ul
#   "X"   unit-length scaling of X (no weights)     -> X_ul' X_ul
#   "RAW" W^(1/2) X (no transformation)             -> X'WX
#
# and its spectral decomposition, needed for the condition number and the
# Belsley-Kuh-Welsch diagnostics. Nothing is ever centered, so the intercept
# column, if any, always takes part in the diagnostic.
.build_design <- function(mod, method = "WS") {
  if (!inherits(mod, "glm")) {
    stop("A 'glm' object was expected.", call. = FALSE)
  }
  method <- .check_method(method)

  X <- stats::model.matrix(mod)

  if (method == "X") {
    M <- .unit_length(X)
  } else if (method == "MP") {
    M <- .refit_unit_glm(mod)
    M <- M$X_ul * sqrt(M$weights)
  } else {
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
    Xw <- X * sqrt(w) # W^(1/2) X without building the diagonal matrix
    M <- if (method == "WS") .unit_length(Xw) else Xw
  }
  colnames(M) <- colnames(X)

  eig <- eigen(crossprod(M), symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)

  eigenvalues <- eig$values[ord]
  eigenvectors <- eig$vectors[, ord, drop = FALSE]
  colnames(eigenvectors) <- paste0("dim", seq_len(ncol(eigenvectors)))
  rownames(eigenvectors) <- colnames(X)

  list(
    model = mod,
    method = method,
    X = X,
    M = M,
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
