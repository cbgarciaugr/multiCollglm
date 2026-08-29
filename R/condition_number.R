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
#'   vector of per-column values to subtract. Centering an IRLS-weighted
#'   design matrix that includes an intercept can leave it exactly
#'   rank-deficient for *some* family/link combinations, but this is not a
#'   general canonical-link fact: it happens whenever `sqrt(w) *
#'   linear.predictors` is constant across all observations, which is a
#'   structural identity of the **Gamma family with its canonical (inverse)
#'   link** (`sqrt(w) = mu`, `eta = 1/mu`, so their product is exactly 1
#'   for every observation) -- and, notably for Gamma with that link, this
#'   holds **regardless of whether an intercept is present**, so dropping
#'   the intercept does not avoid it there (see `method = "OZ"` below). It
#'   does *not* occur in general for other canonical links such as the
#'   logit (binomial can be centered with an intercept without issue).
#'   When it does occur, [bkw_diagnostics()] and [rvif_diagnostics()]
#'   report it as an error rather than silent `NaN`s; this is the classical
#'   reason Belsley, Kuh and Welsch advise against centering for
#'   collinearity diagnostics. Ignored whenever `method` is not `NULL`
#'   (see below).
#' @param scale Either `"unit"` (default; each column is scaled to
#'   Euclidean/L2 unit length, as originally proposed and required for the
#'   classical condition-index thresholds and the RVIF formula) or any
#'   value accepted by [scale()]'s own `scale` argument: `TRUE` for
#'   root-mean-square scaling, `FALSE` for no scaling, or a numeric vector
#'   of per-column divisors. Ignored whenever `method` is not `NULL` (see
#'   below).
#' @param method Optional shortcut selecting one of five condition-number
#'   definitions; when supplied it overrides whatever `center`/`scale` were
#'   passed and labels the result accordingly (`NC_RAW`, `NC_MP`, `NC_MS`,
#'   `NC_WS` or `NC_OZ`):
#'   \describe{
#'     \item{`"RAW"`}{No transformation at all. The condition number is
#'       computed directly on \eqn{X'WX} -- the IRLS-weighted design matrix
#'       (including the intercept column, if any), neither centered nor
#'       scaled (equivalent to `center = FALSE, scale = FALSE`). This is
#'       the "original, uncentered and unscaled data" baseline against
#'       which `"WS"`'s unit-length scaling can be compared: because the
#'       explanatory variables are typically on very different scales,
#'       `NC_RAW` is dominated by that scale disparity and is usually far
#'       larger than `NC_WS`, even when the underlying collinearity is the
#'       same -- which is exactly why Weissfeld and Sereika (1991)
#'       introduced the unit-length rescaling in the first place.}
#'     \item{`"MP"`}{MacKinnon and Puterman (1989). The *original*
#'       explanatory variables (including the intercept column) are first
#'       rescaled to unit Euclidean length with [multiColl::lu()], the GLM
#'       is **refit** on those rescaled variables (with no further
#'       intercept, since the rescaled intercept column already plays that
#'       role), and the condition number is then computed on that refit's
#'       IRLS-weighted information matrix with **no further
#'       transformation** (equivalent to `center = FALSE, scale = FALSE`
#'       applied to the refit). Refitting is mathematically redundant --
#'       a GLM's fitted values, and hence its IRLS weights at convergence,
#'       are invariant to any linear rescaling of the columns of `X` -- but
#'       is carried out explicitly to mirror the original definition (and
#'       requires the \CRANpkg{multiColl} package).}
#'     \item{`"MS"`}{Marx and Smith (1990). Like `"MP"`, a transform-then-fit
#'       method, but the explanatory variables (excluding the intercept
#'       column, if any -- a constant column of ones cannot itself be
#'       centered and rescaled to unit length) are first **centered on their
#'       own mean and then rescaled to unit Euclidean length**, the GLM is
#'       **refit** on those centered-and-scaled variables (with an ordinary
#'       intercept re-added whenever the original model had one), and the
#'       condition number is computed on that refit's IRLS-weighted
#'       information matrix with **no further transformation**. Unlike
#'       `"MP"`'s pure rescaling, centering here generally shifts the fitted
#'       values themselves (the refit intercept absorbs the new baseline),
#'       so the refit is a genuinely different fit from `model`, exactly as
#'       Marx and Smith (1990) prescribe for their weighted multicollinearity
#'       diagnostics in logistic regression.}
#'     \item{`"WS"`}{Weissfeld and Sereika (1991). The GLM is fit on the
#'       original variables and the IRLS-weighted information matrix is
#'       scaled to unit column length *afterwards*; equivalent to
#'       `center = FALSE, scale = "unit"`, the default of this function and
#'       of [bkw_diagnostics()]/[rvif_diagnostics()].}
#'     \item{`"OZ"`}{Ozkale (2019). The Weissfeld-Sereika procedure, but the
#'       IRLS-weighted design matrix is additionally *centered* before
#'       being scaled to unit length (equivalent to `center = TRUE, scale =
#'       "unit"`). If `model` has an intercept, its column is **dropped**
#'       from the design matrix before centering/scaling (the IRLS weights
#'       used are still those of the original, with-intercept fit): every
#'       worked example in Ozkale (2019, Sec. 5) fits the GLM *without* an
#'       intercept term, and for most family/link combinations dropping the
#'       intercept is exactly what avoids the exact rank-deficiency
#'       described under `center` above. **Exception:** for the Gamma
#'       family with its canonical (inverse) link, that rank-deficiency is
#'       structural regardless of whether an intercept is present or has
#'       been dropped (`sqrt(w) * linear.predictors` is exactly 1 for every
#'       observation either way) -- this is exactly what happens in Ozkale
#'       (2019, Sec. 5.3)'s own Gamma example, where the smallest eigenvalue
#'       she reports (4.0371e-17) is already at floating-point noise level;
#'       `NC_OZ` will then legitimately be reported as `Inf` with a warning,
#'       rather than the platform-dependent large-but-nominally-finite
#'       number her software happened to print.}
#'   }
#'   Leave as `NULL` (the default) to use `center`/`scale` directly, exactly
#'   as in earlier versions of this function.
#' @param ... Additional arguments passed to the methods (see `x`, `y`,
#'   `family`, `s`, `weights`, `offset` for `glmnet`/`cv.glmnet`).
#'
#' @return An object of class `multicollglm_cn` with, among others, the
#'   eigenvalues used (of \eqn{X'WX} for `method = "RAW"`, of
#'   \eqn{X_{unit}'X_{unit}} for `"WS"`/`"OZ"`, or, for `method = "MP"`/`"MS"`,
#'   of the respective refit's \eqn{X_{lu}'WX_{lu}} or \eqn{X_{cu}'WX_{cu}}),
#'   the condition number on the
#'   eigenvalue scale (`condition_number`), the classical condition index
#'   on the singular-value scale (`condition_index`, i.e.
#'   `sqrt(condition_number)`), and, when `method` was supplied, `method`
#'   (`"RAW"`/`"MP"`/`"MS"`/`"WS"`/`"OZ"`) and `nc_label`
#'   (`"NC_RAW"`/`"NC_MP"`/`"NC_MS"`/`"NC_WS"`/`"NC_OZ"`).
#'
#' @references
#' Mackinnon, M.J. and Puterman, M.L. (1989). Collinearity in generalized
#' linear models. Communications in Statistics - Theory and Methods, 18(9),
#' 3463-3472. \doi{10.1080/03610928908830102}
#'
#' Marx, B.D. and Smith, E.P. (1990). Weighted multicollinearity in logistic
#' regression: diagnostics and biased estimation techniques with an example
#' from lake acidification. Canadian Journal of Fisheries and Aquatic
#' Sciences, 47(6), 1128-1135. \doi{10.1139/f90-131}
#'
#' Weissfeld, L.A. and Sereika, S.M. (1991). A multicollinearity diagnostic
#' for generalized linear models. Communications in Statistics - Theory and
#' Methods, 20(4), 1183-1198. \doi{10.1080/03610929108830558}
#'
#' Ozkale, M.R. (2019). The red indicator and corrected VIFs in generalized
#' linear models. Communications in Statistics - Simulation and
#' Computation. \doi{10.1080/03610918.2019.1639740}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
#' mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
#' y <- rgamma(n, shape = 5, rate = 5 / mu)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
#' condition_number(mod)
#' condition_number(mod, method = "RAW") # X'WX, no centering or scaling
#' condition_number(mod, method = "WS") # same numbers, labeled NC_WS
#' condition_number(mod, method = "OZ")
#' condition_number(mod, method = "MS") # center + unit length, then refit
#' \dontrun{
#' condition_number(mod, method = "MP") # requires the 'multiColl' package
#' }
#'
#' @export
condition_number <- function(model, ...) {
  UseMethod("condition_number")
}

#' @export
condition_number.glm <- function(model, center = FALSE, scale = "unit", method = NULL, ...) {
  drop_intercept <- FALSE

  if (!is.null(method)) {
    method <- match.arg(method, c("RAW", "MP", "MS", "WS", "OZ"))

    if (identical(method, "MP")) {
      refit <- .refit_lu_glm(model)
      out <- condition_number.glm(refit, center = FALSE, scale = FALSE)
      out$method <- "MP"
      out$nc_label <- "NC_MP"
      return(out)
    }

    if (identical(method, "MS")) {
      refit <- .refit_center_unit_glm(model)
      out <- condition_number.glm(refit, center = FALSE, scale = FALSE)
      out$method <- "MS"
      out$nc_label <- "NC_MS"
      return(out)
    }

    if (identical(method, "RAW")) {
      center <- FALSE
      scale <- FALSE
      drop_intercept <- FALSE
    } else {
      center <- identical(method, "OZ")
      scale <- "unit"
      drop_intercept <- identical(method, "OZ")
    }
  }

  core <- .build_design(model, center = center, scale = scale, drop_intercept = drop_intercept)
  ev_raw <- core$eigenvalues

  # crossprod() of a real matrix is symmetric positive semi-definite, so its
  # true eigenvalues are never negative; a value <= 0 here can only be
  # floating-point noise around an exactly zero eigenvalue. Two known,
  # unrelated causes: (1) centering an IRLS-weighted design matrix that
  # still includes an intercept, for family/link combinations where that
  # is degenerate (method = "OZ" drops the intercept column precisely to
  # avoid this one); and (2) the Gamma family with its canonical (inverse)
  # link, for which sqrt(w) * linear.predictors is exactly 1 for every
  # observation regardless of whether an intercept is present or has been
  # dropped -- method = "OZ" cannot avoid this second cause, since it is a
  # property of the fitted coefficients on whatever columns the model
  # actually uses, not of the intercept specifically (Ozkale 2019, Sec.
  # 5.3's own Gamma example hits exactly this, with her reported smallest
  # eigenvalue, 4.0371e-17, already at floating-point noise level). Report
  # the mathematically correct Inf in either case instead of a nonsensical
  # negative condition number, or a platform-dependent finite number that
  # is really just noise, from dividing by a tiny floating-point value.
  if (min(ev_raw) <= 0) {
    warning(
      "The (weighted, transformed) design matrix is exactly rank-deficient (a zero or ",
      "floating-point-negative eigenvalue was found); the condition number is mathematically ",
      "infinite and is reported as Inf. Two known causes: (1) centering (center = TRUE) an ",
      "IRLS-weighted design matrix that still includes an intercept, for some family/link ",
      "combinations (method = \"OZ\" drops the intercept column to avoid this); and (2) the ",
      "Gamma family with its canonical (inverse) link, where sqrt(w) * linear.predictors is ",
      "exactly 1 for every observation regardless of the intercept, so method = \"OZ\" cannot ",
      "avoid this second cause (see Ozkale 2019, Sec. 5.3, whose own Gamma example hits exactly ",
      "this and reports a large but essentially arbitrary finite number for the same reason).",
      call. = FALSE
    )
  }
  ev <- pmax(ev_raw, 0)
  cn <- if (min(ev) <= 0) Inf else max(ev) / min(ev)

  structure(
    list(
      eigenvalues = ev,
      condition_number = cn,
      condition_index = sqrt(cn),
      var_names = core$var_names,
      method = method,
      nc_label = if (is.null(method)) NULL else paste0("NC_", method)
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
                                     center = FALSE, scale = "unit", method = NULL, ...) {
  if (missing(s) || is.null(s)) {
    stop("You must supply 's' (the lambda value) for a 'glmnet' object.", call. = FALSE)
  }
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- condition_number.glm(refit$model, center = center, scale = scale, method = method)
  out$lambda <- s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @rdname condition_number
#' @export
condition_number.cv.glmnet <- function(model, x, y, family, s = "lambda.min", weights = NULL, offset = NULL,
                                        center = FALSE, scale = "unit", method = NULL, ...) {
  family <- .check_family(family)
  cf <- stats::coef(model, s = s)
  refit <- .refit_active_glm(cf, x = x, y = y, family = family, weights = weights, offset = offset)

  out <- condition_number.glm(refit$model, center = center, scale = scale, method = method)
  out$lambda <- if (is.character(s)) model[[s]] else s
  out$active_vars <- refit$active
  out$dropped_vars <- refit$dropped
  out
}

#' @export
print.multicollglm_cn <- function(x, digits = 4, ...) {
  mat_label <- switch(
    if (is.null(x$method)) "" else x$method,
    "RAW" = "X'WX",
    "MP"  = "X_lu'WX_lu",
    "MS"  = "X_cu'WX_cu",
    "X_unit'X_unit" # default for "WS", "OZ" and method = NULL
  )
  cat(sprintf("Eigenvalues of %s (decreasing order):\n", mat_label))
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

  label <- if (is.null(x$nc_label)) "" else sprintf(" (%s)", x$nc_label)
  cat(sprintf("\nCondition number%s (eigenvalue scale): %.*f\n", label, digits, x$condition_number))
  cat(sprintf("Condition index%s (classical, sqrt): %.*f\n", label, digits, x$condition_index))
  invisible(x)
}
