# Condition number of a model (collinearity diagnostic)

Computes the condition number over the design matrix weighted by the
IRLS weights at convergence and, by default, scaled to unit length by
column (without centering). This generalizes, to any GLM family and link
function, the classical condition number of linear regression.

## Usage

``` r
condition_number(model, ...)

# S3 method for class 'glmnet'
condition_number(
  model,
  x,
  y,
  family,
  s,
  weights = NULL,
  offset = NULL,
  center = FALSE,
  scale = "unit",
  method = NULL,
  ...
)

# S3 method for class 'cv.glmnet'
condition_number(
  model,
  x,
  y,
  family,
  s = "lambda.min",
  weights = NULL,
  offset = NULL,
  center = FALSE,
  scale = "unit",
  method = NULL,
  ...
)
```

## Arguments

- model:

  A `glm`, `glmnet` or `cv.glmnet` object.

- ...:

  Additional arguments passed to the methods (see `x`, `y`, `family`,
  `s`, `weights`, `offset` for `glmnet`/`cv.glmnet`).

- x:

  Predictor matrix (the same one used to fit the `glmnet` model), with
  `colnames`.

- y:

  Response variable used to fit the `glmnet` model.

- family:

  GLM family (a `family` object, e.g. `Gamma(link = "inverse")`) used to
  refit a `glm` on the active set of variables at the given lambda.

- s:

  Lambda value (numeric) for a `glmnet` object. For `cv.glmnet` it can
  also be `"lambda.min"` or `"lambda.1se"`.

- weights:

  Optional prior weights for the refit.

- offset:

  Optional offset for the refit.

- center:

  Either `FALSE` (default; no centering, as originally proposed and
  required to detect non-essential collinearity involving the intercept)
  or any value accepted by
  [`scale()`](https://rdrr.io/r/base/scale.html)'s own `center`
  argument: `TRUE` to center each column on its mean, or a numeric
  vector of per-column values to subtract. Centering an IRLS-weighted
  design matrix that includes an intercept can leave it exactly
  rank-deficient for *some* family/link combinations, but this is not a
  general canonical-link fact: it happens whenever
  `sqrt(w) * linear.predictors` is constant across all observations,
  which is a structural identity of the **Gamma family with its
  canonical (inverse) link** (`sqrt(w) = mu`, `eta = 1/mu`, so their
  product is exactly 1 for every observation) – and, notably for Gamma
  with that link, this holds **regardless of whether an intercept is
  present**, so dropping the intercept does not avoid it there (see
  `method = "OZ"` below). It does *not* occur in general for other
  canonical links such as the logit (binomial can be centered with an
  intercept without issue). When it does occur,
  [`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md)
  and
  [`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md)
  report it as an error rather than silent `NaN`s; this is the classical
  reason Belsley, Kuh and Welsch advise against centering for
  collinearity diagnostics. Ignored whenever `method` is not `NULL` (see
  below).

- scale:

  Either `"unit"` (default; each column is scaled to Euclidean/L2 unit
  length, as originally proposed and required for the classical
  condition-index thresholds and the RVIF formula) or any value accepted
  by [`scale()`](https://rdrr.io/r/base/scale.html)'s own `scale`
  argument: `TRUE` for root-mean-square scaling, `FALSE` for no scaling,
  or a numeric vector of per-column divisors. Ignored whenever `method`
  is not `NULL` (see below).

- method:

  Optional shortcut selecting one of five condition-number definitions;
  when supplied it overrides whatever `center`/`scale` were passed and
  labels the result accordingly (`NC_RAW`, `NC_MP`, `NC_MS`, `NC_WS` or
  `NC_OZ`):

  `"RAW"`

  :   No transformation at all. The condition number is computed
      directly on \\X'WX\\ – the IRLS-weighted design matrix (including
      the intercept column, if any), neither centered nor scaled
      (equivalent to `center = FALSE, scale = FALSE`). This is the
      "original, uncentered and unscaled data" baseline against which
      `"WS"`'s unit-length scaling can be compared: because the
      explanatory variables are typically on very different scales,
      `NC_RAW` is dominated by that scale disparity and is usually far
      larger than `NC_WS`, even when the underlying collinearity is the
      same – which is exactly why Weissfeld and Sereika (1991)
      introduced the unit-length rescaling in the first place.

  `"MP"`

  :   MacKinnon and Puterman (1989). The *original* explanatory
      variables (including the intercept column) are first rescaled to
      unit Euclidean length with
      [`multiColl::lu()`](https://rdrr.io/pkg/multiColl/man/lu.html),
      the GLM is **refit** on those rescaled variables (with no further
      intercept, since the rescaled intercept column already plays that
      role), and the condition number is then computed on that refit's
      IRLS-weighted information matrix with **no further
      transformation** (equivalent to `center = FALSE, scale = FALSE`
      applied to the refit). Refitting is mathematically redundant – a
      GLM's fitted values, and hence its IRLS weights at convergence,
      are invariant to any linear rescaling of the columns of `X` – but
      is carried out explicitly to mirror the original definition (and
      requires the multiColl package).

  `"MS"`

  :   Marx and Smith (1990). Like `"MP"`, a transform-then-fit method,
      but the explanatory variables (excluding the intercept column, if
      any – a constant column of ones cannot itself be centered and
      rescaled to unit length) are first **centered on their own mean
      and then rescaled to unit Euclidean length**, the GLM is **refit**
      on those centered-and-scaled variables (with an ordinary intercept
      re-added whenever the original model had one), and the condition
      number is computed on that refit's IRLS-weighted information
      matrix with **no further transformation**. Unlike `"MP"`'s pure
      rescaling, centering here generally shifts the fitted values
      themselves (the refit intercept absorbs the new baseline), so the
      refit is a genuinely different fit from `model`, exactly as Marx
      and Smith (1990) prescribe for their weighted multicollinearity
      diagnostics in logistic regression. For a model with **no**
      intercept, all columns are centered and the refit also has no
      intercept; this can fail to converge for some family/link
      combinations (e.g. Gamma with the inverse link, which needs
      positive fitted values) since centering removes the constant
      baseline a no-intercept model relies on.

  `"WS"`

  :   Weissfeld and Sereika (1991). The GLM is fit on the original
      variables and the IRLS-weighted information matrix is scaled to
      unit column length *afterwards*; equivalent to
      `center = FALSE, scale = "unit"`, the default of this function and
      of
      [`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md)/[`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md).

  `"OZ"`

  :   Ozkale (2019). The Weissfeld-Sereika procedure, but the
      IRLS-weighted design matrix is additionally *centered* before
      being scaled to unit length (equivalent to
      `center = TRUE, scale = "unit"`). If `model` has an intercept, its
      column is **dropped** from the design matrix before
      centering/scaling (the IRLS weights used are still those of the
      original, with-intercept fit): every worked example in Ozkale
      (2019, Sec. 5) fits the GLM *without* an intercept term, and for
      most family/link combinations dropping the intercept is exactly
      what avoids the exact rank-deficiency described under `center`
      above. **Exception:** for the Gamma family with its canonical
      (inverse) link, that rank-deficiency is structural regardless of
      whether an intercept is present or has been dropped
      (`sqrt(w) * linear.predictors` is exactly 1 for every observation
      either way) – this is exactly what happens in Ozkale (2019, Sec.
      5.3)'s own Gamma example, where the smallest eigenvalue she
      reports (4.0371e-17) is already at floating-point noise level;
      `NC_OZ` will then legitimately be reported as `Inf` with a
      warning, rather than the platform-dependent
      large-but-nominally-finite number her software happened to print.

  Leave as `NULL` (the default) to use `center`/`scale` directly,
  exactly as in earlier versions of this function.

## Value

An object of class `multicollglm_cn` with, among others, the eigenvalues
used (of \\X'WX\\ for `method = "RAW"`, of \\X\_{unit}'X\_{unit}\\ for
`"WS"`/`"OZ"`, or, for `method = "MP"`/`"MS"`, of the respective refit's
\\X\_{lu}'WX\_{lu}\\ or \\X\_{cu}'WX\_{cu}\\), the condition number on
the eigenvalue scale (`condition_number`), the classical condition index
on the singular-value scale (`condition_index`, i.e.
`sqrt(condition_number)`), and, when `method` was supplied, `method`
(`"RAW"`/`"MP"`/`"MS"`/`"WS"`/`"OZ"`) and `nc_label`
(`"NC_RAW"`/`"NC_MP"`/`"NC_MS"`/`"NC_WS"`/`"NC_OZ"`).

## References

Mackinnon, M.J. and Puterman, M.L. (1989). Collinearity in generalized
linear models. Communications in Statistics - Theory and Methods, 18(9),
3463-3472.
[doi:10.1080/03610928908830102](https://doi.org/10.1080/03610928908830102)

Marx, B.D. and Smith, E.P. (1990). Weighted multicollinearity in
logistic regression: diagnostics and biased estimation techniques with
an example from lake acidification. Canadian Journal of Fisheries and
Aquatic Sciences, 47(6), 1128-1135.
[doi:10.1139/f90-131](https://doi.org/10.1139/f90-131)

Weissfeld, L.A. and Sereika, S.M. (1991). A multicollinearity diagnostic
for generalized linear models. Communications in Statistics - Theory and
Methods, 20(4), 1183-1198.
[doi:10.1080/03610929108830558](https://doi.org/10.1080/03610929108830558)

Ozkale, M.R. (2019). The red indicator and corrected VIFs in generalized
linear models. Communications in Statistics - Simulation and
Computation.
[doi:10.1080/03610918.2019.1639740](https://doi.org/10.1080/03610918.2019.1639740)

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
y <- rgamma(n, shape = 5, rate = 5 / mu)
mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
condition_number(mod)
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 3.2659 0.8814 0.6208 0.2315 0.0003
#> 
#> Condition number (eigenvalue scale): 10123.1410
#> Condition index (classical, sqrt): 100.6138
condition_number(mod, method = "RAW") # X'WX, no centering or scaling
#> Eigenvalues of X'WX (decreasing order):
#> [1] 25684.1167  4325.2941  3318.8684  1000.0895     3.3983
#> 
#> Condition number (NC_RAW) (eigenvalue scale): 7557.9499
#> Condition index (NC_RAW) (classical, sqrt): 86.9365
condition_number(mod, method = "WS") # same numbers, labeled NC_WS
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 3.2659 0.8814 0.6208 0.2315 0.0003
#> 
#> Condition number (NC_WS) (eigenvalue scale): 10123.1410
#> Condition index (NC_WS) (classical, sqrt): 100.6138
condition_number(mod, method = "OZ")
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 2.4718 0.8806 0.6473 0.0003
#> 
#> Condition number (NC_OZ) (eigenvalue scale): 7201.4286
#> Condition index (NC_OZ) (classical, sqrt): 84.8612
condition_number(mod, method = "MS") # center + unit length, then refit
#> Eigenvalues of X_cu'WX_cu (decreasing order):
#> [1] 3816.0817   42.5254   18.0071   13.6069    0.0198
#> 
#> Condition number (NC_MS) (eigenvalue scale): 192963.4257
#> Condition index (NC_MS) (classical, sqrt): 439.2760
if (FALSE) { # \dontrun{
condition_number(mod, method = "MP") # requires the 'multiColl' package
} # }
```
