# Full Belsley-Kuh-Welsch collinearity diagnostics

Computes, over the design matrix weighted by the IRLS weights at
convergence and, by default, scaled to unit length (without centering),
`sqrt(CN)` per component and the Belsley, Kuh and Welsch (1980)
variance-decomposition proportion matrix, which identifies which
coefficients are involved in each collinearity relationship.

## Usage

``` r
bkw_diagnostics(model, ...)

# S3 method for class 'glmnet'
bkw_diagnostics(
  model,
  x,
  y,
  family,
  s,
  weights = NULL,
  offset = NULL,
  index_threshold = 10,
  proportion_threshold = 0.5,
  center = FALSE,
  scale = "unit",
  ...
)

# S3 method for class 'cv.glmnet'
bkw_diagnostics(
  model,
  x,
  y,
  family,
  s = "lambda.min",
  weights = NULL,
  offset = NULL,
  index_threshold = 10,
  proportion_threshold = 0.5,
  center = FALSE,
  scale = "unit",
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

- index_threshold:

  `sqrt(CN)` threshold used to flag a dimension as suspicious (default
  10, i.e. a `CN` of 100 before taking the square root).

- proportion_threshold:

  Variance-proportion threshold used to consider a coefficient involved
  in a suspicious dimension (default 0.5).

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
  intercept without issue). When it does occur, `bkw_diagnostics()` and
  [`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md)
  report it as an error rather than silent `NaN`s; this is the classical
  reason Belsley, Kuh and Welsch advise against centering for
  collinearity diagnostics. Ignored whenever `method` is not `NULL` (see
  below).

- scale:

  Either `"unit"` (default; each column is scaled to Euclidean/L2 unit
  length, as originally proposed and required for the classical
  `sqrt(CN)` thresholds and the RVIF formula) or any value accepted by
  [`scale()`](https://rdrr.io/r/base/scale.html)'s own `scale` argument:
  `TRUE` for root-mean-square scaling, `FALSE` for no scaling, or a
  numeric vector of per-column divisors. Ignored whenever `method` is
  not `NULL` (see below).

## Value

An object of class `multicollglm_bkw` with the eigenvalues, singular
values, `sqrt(CN)` per component (`condition_index`), the proportions
matrix (`proportions`, one row per coefficient and one column per
component), and the list of flagged dimensions (`flagged`).

## Details

A dimension is flagged as problematic when its `sqrt(CN)` is greater
than or equal to `index_threshold` and at least two coefficients have a
variance proportion greater than or equal to `proportion_threshold` on
that dimension. The default `index_threshold = 10` corresponds to a `CN`
of 100 on the eigenvalue (non-square-rooted) scale, since each
component's `sqrt(CN)` is simply the square root of that component's
`CN`.

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
y <- rgamma(n, shape = 5, rate = 5 / mu)
mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
bkw_diagnostics(mod)
#> Collinearity diagnostics (Belsley-Kuh-Welsch)
#> 
#>      Eigenvalue Sing.value sqrt(CN)
#> dim1      3.266      1.807    1.000
#> dim2      0.881      0.939    1.925
#> dim3      0.621      0.788    2.294
#> dim4      0.231      0.481    3.756
#> dim5      0.000      0.018  100.614
#> 
#> Variance-decomposition proportions (row = coefficient, column = component):
#>              dim1  dim2  dim3  dim4  dim5
#> (Intercept) 0.023 0.003 0.000 0.973 0.001
#> x1          0.000 0.000 0.000 0.000 1.000
#> x2          0.000 0.000 0.000 0.000 1.000
#> x3          0.029 0.002 0.809 0.141 0.018
#> x4          0.015 0.935 0.043 0.005 0.002
#> 
#> >> Possible collinearity problems (sqrt(CN) >= 10 and proportion >= 0.5 on >= 2 variables):
#>   - dim5 (sqrt(CN) = 100.61): x1, x2
```
