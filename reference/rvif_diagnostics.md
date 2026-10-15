# Redefined Variance Inflation Factor (RVIF) diagnostics

Computes the Redefined Variance Inflation Factor (RVIF) of Salmeron,
Garcia and Garcia (2025), via
[`rvif::rvifs()`](https://rdrr.io/pkg/rvif/man/rvifs.html), on the same
matrix used by
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)
and
[`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md),
chosen with `method`. By default this is the IRLS-weighted design matrix
scaled to unit length by column, without centering (Weissfeld and
Sereika, 1991). Unlike the classical VIF, the RVIF is able to detect
both essential collinearity (near-linear relationships among the
regressors) and non-essential collinearity (near-linear relationships
between the intercept and one or more regressors), and also decomposes,
for each variable, the percentage of near collinearity it is responsible
for.

## Usage

``` r
rvif_diagnostics(model, ...)

# S3 method for class 'glmnet'
rvif_diagnostics(
  model,
  x,
  y,
  family,
  s,
  weights = NULL,
  offset = NULL,
  method = c("WS", "MP", "X", "RAW"),
  tol = 1e-30,
  ...
)

# S3 method for class 'cv.glmnet'
rvif_diagnostics(
  model,
  x,
  y,
  family,
  s = "lambda.min",
  weights = NULL,
  offset = NULL,
  method = c("WS", "MP", "X", "RAW"),
  tol = 1e-30,
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

- method:

  One of `"WS"` (default), `"MP"`, `"X"` or `"RAW"`. See Details.

- tol:

  Tolerance used by
  [`rvif::rvifs()`](https://rdrr.io/pkg/rvif/man/rvifs.html) to decide
  whether the system is computationally singular (default `1e-30`).

## Value

An object of class `multicollglm_rvif` wrapping the data frame returned
by [`rvif::rvifs()`](https://rdrr.io/pkg/rvif/man/rvifs.html) (`table`,
with columns `RVIF` and `%`, one row per column of the design matrix,
named after the model's variables) and the `method` used.

## Details

[`rvif::rvifs()`](https://rdrr.io/pkg/rvif/man/rvifs.html) is always
called with `ul = TRUE`, so it re-scales the matrix to Euclidean unit
length right before computing the RVIF. As a consequence `"WS"`, `"MP"`
and `"RAW"`, which only differ in how the weighted matrix is scaled,
give the same RVIF; what changes the result is whether the weights take
part at all, i.e. `"X"` (the original design matrix, ignoring the
weights) against the other three.

## References

Salmeron, R., Garcia, C.B. and Garcia, J. (2025). A redefined Variance
Inflation Factor: overcoming the limitations of the Variance Inflation
Factor. Computational Economics, 65, 337-363.
[doi:10.1007/s10614-024-10575-8](https://doi.org/10.1007/s10614-024-10575-8)

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
y <- rgamma(n, shape = 5, rate = 5 / mu)
mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
rvif_diagnostics(mod)
#> Redefined Variance Inflation Factor (RVIF) on (X'WX)_ul [method = "WS"]
#> 
#>                 RVIF      %
#> (Intercept)    3.288 69.589
#> x1          1554.971 99.936
#> x2          1546.036 99.935
#> x3             1.581 36.758
#> x4             1.111  9.960
```
