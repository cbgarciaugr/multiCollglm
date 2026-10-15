# Full Belsley-Kuh-Welsch collinearity diagnostics

Computes the condition index of every dimension and the Belsley, Kuh and
Welsch (1980) variance-decomposition proportion matrix, which identifies
which coefficients are involved in each collinearity relationship. The
matrix being decomposed is chosen with `method`, with the same four
definitions as
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md);
by default it is the IRLS-weighted design matrix scaled to unit length
by column, without centering (Weissfeld and Sereika, 1991).

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
  method = c("WS", "MP", "X", "RAW"),
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
  method = c("WS", "MP", "X", "RAW"),
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

- method:

  One of `"WS"` (default), `"MP"`, `"X"` or `"RAW"`. See Details.

## Value

An object of class `multicollglm_bkw` with the eigenvalues, singular
values, `sqrt(CN)` per component (`condition_index`), the proportions
matrix (`proportions`, one row per coefficient and one column per
component), the list of flagged dimensions (`flagged`) and the `method`
used.

## Details

A dimension is flagged as problematic when its `sqrt(CN)` is greater
than or equal to `index_threshold` and at least two coefficients have a
variance proportion greater than or equal to `proportion_threshold` on
that dimension. The default `index_threshold = 10` corresponds to a `CN`
of 100 on the eigenvalue (non-square-rooted) scale, since each
component's `sqrt(CN)` is simply the square root of that component's
`CN`.

The variance-decomposition proportions depend on how the columns are
scaled, so they are only comparable with those of other software when
the same matrix is decomposed. For instance, with `method = "X"` and a
model with an intercept the proportions are those of
`VisCollin::colldiag(mod, add.intercept = TRUE)`, since both decompose
the unit-length (unweighted) design matrix. The thresholds above are
meant for the unit-length definitions (`"WS"`, `"MP"` and `"X"`); with
`method = "RAW"` the condition indices are dominated by the units of the
variables.

## References

Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). *Regression Diagnostics:
Identifying Influential Data and Sources of Collinearity*. Wiley.

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
y <- rgamma(n, shape = 5, rate = 5 / mu)
mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
bkw_diagnostics(mod)
#> Collinearity diagnostics (Belsley-Kuh-Welsch) on (X'WX)_ul [method = "WS"]
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

# Same decomposition on the original (unweighted) design matrix
bkw_diagnostics(mod, method = "X")
#> Collinearity diagnostics (Belsley-Kuh-Welsch) on X_ul'X_ul [method = "X"]
#> 
#>      Eigenvalue Sing.value sqrt(CN)
#> dim1      2.017      1.420    1.000
#> dim2      1.091      1.044    1.360
#> dim3      1.020      1.010    1.406
#> dim4      0.871      0.933    1.522
#> dim5      0.001      0.038   36.973
#> 
#> Variance-decomposition proportions (row = coefficient, column = component):
#>              dim1  dim2  dim3  dim4  dim5
#> (Intercept) 0.001 0.477 0.071 0.450 0.001
#> x1          0.001 0.000 0.000 0.000 0.999
#> x2          0.001 0.000 0.000 0.000 0.999
#> x3          0.004 0.002 0.782 0.212 0.001
#> x4          0.003 0.424 0.117 0.455 0.000
#> 
#> >> Possible collinearity problems (sqrt(CN) >= 10 and proportion >= 0.5 on >= 2 variables):
#>   - dim5 (sqrt(CN) = 36.97): x1, x2
```
