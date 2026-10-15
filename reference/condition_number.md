# Condition number of a model (collinearity diagnostic)

Computes the condition number of a generalized linear model, following
one of four definitions chosen with `method`. By default
(`method = "WS"`) it follows Weissfeld and Sereika (1991): the model is
fitted on the original design matrix and the information matrix
\\X'\hat{W}X\\ at convergence is rescaled so that its diagonal is all
ones, which is equivalent to scaling the IRLS-weighted design matrix
\\\hat{W}^{1/2}X\\ to unit length by column, *without* centering. This
generalizes to any GLM family and link function the classical condition
number of linear regression. The other definitions are the one of
MacKinnon and Puterman (1989), the one of the original (unweighted)
design matrix \\X'X\\, and the one of the untransformed information
matrix \\X'\hat{W}X\\.

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
  method = c("WS", "MP", "X", "RAW"),
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

- method:

  One of `"WS"` (default), `"MP"`, `"X"` or `"RAW"`. See Details.

## Value

An object of class `multicollglm_cn` with the eigenvalues of the matrix
defined by `method` (`eigenvalues`), the condition number as the ratio
between the largest and the smallest eigenvalue (`condition_number`),
its square root on the classical singular-value scale
(`condition_index`), the variable names (`var_names`), the method used
(`method`) and its label (`nc_label`, e.g. `"NC_WS"`). If the matrix is
exactly rank-deficient (perfect collinearity), the condition number is
`Inf` and a warning is issued.

## Details

### Methods

Let \\X\\ be the design matrix (including the intercept column, if any),
\\\hat{W}\\ the diagonal matrix of IRLS weights at convergence, and let
the subscript \\ul\\ denote scaling every column to unit Euclidean
length. None of the four methods centers the data, so the intercept
always takes part in the diagnostic and non-essential collinearity
(between the intercept and the explanatory variables) can be detected.

|  |  |  |
|----|----|----|
| `method` | Matrix whose eigenvalues are used | Reference |
| `"WS"` (default) | \\(X'\hat{W}X)\_{ul} = S\\X'\hat{W}X\\S\\, with \\S = \mathrm{diag}(X'\hat{W}X)^{-1/2}\\ | Weissfeld and Sereika (1991) |
| `"MP"` | \\X\_{ul}'\hat{W}X\_{ul}\\, with \\X\\ scaled to unit length *before* refitting the model | MacKinnon and Puterman (1989) |
| `"X"` | \\X\_{ul}'X\_{ul}\\: the original design matrix, no weights | Belsley, Kuh and Welsch (1980); \\\kappa_X\\ of Lesaffre and Marx (1993) |
| `"RAW"` | \\X'\hat{W}X\\: no transformation at all | – |

- `"WS"` fits the GLM on the raw `X` and rescales afterwards. Under
  orthogonality in the GLM metric (\\X'\hat{W}X\\ diagonal) it attains
  its minimum value, 1, for any weights, exactly as the unit-length
  condition number does in the linear model. This is the recommended
  definition.

- `"MP"` scales `X` to unit length first, refits the GLM on that matrix
  (this is done internally) and reads the information matrix off the
  refit. Because the weights enter after the scaling, the result is 1
  under orthogonality only if all the weights are equal.

- `"X"` ignores the model's weights, so it measures the collinearity of
  the explanatory variables themselves, not of the information matrix.
  It is what a linear-model diagnostic computes on the same predictors,
  e.g. `VisCollin::colldiag(mod, add.intercept = TRUE)`.

- `"RAW"` is the baseline without any transformation. It depends on the
  units of the explanatory variables and is usually far larger than the
  others, even when the underlying collinearity is the same.

The condition number is returned on two scales, as in the literature:
`condition_number`, the ratio \\\lambda\_{max}/\lambda\_{min}\\ of the
eigenvalues, and `condition_index`, its square root, which is the
classical condition number computed from the singular values of the
matrix. Thresholds such as 10 (weak), 30 (moderate to strong) or 100
(Belsley, Kuh and Welsch, 1980) are stated on the square-root scale.

### Removed in version 0.4.0

`method = "OZ"` (Ozkale, 2021) and `method = "MS"` (Marx and Smith,
1990), as well as the arguments `center` and `scale`, no longer exist.
All of them centered the data, which removes the intercept from the
analysis and so rules out the detection of non-essential collinearity.

## References

Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). *Regression Diagnostics:
Identifying Influential Data and Sources of Collinearity*. Wiley.

Lesaffre, E. and Marx, B.D. (1993). Collinearity in generalized linear
regression. *Communications in Statistics - Theory and Methods*, 22(7),
1933-1952.
[doi:10.1080/03610929308831126](https://doi.org/10.1080/03610929308831126)

MacKinnon, M.J. and Puterman, M.L. (1989). Collinearity in generalized
linear models. *Communications in Statistics - Theory and Methods*,
18(9), 3463-3472.
[doi:10.1080/03610928908830102](https://doi.org/10.1080/03610928908830102)

Weissfeld, L.A. and Sereika, S.M. (1991). A multicollinearity diagnostic
for generalized linear models. *Communications in Statistics - Theory
and Methods*, 20(4), 1183-1198.
[doi:10.1080/03610929108830558](https://doi.org/10.1080/03610929108830558)

## See also

[`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md),
[`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md),
[`ml_collinearity()`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md).

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n); x2 <- x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n); x4 <- rnorm(n)
mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
y <- rgamma(n, shape = 5, rate = 5 / mu)
mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))

condition_number(mod) # Weissfeld and Sereika (1991), the default
#> Eigenvalues of (X'WX)_ul (decreasing order):
#> [1] 3.2659 0.8814 0.6208 0.2315 0.0003
#> 
#> CN (NC_WS, eigenvalue scale): 10123.1410
#> sqrt(CN) (NC_WS, classical, singular-value scale): 100.6138

# The four definitions side by side (sqrt scale)
sapply(c("WS", "MP", "X", "RAW"),
       function(m) condition_number(mod, method = m)$condition_index)
#>        WS        MP         X       RAW 
#> 100.61382  85.27354  36.97273  86.93647 
```
