# Introduction to multiCollglm

`multiCollglm` diagnoses collinearity in generalized linear models
(GLMs), generalizing the classical linear-regression diagnostics
(condition number, Belsley-Kuh-Welsch decomposition, and the Redefined
Variance Inflation Factor, RVIF (Salmerón et al. 2025)) to any family
and link function. All three are computed on the design matrix weighted
by the IRLS weights at convergence and, by default, scaled to unit
length by column (without centering).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("cbgarciaugr/multiCollglm")
```

## Usage with `glm`

``` r

library(multiCollglm)

set.seed(1)
n <- 200
x1 <- rnorm(n)
x2 <- x1 + rnorm(n, sd = 0.05) # almost collinear with x1
x3 <- rnorm(n)
x4 <- rnorm(n)
mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3 + 0.1 * x4)
y <- rgamma(n, shape = 5, rate = 5 / mu)

mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"))
```

``` r

condition_number(mod)
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 3.2659 0.8814 0.6208 0.2315 0.0003
#> 
#> CN (eigenvalue scale): 10123.1410
#> sqrt(CN) (classical, singular-value scale): 100.6138
```

``` r

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

``` r

rvif_diagnostics(mod)
#> Redefined Variance Inflation Factor (RVIF)
#> 
#>                 RVIF      %
#> (Intercept)    3.288 69.589
#> x1          1554.971 99.936
#> x2          1546.036 99.935
#> x3             1.581 36.758
#> x4             1.111  9.960
```

## Condition-number methods

[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)
also accepts a `method` argument that switches between five published
(or, for `"RAW"`, textbook-baseline) definitions of CN for a GLM:

``` r

condition_number(mod, method = "RAW") # no transformation at all
#> Eigenvalues of X'WX (decreasing order):
#> [1] 25684.1167  4325.2941  3318.8684  1000.0895     3.3983
#> 
#> CN (NC_RAW) (eigenvalue scale): 7557.9499
#> sqrt(CN) (NC_RAW) (classical, singular-value scale): 86.9365
condition_number(mod, method = "MP")  # requires multiColl package
#> Eigenvalues of X_ul'WX_ul (decreasing order):
#> [1] 143.5803  18.7623  14.5861   5.0021   0.0197
#> 
#> CN (NC_MP) (eigenvalue scale): 7271.5771
#> sqrt(CN) (NC_MP) (classical, singular-value scale): 85.2735
condition_number(mod, method = "WS")  # same as the default above, explicitly labeled
#> Eigenvalues of (X'WX)_ul (decreasing order):
#> [1] 3.2659 0.8814 0.6208 0.2315 0.0003
#> 
#> CN (NC_WS) (eigenvalue scale): 10123.1410
#> sqrt(CN) (NC_WS) (classical, singular-value scale): 100.6138
condition_number(mod, method = "OZ")  # centered and unit-scaled [@ozkale2021]
#> Eigenvalues of (X'WX)_ulc (decreasing order):
#> [1] 2.4718 0.8806 0.6473 0.0003
#> 
#> CN (NC_OZ) (eigenvalue scale): 7201.4286
#> sqrt(CN) (NC_OZ) (classical, singular-value scale): 84.8612
condition_number(mod, method = "MS")  # centered, unit-scaled, then refit [@marx1990]
#> Eigenvalues of X_ulc'WX_ulc (decreasing order):
#> [1] 3816.0817   42.5254   18.0071   13.6069    0.0198
#> 
#> CN (NC_MS) (eigenvalue scale): 192963.4257
#> sqrt(CN) (NC_MS) (classical, singular-value scale): 439.2760
```

`"MP"` (Mackinnon and Puterman 1989) additionally requires the
`multiColl` package. See
[`?condition_number`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)
and the package’s
[README](https://github.com/cbgarciaugr/multiCollglm#condition-number-methods-method-argument)
for the full description of each method, and for the caveats around
`"MS"` (which can fail to converge for a no-intercept model) and `"OZ"`
(which can legitimately return `Inf` for a Gamma model with the
canonical inverse link).

## Changing the data transformation

All three functions accept `center` and `scale` arguments, applied to
the IRLS-weighted design matrix before computing any diagnostic (ignored
whenever `method` is supplied). They accept the same values as base R’s
own [`scale()`](https://rdrr.io/r/base/scale.html):

- `center`: `FALSE` (default), `TRUE` (subtract each column’s mean), or
  a numeric vector of per-column values to subtract.
- `scale`: `"unit"` (default; Euclidean/L2 unit-length columns, the
  originally proposed transformation), `TRUE` (root-mean-square
  scaling), `FALSE` (no scaling), or a numeric vector of per-column
  divisors.

``` r

condition_number(mod, scale = TRUE) # root-mean-square scaling instead of unit length
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 649.9211 175.4071 123.5469  46.0607   0.0642
#> 
#> CN (eigenvalue scale): 10123.1410
#> sqrt(CN) (classical, singular-value scale): 100.6138
```

**Do not use `center = TRUE`** on an intercept model unless you know
what you are doing: for a canonical-link GLM, the IRLS-weighted design
matrix satisfies `Xw %*% coef(mod) == 1` exactly, so centering the
intercept column is guaranteed to drop the matrix’s rank by exactly one.
This is precisely why Belsley, Kuh and Welsch (1980) recommend against
centering for collinearity diagnostics;
[`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md)
and
[`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md)
detect this and raise an informative error instead of returning `NaN`.

``` r

bkw_diagnostics(mod, center = TRUE)
#> Error:
#> ! The (weighted, transformed) design matrix is exactly rank-deficient, so the variance-decomposition proportions are undefined. Two known causes: (1) centering (center = TRUE) an IRLS-weighted design matrix that includes an intercept, for some family/link combinations (this is not a general canonical-link fact -- e.g. binomial with the logit link is unaffected); and (2) the Gamma family with its canonical (inverse) link, for which sqrt(w) * linear.predictors is exactly 1 for every observation regardless of whether an intercept is present, so removing the intercept does not avoid it there. This is exactly why Belsley, Kuh and Welsch recommend against centering for collinearity diagnostics; use the default center = FALSE, or remove the exact collinearity from the data.
```

## Usage with `glmnet` / `cv.glmnet`

When the model is a regularized fit (`glmnet` or `cv.glmnet`), there is
no IRLS fit and no convergence weights of its own. `multiCollglm`
identifies the active set of variables at the given `lambda` and refits
a standard [`glm()`](https://rdrr.io/r/stats/glm.html) on them so the
same diagnostics can be computed:

``` r

library(glmnet)

x <- as.matrix(data.frame(x1, x2, x3, x4))
cvfit <- cv.glmnet(x, y, family = Gamma(link = "inverse"))

condition_number(cvfit, x = x, y = y, family = Gamma(link = "inverse"), s = "lambda.min")
bkw_diagnostics(cvfit, x = x, y = y, family = Gamma(link = "inverse"), s = "lambda.min")
rvif_diagnostics(cvfit, x = x, y = y, family = Gamma(link = "inverse"), s = "lambda.min")
```

The result also includes `active_vars` and `dropped_vars`, so you know
exactly which sub-model the diagnostic was computed on.

## BKW rule of thumb

[`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md)
flags as suspicious any dimension whose `sqrt(CN)` is `>= 10` and on
which at least two coefficients have a variance proportion `>= 0.5`
(Belsley et al. 1980). Both thresholds are adjustable via
`index_threshold` and `proportion_threshold`.

## ML-collinearity

[`ml_collinearity()`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
implements the diagnostic of Lesaffre and Marx (1993), which compares
`sqrt(CN)` of the *original* design matrix (via
[`multiColl::CNs()`](https://rdrr.io/pkg/multiColl/man/CNs.html))
against `sqrt(CN)` of the IRLS-weighted information matrix under
`method = "MP"`, to distinguish ordinary collinearity among the
explanatory variables from *ML-collinearity*: ill-conditioning that
instead arises from the combination of the response, the fitted
coefficients and the link function.

``` r

ml_collinearity(mod)
#> ML-collinearity diagnostic (Lesaffre and Marx, 1993)
#> 
#> sqrt(CN) of X (original design matrix, multiColl::CNs(X)[1]): 36.9284
#> sqrt(CN) of W (MacKinnon-Puterman information matrix, method = "MP"): 85.2735
#> ratio (sqrt(CN)_W / sqrt(CN)_X): 2.3092
#> 
#> Collinearity in X (sqrt(CN)_X > 30): YES
#> ML-collinearity (ratio_wx > 5 and sqrt(CN)_W > 30): no
```

A model is flagged with `ml_collinearity = TRUE` only when both the
ratio and `sqrt(CN)` of the information matrix exceed their thresholds
(5 and 30 by default); see
[`?ml_collinearity`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
and the package’s
[README](https://github.com/cbgarciaugr/multiCollglm#ml-collinearity-lesaffre-and-marx-1993)
for the full rationale, and the [`LeeCancer` reproduction
article](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-lee-cancer.html)
for Lesaffre and Marx’s own worked example, reproduced exactly.

## Beyond this guide

This vignette covers the package’s basic usage. The [package
website](https://cbgarciaugr.github.io/multiCollglm/articles/index.html)
also includes a series of **reproduction articles**: real cases from the
literature where collinearity was diagnosed without applying the correct
data transformation, together with a comparison against the results
`multiCollglm` produces.

## References

Belsley, D. A., E. Kuh, and R. E. Welsch. 1980. *Regression Diagnostics:
Identifying Influential Data and Sources of Collinearity*. Wiley.

Lesaffre, E., and B. D. Marx. 1993. “Collinearity in Generalized Linear
Regression.” *Communications in Statistics - Theory and Methods* 22 (7):
1933–52. <https://doi.org/10.1080/03610929308831126>.

Mackinnon, M. J., and M. L. Puterman. 1989. “Collinearity in Generalized
Linear Models.” *Communications in Statistics - Theory and Methods* 18
(9): 3463–72. <https://doi.org/10.1080/03610928908830102>.

Salmerón, R., C. B. García, and J. García. 2025. “A Redefined Variance
Inflation Factor: Overcoming the Limitations of the Variance Inflation
Factor.” *Computational Economics* 65: 337–63.
<https://doi.org/10.1007/s10614-024-10575-8>.
