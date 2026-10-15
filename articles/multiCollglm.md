# Introduction to multiCollglm

`multiCollglm` diagnoses collinearity in generalized linear models
(GLMs), generalizing the classical linear-regression diagnostics
(condition number, Belsley-Kuh-Welsch decomposition, and the Redefined
Variance Inflation Factor, RVIF (Salmerón et al. 2025)) to any family
and link function. By default all three follow Weissfeld and Sereika
(1991): they are computed on the design matrix weighted by the IRLS
weights at convergence and scaled to unit length by column, **without
centering**, so the intercept takes part in the diagnosis and
non-essential collinearity (between the intercept and the explanatory
variables) can be detected. The `method` argument gives access to three
more definitions, including the one of the original design matrix
\\X'X\\.
[`ml_collinearity()`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
adds a fourth, complementary diagnostic: Lesaffre and Marx’s (1993)
comparison between the condition number of the *original* design matrix
and that of the IRLS-weighted information matrix, used to tell apart
ordinary collinearity among the explanatory variables from
*ML-collinearity* (see below). This vignette is the package’s complete
guide; see
[`vignette("multiCollglm")`](https://cbgarciaugr.github.io/multiCollglm/articles/multiCollglm.md)
any time you need to come back to it.

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
#> Eigenvalues of (X'WX)_ul (decreasing order):
#> [1] 3.2659 0.8814 0.6208 0.2315 0.0003
#> 
#> CN (NC_WS, eigenvalue scale): 10123.1410
#> sqrt(CN) (NC_WS, classical, singular-value scale): 100.6138
```

``` r

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
```

``` r

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

## Condition-number methods

[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md),
[`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md)
and
[`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md)
all take a `method` argument with four possible values. The default is
`"WS"`:

``` r

condition_number(mod, method = "WS")  # the default, written out
#> Eigenvalues of (X'WX)_ul (decreasing order):
#> [1] 3.2659 0.8814 0.6208 0.2315 0.0003
#> 
#> CN (NC_WS, eigenvalue scale): 10123.1410
#> sqrt(CN) (NC_WS, classical, singular-value scale): 100.6138
condition_number(mod, method = "MP")
#> Eigenvalues of X_ul'WX_ul (decreasing order):
#> [1] 143.5803  18.7623  14.5861   5.0021   0.0197
#> 
#> CN (NC_MP, eigenvalue scale): 7271.5771
#> sqrt(CN) (NC_MP, classical, singular-value scale): 85.2735
condition_number(mod, method = "X")   # original design matrix, no weights
#> Eigenvalues of X_ul'X_ul (decreasing order):
#> [1] 2.0166 1.0907 1.0201 0.8710 0.0015
#> 
#> CN (NC_X, eigenvalue scale): 1366.9826
#> sqrt(CN) (NC_X, classical, singular-value scale): 36.9727
condition_number(mod, method = "RAW") # X'WX with no transformation
#> Eigenvalues of X'WX (decreasing order):
#> [1] 25684.1167  4325.2941  3318.8684  1000.0895     3.3983
#> 
#> CN (NC_RAW, eigenvalue scale): 7557.9499
#> sqrt(CN) (NC_RAW, classical, singular-value scale): 86.9365
```

Let \\X\\ be the design matrix (intercept included), \\\hat W\\ the IRLS
weights at convergence, and let the subscript \\ul\\ mean “every column
scaled to unit Euclidean length”:

| `method` | Matrix whose eigenvalues are used | Comments |
|----|----|----|
| `"WS"` (default) | \\(X'\hat W X)\_{ul} = S\\X'\hat W X\\S\\, with \\S=\mathrm{diag}(X'\hat W X)^{-1/2}\\ (Weissfeld and Sereika 1991) | The GLM is fitted on the original variables and the information matrix is rescaled afterwards (equivalently: unit-length columns of \\\hat W^{1/2}X\\). Equals 1 whenever \\X'\hat W X\\ is diagonal. |
| `"MP"` | \\X\_{ul}'\hat W X\_{ul}\\, with \\X\\ scaled **before** refitting (Mackinnon and Puterman 1989) | The weights enter after the scaling, so it equals 1 under orthogonality only if all the weights are equal. |
| `"X"` | \\X\_{ul}'X\_{ul}\\ | The original design matrix, ignoring the weights: the classical linear-model diagnostic of the same predictors (\\\kappa_X\\ of Lesaffre and Marx (1993)). |
| `"RAW"` | \\X'\hat W X\\ | No transformation at all. Dominated by the units of the variables; a baseline only. |

None of the four centers the data. This is deliberate: centering removes
the intercept column from the analysis, so a near-linear relationship
between the intercept and the regressors (non-essential collinearity)
can no longer be seen (Belsley et al. 1980), and for some family/link
combinations it is not even possible – for the Gamma family with its
canonical (inverse) link, \\\sqrt{\hat w_i}\\\hat\eta_i = 1\\ for every
observation, so a centered, weighted design matrix is exactly singular.
The [`Nitrogen` reproduction
article](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-nitrogen.html)
shows it on real data.

Some things worth knowing:

- All of them are on the classical scale: `condition_index` is
  \\\sqrt{\lambda\_{max}/\lambda\_{min}}\\ (the one that the thresholds
  10, 30 and 100 refer to) and `condition_number` is the ratio
  \\\lambda\_{max}/\lambda\_{min}\\ itself.
- With a Gamma family and *log* link the IRLS weights are identically 1,
  so `"WS"`, `"MP"` and `"X"` coincide (this happens in the `Bodyfat`
  application). The same happens for a Gaussian model with identity
  link.
- [`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md)
  calls [`rvif::rvifs()`](https://rdrr.io/pkg/rvif/man/rvifs.html) with
  `ul = TRUE`, which rescales the matrix to unit length once more. Hence
  `"WS"`, `"MP"` and `"RAW"` give the *same* RVIF; only `"X"` (no
  weights) differs.
- `bkw_diagnostics(mod, method = "X")` is, for a model with intercept,
  exactly `VisCollin::colldiag(mod, add.intercept = TRUE)`. See the
  [comparison with
  VisCollin](https://cbgarciaugr.github.io/multiCollglm/articles/comparison-viscollin.html).

### Changes with respect to version 0.3.0

Version 0.4.0 removes `method = "OZ"` (Ozkale), `method = "MS"` (Marx
and Smith) and the arguments `center` and `scale`, because all of them
centered the data. Using them now stops with an explanatory message. The
`multiColl` package is no longer needed for `method = "MP"` or for
[`ml_collinearity()`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md).

``` r

condition_number(mod, method = "OZ")
#> Error:
#> ! method = "OZ" was removed in multiCollglm 0.4.0. Available methods: "WS", "MP", "X", "RAW".
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
`sqrt(CN)` of the *original* design matrix,
`condition_number(mod, method = "X")`, against `sqrt(CN)` of the
IRLS-weighted information matrix under `method = "MP"`, to distinguish
ordinary collinearity among the explanatory variables from
*ML-collinearity*: ill-conditioning that instead arises from the
combination of the response, the fitted coefficients and the link
function. Both are computed natively; no other package is needed.

``` r

ml_collinearity(mod)
#> ML-collinearity diagnostic (Lesaffre and Marx, 1993)
#> 
#> sqrt(CN) of X (original design matrix, unit length, with intercept; method = "X"): 36.9727
#> sqrt(CN) of W (MacKinnon-Puterman information matrix, method = "MP"): 85.2735
#> ratio (sqrt(CN)_W / sqrt(CN)_X): 2.3064
#> 
#> Collinearity in X (sqrt(CN)_X > 30): YES
#> ML-collinearity (ratio_wx > 5 and sqrt(CN)_W > 30): no
```

The example below is Lesaffre and Marx’s own Lee cancer-remission model
(`TEMP` is expressed in hundredths of a degree in the packaged data,
hence the factor 100), whose published \\\kappa_X = 190.78\\ is
reproduced exactly (see also the [`LeeCancer` reproduction
article](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-lee-cancer.html)):

``` r

lee <- LeeCancer
lee$x6100 <- lee$x6 * 100
mod_lee <- glm(y ~ x1 + x4 + x6100, family = binomial(), data = lee)
ml_collinearity(mod_lee)
#> ML-collinearity diagnostic (Lesaffre and Marx, 1993)
#> 
#> sqrt(CN) of X (original design matrix, unit length, with intercept; method = "X"): 190.7765
#> sqrt(CN) of W (MacKinnon-Puterman information matrix, method = "MP"): 329.9542
#> ratio (sqrt(CN)_W / sqrt(CN)_X): 1.7295
#> 
#> Collinearity in X (sqrt(CN)_X > 30): YES
#> ML-collinearity (ratio_wx > 5 and sqrt(CN)_W > 30): no
```

Following Lesaffre and Marx (1993, Sec. 4.2), a model is flagged with
`ml_collinearity = TRUE` only when **both** the ratio
`sqrt(CN)_W / sqrt(CN)_X` exceeds `ratio_threshold` (default 5) **and**
`sqrt(CN)_W` itself exceeds `kappa_threshold` (default 30) – a high
ratio by itself, with a small `sqrt(CN)_W`, does not indicate an
ill-conditioned information matrix to begin with. Separately,
`sqrt(CN)_X > kappa_threshold` flags ordinary collinearity among the
explanatory variables (`collinearity_x`); if both flags are set and the
ratio is close to 1, both types of collinearity are simultaneously
present.
[`ml_collinearity()`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
also has `glmnet`/`cv.glmnet` methods, refitting on the active set
exactly as the other three diagnostics do. See
[`?ml_collinearity`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
for the full argument reference.

## Beyond this guide

This vignette covers the package’s basic usage. The [package
website](https://cbgarciaugr.github.io/multiCollglm/articles/index.html)
also includes a series of **reproduction articles**: real cases from the
literature where collinearity was diagnosed without applying the correct
data transformation, together with a comparison against the results
`multiCollglm` produces, and an article comparing the package with
**VisCollin**, the reference R implementation of the linear-model
collinearity diagnostics and their visualization.

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

Weissfeld, L. A., and S. M. Sereika. 1991. “A Multicollinearity
Diagnostic for Generalized Linear Models.” *Communications in
Statistics - Theory and Methods* 20 (4): 1183–98.
<https://doi.org/10.1080/03610929108830558>.
