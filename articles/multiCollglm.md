# Introduction to multiCollglm

`multiCollglm` diagnoses collinearity in generalized linear models
(GLMs), generalizing the classical linear-regression diagnostics
(condition number, Belsley-Kuh-Welsch decomposition, and the Redefined
Variance Inflation Factor, RVIF (Salmerón et al. 2025)) to any family
and link function. All three are computed on the design matrix weighted
by the IRLS weights at convergence and, by default, scaled to unit
length by column (without centering).
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

| `method` | Transformation | Comments |
|----|----|----|
| `"RAW"` | None at all: eigenvalues of $`X'\hat{W}X`$ directly (intercept included, no centering, no scaling). The plain baseline against which the other four can be judged – usually dominated by however differently-scaled the explanatory variables happen to be. |  |
| `"MP"` | Mackinnon and Puterman (1989). The *original* explanatory variables (intercept included) are rescaled to unit Euclidean length with [`multiColl::lu()`](https://cran.r-project.org/package=multiColl) **before** fitting. | A GLM’s IRLS weights at convergence are invariant to any linear rescaling of `X` – but done explicitly to mirror the original definition. |
| `"MS"` | Marx and Smith (1990). The explanatory variables (excluding the intercept, if any) are **centered on their own mean and rescaled to unit Euclidean length before fitting**, with an ordinary intercept re-added if the original model had one. | Unlike `"MP"`, centering shifts the fitted values themselves, so the refit is genuinely different from `model`. |
| `"WS"` | Weissfeld and Sereika (1991), the package default. The GLM is fit on the original variables and the IRLS-weighted information matrix is scaled to unit column length *afterwards* (no centering). |  |
| `"OZ"` | Ozkale (2021). Same as `"WS"`, but the IRLS-weighted design matrix is additionally **centered** before being scaled to unit length. If `model` has an intercept, its column is **dropped** before centering/scaling (every worked example in Ozkale 2019, Sec. 5, fits without an intercept term). |  |

A few things worth knowing before picking one:

- **`"MP"` requires the `multiColl` package** (`Suggests`, not a hard
  dependency); it errors with an informative message if it isn’t
  installed.
- **`"MS"` can fail to converge for a model with no intercept.**
  Centering removes every column’s mean. For some family/link
  combinations (most notably Gamma with the inverse link, which needs
  strictly positive fitted values) the refit’s IRLS iteration can
  diverge. When that happens you get an informative error explaining why
  this is expected for that combination and is not a bug: Marx and
  Smith (1990) developed the method for logistic regression models,
  which ordinarily do include an intercept. If you hit this, use `"WS"`
  or `"OZ"` for that particular model instead.
- **Centering an IRLS-weighted design matrix that includes an intercept
  can leave it exactly rank-deficient** for *some* family/link
  combinations (for example the Gamma family with its canonical
  (inverse) link specifically – `sqrt(w) * eta` is exactly 1 for every
  observation, regardless of whether an intercept is present). `"OZ"`
  sidesteps this for most families by dropping the intercept before
  centering, but for Gamma-inverse the singularity is structural either
  way and `NC_OZ` will then legitimately come back as `Inf` with a
  warning. See the [`Nitrogen` reproduction
  article](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-nitrogen.html)
  for a worked example.
- Leave `method = NULL` (the default) to keep using `center`/`scale`
  directly, exactly as in earlier versions of this function – see the
  next section.

`"MP"` (Mackinnon and Puterman 1989) additionally requires the
`multiColl` package. See
[`?condition_number`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)
for the full argument reference.

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

`scale = TRUE` rescales every eigenvalue by the same factor
(`nrow(X) - 1`) relative to the default `scale = "unit"`, so `CN` and
`sqrt(CN)` are unaffected; `scale = FALSE` or a custom numeric vector
generally do change them, since the columns are no longer put on a
comparable scale.

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
#> sqrt(CN) of X (original design matrix, multiColl::CNs(X), with intercept): 36.9727
#> sqrt(CN) of W (MacKinnon-Puterman information matrix, method = "MP"): 85.2735
#> ratio (sqrt(CN)_W / sqrt(CN)_X): 2.3064
#> 
#> Collinearity in X (sqrt(CN)_X > 30): YES
#> ML-collinearity (ratio_wx > 5 and sqrt(CN)_W > 30): no
```

`kappa_x` (`sqrt(CN)` of `X`) is obtained via
[`multiColl::CNs()`](https://cran.r-project.org/package=multiColl):
`CNs(X)` returns a list with `CN1` (“Condition Number without
intercept”) and `CN2` (“Condition Number with intercept”); since
`kappa_X` standardizes `X` *including* the constant vector, `CNs(X)$CN2`
is the one used.

    ML-collinearity diagnostic (Lesaffre and Marx, 1993)

    sqrt(CN) of X (original design matrix, multiColl::CNs(X)$CN2): 190.7765
    sqrt(CN) of W (MacKinnon-Puterman information matrix, method = "MP"): 329.9542
    ratio (sqrt(CN)_W / sqrt(CN)_X): 1.7295

    Collinearity in X (sqrt(CN)_X > 30): YES
    ML-collinearity (ratio_wx > 5 and sqrt(CN)_W > 30): no

(the printed example above is Lesaffre and Marx’s own Lee
cancer-remission model, reproduced exactly in the [`LeeCancer`
reproduction
article](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-lee-cancer.html).)

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
exactly as the other three diagnostics do. **Requires the `multiColl`
package** (`Suggests`, not a hard dependency). See
[`?ml_collinearity`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
for the full argument reference.

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

Marx, B. D., and E. P. Smith. 1990. “Weighted Multicollinearity in
Logistic Regression: Diagnostics and Biased Estimation Techniques with
an Example from Lake Acidification.” *Canadian Journal of Fisheries and
Aquatic Sciences* 47 (6): 1128–35. <https://doi.org/10.1139/f90-131>.

Ozkale, M. R. 2021. “The Red Indicator and Corrected VIFs in Generalized
Linear Models.” *Communications in Statistics - Simulation and
Computation* 50 (12): 4144–70.
<https://doi.org/10.1080/03610918.2019.1639740>.

Salmerón, R., C. B. García, and J. García. 2025. “A Redefined Variance
Inflation Factor: Overcoming the Limitations of the Variance Inflation
Factor.” *Computational Economics* 65: 337–63.
<https://doi.org/10.1007/s10614-024-10575-8>.

Weissfeld, L. A., and S. M. Sereika. 1991. “A Multicollinearity
Diagnostic for Generalized Linear Models.” *Communications in
Statistics - Theory and Methods* 20 (4): 1183–98.
<https://doi.org/10.1080/03610929108830558>.
