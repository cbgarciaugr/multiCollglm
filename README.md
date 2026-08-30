# multiCollglm

Collinearity diagnostics for generalized linear models (GLMs), on `glm`
objects or on regularized models fitted with `glmnet` / `cv.glmnet`.

All three diagnostics (`condition_number()`, `bkw_diagnostics()` and 
`rvif_diagnostics()`) are computed, generically for any family and link
function, on the design matrix weighted by the IRLS weights at
convergence (`mod$weights`). By default that matrix is scaled to unit
length by column and **not** centered — this is the Weissfeld and Sereika
(1991) definition, described in detail below. `condition_number()` also
lets you pick one of five published condition-number definitions
directly via its `method` argument, or build a custom transformation
yourself with `center`/`scale`.

For `glmnet`/`cv.glmnet`, since there is no IRLS fit and no convergence
weights, the package identifies the active set of variables at the given
`lambda` and refits a standard `glm()` on them so the same diagnostics
can be computed.

See `vignette("multiCollglm")` for a getting-started guide, and the
[package website](https://cbgarciaugr.github.io/multiCollglm/articles/index.html) for
worked reproductions of published works, using four example datasets
bundled with the package (`LeeCancer`, `Nitrogen`, `Mine`,`Plastic`) and a fifth (`Bodyfat` included in `TH.data` package). 

## Installation

```r
# from the local package folder
install.packages("remotes")
remotes::install_local("C:/Users/Usuario/Documents/R/multiCollglm")
```

## Usage with `glm`

```r
library(multiCollglm)

mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = data)

condition_number(mod)
bkw_diagnostics(mod)
rvif_diagnostics(mod)
```

`condition_number(mod)`, with no `method` argument, reproduces exactly
the calculation proposed by Weissfeld and Sereika (1991):

```r
X <- model.matrix(mod)
W <- diag(mod$weights)
Wsqrt <- diag(sqrt(diag(W)))
Xw <- Wsqrt %*% X
X_unit <- scale(Xw, center = FALSE, scale = sqrt(colSums(Xw^2)))
M_star <- t(X_unit) %*% X_unit
eig_star <- sort(eigen(M_star, symmetric = TRUE)$values, decreasing = TRUE)
CN <- max(eig_star) / min(eig_star)
sqrt(CN)
```

`rvif_diagnostics()` calls [`rvif::rvifs()`](https://cran.r-project.org/package=rvif)
on that same `X_unit` matrix (with `ul = TRUE`, so it is re-scaled to
unit length right before the RVIF is computed, a no-op under the default
transformation) to get the Redefined Variance Inflation Factor of
@salmeron2025, which is able to detect both
essential and non-essential collinearity and decomposes the percentage
of near collinearity attributable to each variable.

## Condition-number methods (`method` argument)

`condition_number()` accepts a `method` argument that is a shortcut for
four published  definitions of the
CN for a GLM (`nc_label`: `NC_MP`, `NC_MS`, `NC_WS` or `NC_OZ`). It is also included a method that does not include transformation (`NC_RAW`). Passing `method` overrides whatever
`center`/`scale` you supplied and labels the result accordingly:

```r
condition_number(mod, method = "RAW") # no transformation at all
condition_number(mod, method = "MP")  # Mackinnon and Puterman (1989)
condition_number(mod, method = "MS")  # Marx and Smith (1990)
condition_number(mod, method = "WS")  # Weissfeld and Sereika (1991) -- the default
condition_number(mod, method = "OZ")  # Ozkale (2019)
```

| `method` | Transformation | Comments|
|---|---|---|
| `"RAW"` | None at all: eigenvalues of $X'\hat{W}X$ directly (intercept included, no centering, no scaling). The plain baseline against which the other four can be judged — usually dominated by however differently-scaled the explanatory variables happen to be. |  |
| `"MP"` | @mackinnon1989. The *original* explanatory variables (intercept included) are rescaled to unit Euclidean length with [`multiColl::lu()`](https://cran.r-project.org/package=multiColl) **before** fitting. | A GLM's IRLS weights at convergence are invariant to any linear rescaling of `X` — but done explicitly to mirror the original definition). |
| `"MS"` | @marx1990. The explanatory variables (excluding the intercept, if any) are **centered on their own mean and rescaled to unit Euclidean length before fitting**, with an ordinary intercept re-added if the original model had one. | Unlike `"MP"` centering shifts the fitted values themselves, so the refit is genuinely different from `model`. |
| `"WS"` | @weissfeld1991, the package default. The GLM is fit on the original variables and the IRLS-weighted information matrix is scaled to unit column length *afterwards* (no centering). |  |
| `"OZ"` | @ozkale2021. Same as `"WS"`, but the IRLS-weighted design matrix is additionally **centered** before being scaled to unit length. If `model` has an intercept, its column is **dropped** before centering/scaling (every worked example in Ozkale 2019, Sec. 5, fits without an intercept term). |  |

A few things worth knowing before picking one:

- **`"MP"` requires the `multiColl` package** (`Suggests`, not a hard
  dependency); it errors with an informative message if it isn't
  installed.
- **`"MS"` can fail to converge for a model with no intercept.**
  Centering removes every column's mean. For some family/link combinations (most notably Gamma
  with the inverse link, which needs strictly positive fitted values)
  the refit's IRLS iteration can diverge. When that happens you get an
  informative error explaining why this is expected for that combination and is not a bug: @marx1990
  developed the method for logistic regression models, which
  ordinarily do include an intercept. If you hit this, use `"WS"` or
  `"OZ"` for that particular model instead.
- **Centering an IRLS-weighted design matrix that includes an intercept can leave it exactly rank-deficient** for *some* family/link
  combinations (for example the Gamma family with its canonical (inverse) link
  specifically (`sqrt(w) * eta` is exactly 1 for every observation,
  regardless of whether an intercept is present). `"OZ"` sidesteps this
  for most families by dropping the intercept before centering, but for
  Gamma-inverse the singularity is structural either way and `NC_OZ` will
  then legitimately come back as `Inf` with a warning. See [`Nitrogen` reproduction article]() on the
  package website for a worked example.
- Leave `method = NULL` (the default) to keep using `center`/`scale`
  directly, exactly as in earlier versions of this function — see the
  next section.

## Changing the data transformation directly (without `method`)

All three functions also accept `center` and `scale` arguments, applied
to the IRLS-weighted design matrix before computing any diagnostic (these
are ignored whenever `method` is supplied). They accept the same values
as base R's own [`scale()`](https://rdrr.io/r/base/scale.html):

- `center`: `FALSE` (default), `TRUE` (subtract column means), or a
  numeric vector of per-column values to subtract.
- `scale`: `TRUE` (root-mean-square scaling),
  `FALSE` (no scaling), or a numeric vector of per-column divisors.

In addition, `scale` accept also `"unit"`, for unit legth transformation:
```r
condition_number(mod, scale="unit") # unit length transformation
condition_number(mod, scale = TRUE) # root-mean-square scaling instead of unit length
condition_number(mod, scale = FALSE) # no scaling at all
```
scale = TRUE rescales every eigenvalue by the same factor (nrow(X) - 1) relative to the default scale = "unit", 
so CN and sqrt(CN) are unaffected; scale = FALSE or a custom numeric vector generally do change them, 
since the columns are no longer put on a comparable scale.

**Do not use `center = TRUE`** with an intercept model unless you know
what you are doing, for the same reason described above for `"OZ"`:
`bkw_diagnostics()` and `rvif_diagnostics()` detect the resulting
rank-deficiency and raise an informative error rather than returning
`NaN`; this is precisely why Belsley, Kuh and Welsch recommend against
centering for collinearity diagnostics. If you specifically want a
centered diagnostic, prefer `method = "OZ"` or `method = "MS"`, which
handle the intercept column correctly instead of centering it directly.

## Usage with `glmnet` / `cv.glmnet`

```r
library(glmnet)
library(multiCollglm)

x <- as.matrix(data[, c("x1", "x2", "x3", "x4")])
y <- data$y

fit <- glmnet(x, y, family = Gamma(link = "inverse"))
cvfit <- cv.glmnet(x, y, family = Gamma(link = "inverse"))

condition_number(fit, x = x, y = y, family = Gamma(link = "inverse"), s = 0.01)
condition_number(cvfit, x = x, y = y, family = Gamma(link = "inverse"), s = "lambda.min")

bkw_diagnostics(cvfit, x = x, y = y, family = Gamma(link = "inverse"), s = "lambda.min")
rvif_diagnostics(cvfit, x = x, y = y, family = Gamma(link = "inverse"), s = "lambda.min")
```

In these cases, `condition_number()` / `bkw_diagnostics()` /
`rvif_diagnostics()`:

1. Take the coefficient vector at the given `lambda` (`coef(fit, s = s)`).
2. Keep the variables with a coefficient different from zero (the active
   set; the intercept does not count for this).
3. Refit `glm(y ~ active_variables, family = family, weights = weights)`
   to obtain `mod$weights` (real IRLS weights at convergence).
4. Apply the same calculation as in the direct `glm` case to that model
   (including, if you pass one, the `method` shortcut described above).

The result also includes `active_vars` and `dropped_vars` so you know
exactly on which sub-model the diagnostic was computed.

## BKW rule of thumb

`bkw_diagnostics()` flags as suspicious any dimension whose `sqrt(CN)`
is `>= 10` and on which at least two coefficients have a variance
proportion `>= 0.5` (@belsley1980). `sqrt(CN)`
is on the singular-value (square-rooted) scale, so this threshold of 10
corresponds to a `CN` of 100 on the eigenvalue scale used by
`condition_number()` (`sqrt(CN)` is simply the square root of `CN`, i.e.
`condition_index = sqrt(condition_number)`). Both
thresholds are adjustable via `index_threshold` and
`proportion_threshold`.

## Example datasets

Four real datasets used throughout the package's reproduction articles
ship with the package (`data(Nitrogen)`, `data(Mine)`, `data(LeeCancer)`,
`data(Plastic)`. A fifth, `bodyfat`, is used  but not
bundled here since it already ships in the
[`TH.data`](https://cran.r-project.org/package=TH.data) package.

## References


