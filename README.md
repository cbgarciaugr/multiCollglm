# multiCollglm

Collinearity diagnostics for generalized linear models (GLMs), on `glm`
objects or on regularized models fitted with `glmnet` / `cv.glmnet`.

All three diagnostics are computed, generically for any family and link
function, on the design matrix weighted by the IRLS weights at
convergence (`mod$weights`) and scaled to unit length by column (without
centering). For `glmnet`/`cv.glmnet`, since there is no IRLS fit and no
convergence weights, the package identifies the active set of variables
at the given `lambda` and refits a standard `glm()` on them so the same
diagnostics can be computed.

See `vignette("multiCollglm")` for a getting-started guide, and the
[package website](https://<usuario-github>.github.io/multiCollglm/) for
worked reproductions of published works where collinearity was diagnosed
without the correct data transformation, compared against the results
from this package.

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

`condition_number()` reproduces exactly the calculation:

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
Salmeron, Garcia and Garcia (2025), which is able to detect both
essential and non-essential collinearity and decomposes the percentage
of near collinearity attributable to each variable.

## Changing the data transformation

All three functions accept `center` and `scale` arguments, applied to the
IRLS-weighted design matrix before computing any diagnostic. They accept
the same values as base R's own [`scale()`](https://rdrr.io/r/base/scale.html):

- `center`: `FALSE` (default), `TRUE` (subtract column means), or a
  numeric vector of per-column values to subtract.
- `scale`: `"unit"` (default; Euclidean/L2 unit-length columns, the
  transformation originally proposed), `TRUE` (root-mean-square scaling),
  `FALSE` (no scaling), or a numeric vector of per-column divisors.

```r
condition_number(mod, scale = TRUE) # root-mean-square scaling instead of unit length
condition_number(mod, scale = FALSE) # no scaling at all
```

`scale = TRUE` rescales every eigenvalue by the same factor
(`nrow(X) - 1`) relative to the default `scale = "unit"`, so the
condition number and condition index are unaffected; `scale = FALSE` or a
custom numeric vector generally do change them, since the columns are no
longer put on a comparable scale.

**Do not use `center = TRUE`** with an intercept model unless you know
what you are doing: for a canonical-link GLM (e.g. Gamma with the inverse
link), the IRLS-weighted design matrix satisfies `Xw %*% coef(mod) == 1`
exactly, so centering an intercept column is mathematically guaranteed to
drop the matrix's rank by exactly one. `bkw_diagnostics()` and
`rvif_diagnostics()` detect this and raise an informative error rather
than returning `NaN`; this is precisely why Belsley, Kuh and Welsch
recommend against centering for collinearity diagnostics.

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
4. Apply the same calculation as in the direct `glm` case to that model.

The result also includes `active_vars` and `dropped_vars` so you know
exactly on which sub-model the diagnostic was computed.

## BKW rule of thumb

`bkw_diagnostics()` flags as suspicious any dimension whose condition
index is `>= 10` and on which at least two coefficients have a variance
proportion `>= 0.5` (Belsley, Kuh and Welsch, 1980). The condition index
is on the singular-value (square-rooted) scale, so this threshold of 10
corresponds to a condition number of 100 on the eigenvalue scale used by
`condition_number()` (`condition_index = sqrt(condition_number)`). Both
thresholds are adjustable via `index_threshold` and
`proportion_threshold`.

## References

- Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). *Regression
  Diagnostics: Identifying Influential Data and Sources of
  Collinearity*. Wiley.
- Salmeron, R., Garcia, C.B. and Garcia, J. (2025). A redefined Variance
  Inflation Factor: overcoming the limitations of the Variance Inflation
  Factor. *Computational Economics*, 65, 337-363.
  <https://doi.org/10.1007/s10614-024-10575-8>
