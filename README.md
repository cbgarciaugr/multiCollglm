# multiCollglm

Collinearity diagnostics for generalized linear models (GLMs), on `glm`
objects or on regularized models fitted with `glmnet` / `cv.glmnet`.

`condition_number()`, `bkw_diagnostics()` and `rvif_diagnostics()`
generalize the classical linear-regression collinearity diagnostics to
any family and link function. By default (`method = "WS"`, Weissfeld and
Sereika, 1991) they are computed on the design matrix weighted by the
IRLS weights at convergence and scaled to unit length by column, without
centering, so the intercept takes part in the diagnosis. The `method`
argument also gives the definition of MacKinnon and Puterman (`"MP"`),
the one of the original design matrix $X'X$, unweighted (`"X"`), and the
raw information matrix $X'WX$ (`"RAW"`). `ml_collinearity()` adds a fourth,
complementary diagnostic, comparing the condition number of the
*original* design matrix against that of the IRLS-weighted information
matrix (Lesaffre and Marx, 1993) to tell apart ordinary collinearity
among the explanatory variables from *ML-collinearity*.

## Installation

```r
install.packages("remotes")
remotes::install_github("cbgarciaugr/multiCollglm")
# or, from a local copy of the package folder:
# remotes::install_local("path/to/multiCollglm_0.4.0")
```

## Quick example

```r
library(multiCollglm)

mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = data)

condition_number(mod)                  # Weissfeld and Sereika (default)
condition_number(mod, method = "MP")   # MacKinnon and Puterman
condition_number(mod, method = "X")    # original design matrix X'X (unit length)
condition_number(mod, method = "RAW")  # X'WX, no transformation
bkw_diagnostics(mod)
rvif_diagnostics(mod)
ml_collinearity(mod)
```

## Learn more

`vignette("multiCollglm")` is the package's complete guide: every
function, the four condition-number definitions available via `method`,
the `glmnet`/`cv.glmnet` interface, and the full ML-collinearity
rationale.

The [package website](https://cbgarciaugr.github.io/multiCollglm/) also
has the [full function reference](https://cbgarciaugr.github.io/multiCollglm/reference/index.html)
an article [comparing the package with VisCollin](https://cbgarciaugr.github.io/multiCollglm/articles/comparison-viscollin.html)
and a series of [reproduction articles](https://cbgarciaugr.github.io/multiCollglm/articles/index.html):
real published cases reproduced with `multiCollglm`, using four example
datasets bundled with the package (`LeeCancer`, `Nitrogen`, `Mine`,
`Plastic`) and a fifth (`Bodyfat`, from the `TH.data` package).
