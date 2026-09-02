# multiCollglm

Collinearity diagnostics for generalized linear models (GLMs), on `glm`
objects or on regularized models fitted with `glmnet` / `cv.glmnet`.

`condition_number()`, `bkw_diagnostics()` and `rvif_diagnostics()`
generalize the classical linear-regression collinearity diagnostics to
any family and link function, computed on the design matrix weighted by
the IRLS weights at convergence. `ml_collinearity()` adds a fourth,
complementary diagnostic, comparing the condition number of the
*original* design matrix against that of the IRLS-weighted information
matrix (Lesaffre and Marx, 1993) to tell apart ordinary collinearity
among the explanatory variables from *ML-collinearity*.

## Installation

```r
# from the local package folder
install.packages("remotes")
remotes::install_local("C:/Users/Usuario/Documents/R/multiCollglm")
```

## Quick example

```r
library(multiCollglm)

mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), data = data)

condition_number(mod)
bkw_diagnostics(mod)
rvif_diagnostics(mod)
ml_collinearity(mod)
```

## Learn more

`vignette("multiCollglm")` is the package's complete guide: every
function, the five condition-number definitions available via `method`,
the `glmnet`/`cv.glmnet` interface, and the full ML-collinearity
rationale.

The [package website](https://cbgarciaugr.github.io/multiCollglm/) also
has the [full function reference](https://cbgarciaugr.github.io/multiCollglm/reference/index.html)
and a series of [reproduction articles](https://cbgarciaugr.github.io/multiCollglm/articles/index.html):
real published cases reproduced with `multiCollglm`, using four example
datasets bundled with the package (`LeeCancer`, `Nitrogen`, `Mine`,
`Plastic`) and a fifth (`Bodyfat`, from the `TH.data` package).
