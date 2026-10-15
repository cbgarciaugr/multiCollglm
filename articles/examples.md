# Examples: real-data applications

Every example below applies the functions of `multiCollglm` to a data
set that was used in the literature to illustrate collinearity in a GLM,
and compares the result with what was published. They follow the same
structure: context, data, diagnostic with `multiCollglm`, comparison
with the literature. The code of each article is run when the website is
built, so all outputs can be reproduced.

## Index

| Example | Data | Family and link | What it shows |
|:---|:---|:---|:---|
| [Lee cancer remission](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-lee-cancer.md) | `LeeCancer` | Binomial, logit | Reproduces the values of Lesaffre and Marx (1993) for the condition indices of the design matrix (190.78) and of the MacKinnon-Puterman matrix (329.95), and the ML-collinearity diagnostic. |
| [Nitrogen dioxide](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-nitrogen.md) | `Nitrogen` | Gamma, inverse | Why centering the weighted matrix makes it singular for the canonical Gamma link. |
| [Coal-mine fractures](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-mines.md) | `Mine` | Poisson, log | Reproduces the published eigenvalues of the information matrix (condition index 2093.32). |
| [Body fat](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-bodyfat.md) | `bodyfat` (`TH.data`) | Gamma, log | IRLS weights identically 1: WS, MP and X coincide. |
| [Plastic injection](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-plastic.md) | `Plastic` | Poisson, log | Widely cited values that cannot be reproduced; large effect of the weights. |
| [Comparison with VisCollin](https://cbgarciaugr.github.io/multiCollglm/articles/comparison-viscollin.md) | the five above | several | [`colldiag()`](https://rdrr.io/pkg/VisCollin/man/colldiag.html) equals `method = "X"`; where the IRLS weights make them differ; factors, weights and plots. |

## All data sets at a glance

The largest condition index (square-root scale) of each model for the
four definitions available through the `method` argument:

``` r

library(multiCollglm)

lee <- LeeCancer
lee$x6100 <- lee$x6 * 100

models <- list(
  "Lee cancer (binomial, logit)"  = glm(y ~ x1 + x4 + x6100, binomial(), lee),
  "Nitrogen (Gamma, inverse)"     = glm(y ~ x1 + x2 + x3 + x4, Gamma("inverse"), Nitrogen),
  "Mine (Poisson, log)"           = glm(y ~ x1 + x2 + x3 + x4, poisson(), Mine),
  "Plastic (Poisson, log)"        = glm(Y ~ X1 + X2 + X3 + X4, poisson(), Plastic)
)
if (requireNamespace("TH.data", quietly = TRUE)) {
  data("bodyfat", package = "TH.data")
  models[["Body fat (Gamma, log)"]] <- glm(
    DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4,
    family = Gamma("log"), data = bodyfat)
}

tab <- t(sapply(models, function(m)
  sapply(c(WS = "WS", MP = "MP", X = "X", RAW = "RAW"),
         function(k) condition_number(m, method = k)$condition_index)))
knitr::kable(tab, digits = 2)
```

|                              |     WS |     MP |      X |      RAW |
|:-----------------------------|-------:|-------:|-------:|---------:|
| Lee cancer (binomial, logit) | 315.14 | 329.95 | 190.78 | 10808.13 |
| Nitrogen (Gamma, inverse)    |  77.98 |  75.91 |  67.71 |  8790.65 |
| Mine (Poisson, log)          |  28.58 |  26.45 |  19.90 |  2093.32 |
| Plastic (Poisson, log)       | 876.84 | 899.04 | 595.74 | 73561.62 |
| Body fat (Gamma, log)        | 352.04 | 352.04 | 352.04 |  3410.63 |

These are the values of Table 1 of the manuscript on multicollinearity
diagnosis in GLMs (with intercept). `WS` is the default of the package.
`X` ignores the IRLS weights (it is what
[`VisCollin::colldiag()`](https://rdrr.io/pkg/VisCollin/man/colldiag.html)
computes), `MP` is the definition of MacKinnon and Puterman, and `RAW`
is the information matrix without any transformation.

## Where to go next

- [Complete
  guide](https://cbgarciaugr.github.io/multiCollglm/articles/multiCollglm.md)
  to the functions and the `method` argument.
- [Function
  reference](https://cbgarciaugr.github.io/multiCollglm/reference/index.md).
- To add your own example, copy `reproduction-template.Rmd` (see the
  `README.md` in `vignettes/articles`) and list it in `_pkgdown.yml`.
