# Reproduction: coal mine roof-fracture data (Mine)

## Context

The `Mine` dataset records 44 observations of roof fractures in coal
mines in the Appalachian region of West Virginia, with four continuous
predictors: inner burden thickness, percentage extraction of the
previously mined lower seam, lower seam height, and time since the mine
opened. It is the real-data application shared by five reviewed works,
each fitting a different count-data GLM: Kurtoglu and Ozkale (2019)
\[`review55`\] and Kurtoglu and Ozkale (2019) \[`review118`\] (same
authors, a restricted ridge and a restricted Liu estimator respectively,
both with Poisson and a log link); Sami, Amin and Butt (2022)
\[`review8`\] (Conway-Maxwell-Poisson, a “COMPRR” ridge estimator);
Hussein and Algamal (2026) \[`review115`\] (Bell regression, a
Kibria-Lukman-type estimator); and Sami et al. (2026) \[`review117`\]
(discrete Weibull, a “DWRE” ridge estimator).

`review8` states that it analyzes a different distribution but provides
no diagnostic measures of its own (which should differ across
distributions), instead citing a previous paper to claim that the
dataset exhibits collinearity. `review117` does the same, citing
`review8` despite also using a different distribution; moreover, based
on the code in `review117`’s appendix, the diagnostic appears to be
computed on the unweighted $`X'X`$, and although the code includes a
standardization step, it is never actually applied.

Note on the data: two transcription errors present in some circulating
copies of this dataset (row 30: `y = 22` instead of `y = 2`; row 36:
`x1 = 0` instead of `x1 = 80`) have been corrected in `data(Mine)` after
cross-checking against the appendix data table of a second paper.

## Data

``` r

library(multiCollglm)
data(Mine)
str(Mine)
#> 'data.frame':    44 obs. of  5 variables:
#>  $ y : int  2 1 0 4 1 2 0 0 4 4 ...
#>  $ x1: int  50 230 125 75 70 65 65 350 350 160 ...
#>  $ x2: int  70 65 70 65 65 70 60 60 90 80 ...
#>  $ x3: int  52 42 45 68 53 46 62 54 54 38 ...
#>  $ x4: num  1 6 1 0.5 0.5 3 1 0.5 0.5 0 ...

mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(link = "log"),
           data = Mine)
summary(mod)
#> 
#> Call:
#> glm(formula = y ~ x1 + x2 + x3 + x4, family = poisson(link = "log"), 
#>     data = Mine)
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -3.5930896  1.0256803  -3.503  0.00046 ***
#> x1          -0.0014066  0.0008358  -1.683  0.09240 .  
#> x2           0.0623458  0.0122862   5.074 3.89e-07 ***
#> x3          -0.0020803  0.0050661  -0.411  0.68134    
#> x4          -0.0308135  0.0162648  -1.894  0.05816 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 74.984  on 43  degrees of freedom
#> Residual deviance: 37.856  on 39  degrees of freedom
#> AIC: 144.13
#> 
#> Number of Fisher Scoring iterations: 5
```

## Original diagnostic (a tour of the literature)

Kurtoglu and Ozkale (2019), fitting Poisson with a log link, obtain a
condition number of **2093.32** that at first could not be reproduced
with any of the package’s centered or scaled transformations. Hussein
and Algamal (2026), with a Bell distribution, obtain **2139.89** (same
order of magnitude, consistent with a different but Poisson-like count
distribution for these data).

## Diagnostic with multiCollglm

``` r

methods <- c("RAW", "MP", "MS", "WS", "OZ")

tbl <- do.call(rbind, lapply(methods, function(m) {
  res <- tryCatch(condition_number(mod, method = m),
                   error = function(e) NULL)
  if (is.null(res)) {
    data.frame(method = m, nc_label = NA, condition_number = NA,
               condition_index = NA)
  } else {
    data.frame(method = m, nc_label = res$nc_label,
               condition_number = res$condition_number,
               condition_index = res$condition_index)
  }
}))
knitr::kable(tbl, digits = 4)
```

| method | nc_label | condition_number | condition_index |
|:-------|:---------|-----------------:|----------------:|
| RAW    | NC_RAW   |     4381970.9542 |       2093.3158 |
| MP     | NC_MP    |         699.8016 |         26.4538 |
| MS     | NC_MS    |         101.5642 |         10.0779 |
| WS     | NC_WS    |         816.9398 |         28.5822 |
| OZ     | NC_OZ    |           4.3874 |          2.0946 |

## Comparison

The condition index (square-root scale) of `method = "RAW"` – no
transformation at all, directly on $`X'\hat{W}X`$ – is **2093.3158**,
which matches, up to rounding, the **2093.32** reported by Kurtoglu and
Ozkale (2019) \[`review55`, `review118`\]. In other words: the value the
text initially describes as “not reproducible” **is exactly reproduced**
once the `"RAW"` method is added to the package – what those two works
compute is the condition index (square root of the eigenvalue ratio) of
the IRLS-weighted information matrix, with no centering or scaling at
all.

When the data are rescaled to unit length without centering (`"WS"`,
which does include the intercept), the condition index drops to
**28.58** – of the same order as the reference value (~28.8) obtained
with this transformation; the conclusion of severe collinearity holds,
but with a much more moderate and stable value, no longer dominated by
the very different scales of the four variables (thickness, percentage,
height, years).

We could not independently verify Hussein and Algamal’s (2026) value for
the Bell regression, since that family is not part of R base’s
[`glm()`](https://rdrr.io/r/stats/glm.html) families and would require
an additional package; it is left as reported by the original source.

|  | Literature | multiCollglm |
|----|----|----|
| Transformation | None (Poisson, log link, with intercept) | `"RAW"`: untransformed |
| Value (condition index, square root) | 2093.32 (Kurtoglu and Ozkale 2019) | 2093.3158 |
| Rescaled to unit length, no centering | – | 28.5822 |
| Conclusion about collinearity | Severe | Severe with `"RAW"`; still relevant but less extreme with `"WS"` |

## Conclusion

This dataset is a good example of why it’s worth including
`method = "RAW"` as an explicit baseline: the most-cited reference value
for the mine data turns out to be exactly the condition index with no
transformation whatsoever – not an error or an undocumented calculation,
but simply the plainest possible case, which wasn’t among this package’s
published methods until it was added. Compared with `"WS"` (rescaled to
unit length, no centering), the diagnostic’s magnitude changes
substantially because the four regressors are on very different units
and `"RAW"` is dominated by that scale disparity, not only by the
underlying collinearity.

## References

- Kurtoglu, F. and Ozkale, M.R. (2019). Restricted ridge estimator in
  generalized linear models: Monte Carlo simulation studies on Poisson
  and binomial distributed responses. *Communications in Statistics -
  Simulation and Computation*, 48(4), 1191-1218.
  <https://doi.org/10.1080/03610918.2017.1408822>
- Kurtoglu, F. and Ozkale, M.R. (2019). Restricted Liu estimator in
  generalized linear models: Monte Carlo simulation studies on gamma and
  Poisson distributed responses. *Hacettepe Journal of Mathematics and
  Statistics*, 48(4), 1250-1276.
  <https://doi.org/10.15672/HJMS.2018.624>
- Sami, F., Amin, M. and Butt, M.M. (2022). On the ridge estimation of
  the Conway-Maxwell Poisson regression model with multicollinearity:
  Methods and applications. *Concurrency and Computation: Practice and
  Experience*, 34(1), e6477. <https://doi.org/10.1002/cpe.6477>
- Hussein, S.M.J.A. and Algamal, Z.Y. (2026). A New Shrinkage Estimator
  for the Bell Regression Model. *Thailand Statistician*, 24(2),
  462-473.
- Sami, F., Aljohani, H.M., Iqbal, K., Ghorashi, S. and Abu Bakar
  Siddique, M. (2026). A new estimator for discrete Weibull regression
  model with correlated variables: an illustrative example. *Journal of
  Statistical Computation and Simulation*, 96(11), 2737-2757.
  <https://doi.org/10.1080/00949655.2026.2636232>
- Myers, R.H. (1990). *Classical and Modern Regression with
  Applications* (2nd ed.). Duxbury Press.
