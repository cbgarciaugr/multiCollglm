# Reproduction: cancer remission data (Lee 1974)

## Context

Lee’s (1974) cancer remission dataset records clinical observations from
27 patients, of whom 9 had complete remission (`y = 1`) and 18 had
incomplete remission (`y = 0`), together with six quantitative
covariates: cellularity index (CELL), smear index (SMEAR), infiltrate
index (INFIL), labeling index (LI), percentage of blast cells (BLAST)
and body temperature (TEMP). It is a classic example for collinearity
diagnostics in binomial logistic regression, used with different subsets
of covariates by Lesaffre and Marx (1993) \[`review133`\], Ozkale and
Arican (2016) \[`review35`\] together with Ozkale (2021) \[`review42`\],
and Huang, Jou and Cho (2015) \[`review137`\].

Lesaffre and Marx (1993) are the first to use this dataset, with only
three covariates (LI, TEMP, CELL) that show severe collinearity – a
classic scenario of extreme variance inflation and numerical
instability. The remaining reviewed works use all covariates except, in
some cases, LI.

## Data

``` r

library(multiCollglm)
data(LeeCancer)
str(LeeCancer)
#> 'data.frame':    27 obs. of  7 variables:
#>  $ y : int  1 1 0 0 1 0 1 0 0 0 ...
#>  $ x1: num  0.8 0.9 0.8 1 0.9 1 0.95 0.95 1 0.95 ...
#>  $ x2: num  0.83 0.36 0.88 0.87 0.75 0.65 0.97 0.87 0.45 0.36 ...
#>  $ x3: num  0.66 0.32 0.7 0.87 0.68 0.65 0.92 0.83 0.45 0.34 ...
#>  $ x4: num  1.9 1.4 0.8 0.7 1.3 0.6 1 1.9 0.8 0.5 ...
#>  $ x5: num  1.1 0.74 0.176 1.053 0.519 ...
#>  $ x6: num  0.996 0.992 0.982 0.986 0.98 ...
# TEMP is recorded in the original source on a /100 scale; it is
# multiplied back by 100 to match the scale used by the cited works.
LeeCancer$x6100 <- LeeCancer$x6 * 100
```

## Model 1: Lesaffre and Marx (1993) – CELL, LI, TEMP

``` r

mod133 <- glm(y ~ x1 + x4 + x6100, family = binomial(), data = LeeCancer)
summary(mod133)
#> 
#> Call:
#> glm(formula = y ~ x1 + x4 + x6100, family = binomial(), data = LeeCancer)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  67.6339    56.8875   1.189   0.2345  
#> x1            9.6522     7.7511   1.245   0.2130  
#> x4            3.8671     1.7783   2.175   0.0297 *
#> x6100        -0.8207     0.6171  -1.330   0.1835  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 34.372  on 26  degrees of freedom
#> Residual deviance: 21.953  on 23  degrees of freedom
#> AIC: 29.953
#> 
#> Number of Fisher Scoring iterations: 7
```

## Model 2: Ozkale and Arican (2016) / Ozkale (2021) – CELL, SMEAR, INFIL, BLAST, TEMP (no LI)

``` r

mod3542 <- glm(y ~ x1 + x2 + x3 + x5 + x6100, family = binomial(),
                data = LeeCancer)
summary(mod3542)
#> 
#> Call:
#> glm(formula = y ~ x1 + x2 + x3 + x5 + x6100, family = binomial(), 
#>     data = LeeCancer)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  56.9903    56.2314   1.013   0.3108  
#> x1           23.4507    31.1197   0.754   0.4511  
#> x2           27.6185    38.6541   0.715   0.4749  
#> x3          -31.6856    41.5318  -0.763   0.4455  
#> x5            2.9396     1.6240   1.810   0.0703 .
#> x6100        -0.8091     0.5367  -1.508   0.1317  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 34.372  on 26  degrees of freedom
#> Residual deviance: 26.095  on 21  degrees of freedom
#> AIC: 38.095
#> 
#> Number of Fisher Scoring iterations: 7
```

Ozkale and Arican (2016) standardize and center the explanatory
variables before fitting, explicitly stating that *“standardization
removes any non-essential ill-conditioning resulting from the
intercept”* – that is, they compute a **classical condition number,
without any IRLS weighting**, on the standardized data. They obtain
$`CN = 295.7026`$ (eigenvalue ratio, no square root), which on the
square-root scale would be $`17.196`$.

Ozkale (2021), with the same model, reports the untransformed
eigenvalues $`15.6562, 0.4601, 0.1557, 0.0145, 0.0001`$ and from those a
condition number of $`1.3032\times10^{5}`$ – a value that is obtained
exactly neither as the ratio ($`15.6562/0.0001 = 156,562`$) nor as its
square root ($`\approx 395.7`$) from those very eigenvalues she reports.
The same happens with her centered-and-scaled scenario: eigenvalues
$`4.4133,
0.3657, 0.2009, 0.0200, 0.0002`$, from which the ratio would be
$`22,066.5`$, not the $`2.6099\times10^{4}`$ she reports.

This same model (5 covariates, no LI, no intercept) was the subject of
an additional check earlier in this project: refitting the model from
four different starting points (least squares, zeros, sign-reversed
least squares, and [`glm()`](https://rdrr.io/r/stats/glm.html)’s own
defaults), all four converge to the same unique point
$`(38.9501, 46.4911, -50.9001, 2.2396, -38.3061)`$ – which, since
logistic regression is a convex problem, confirms it is the global
maximum-likelihood estimate for these data. That point does **not**
match the coefficients published by Ozkale (2019, Sec. 5.2) for the same
model ($`-45.5605, 46.2072, -47.9364, 3.5485, -48.5129`$), which
suggests an error in Ozkale’s own paper rather than in the data
(independently verified against SAS documentation and against the
original SAS `datalines` for the `Remission` dataset).

## Model 3: Huang, Jou and Cho (2015) – all covariates

``` r

mod137 <- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6100, family = binomial(),
              data = LeeCancer)
summary(mod137)
#> 
#> Call:
#> glm(formula = y ~ x1 + x2 + x3 + x4 + x5 + x6100, family = binomial(), 
#>     data = LeeCancer)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  58.0385    71.2364   0.815   0.4152  
#> x1           24.6615    47.8377   0.516   0.6062  
#> x2           19.2936    57.9500   0.333   0.7392  
#> x3          -19.6013    61.6815  -0.318   0.7507  
#> x4            3.8960     2.3371   1.667   0.0955 .
#> x5            0.1511     2.2786   0.066   0.9471  
#> x6100        -0.8743     0.6757  -1.294   0.1957  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 34.372  on 26  degrees of freedom
#> Residual deviance: 21.751  on 20  degrees of freedom
#> AIC: 35.751
#> 
#> Number of Fisher Scoring iterations: 8
```

Huang, Jou and Cho (2015) use all six covariates and report $`CN =
418.13`$.

## Diagnostic with multiCollglm

``` r

methods <- c("RAW", "MP", "MS", "WS", "OZ")

calc_table <- function(mod) {
  do.call(rbind, lapply(methods, function(m) {
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
}

cat("Model 1 (Lesaffre and Marx), with intercept:\n")
#> Model 1 (Lesaffre and Marx), with intercept:
knitr::kable(calc_table(mod133), digits = 4)
```

| method | nc_label | condition_number | condition_index |
|:-------|:---------|-----------------:|----------------:|
| RAW    | NC_RAW   |     1.168157e+08 |      10808.1330 |
| MP     | NC_MP    |     1.088698e+05 |        329.9542 |
| MS     | NC_MS    |     2.567380e+02 |         16.0230 |
| WS     | NC_WS    |     9.931110e+04 |        315.1366 |
| OZ     | NC_OZ    |     2.788063e+02 |         16.6975 |

``` r


cat("Model 2 (Ozkale), with intercept:\n")
#> Model 2 (Ozkale), with intercept:
knitr::kable(calc_table(mod3542), digits = 4)
```

| method | nc_label | condition_number | condition_index |
|:-------|:---------|-----------------:|----------------:|
| RAW    | NC_RAW   |    219838806\.31 |      14826.9621 |
| MP     | NC_MP    |        163327.20 |        404.1376 |
| MS     | NC_MS    |         23296.31 |        152.6313 |
| WS     | NC_WS    |        143044.11 |        378.2117 |
| OZ     | NC_OZ    |         14371.56 |        119.8814 |

``` r

mod3542b <- glm(y ~ x1 + x2 + x3 + x5 + x6100 - 1, family = binomial(),
                 data = LeeCancer)
cat("Model 2 (Ozkale), without intercept (their own specification):\n")
#> Model 2 (Ozkale), without intercept (their own specification):
knitr::kable(calc_table(mod3542b), digits = 4)
```

| method | nc_label | condition_number | condition_index |
|:-------|:---------|-----------------:|----------------:|
| RAW    | NC_RAW   |     1.942571e+08 |      13937.6155 |
| MP     | NC_MP    |     7.769953e+04 |        278.7464 |
| MS     | NC_MS    |     4.202298e+02 |         20.4995 |
| WS     | NC_WS    |     7.476822e+04 |        273.4378 |
| OZ     | NC_OZ    |     1.643415e+04 |        128.1957 |

``` r


cat("Model 3 (Huang, Jou and Cho), with intercept:\n")
#> Model 3 (Huang, Jou and Cho), with intercept:
knitr::kable(calc_table(mod137), digits = 4)
```

| method | nc_label | condition_number | condition_index |
|:-------|:---------|-----------------:|----------------:|
| RAW    | NC_RAW   |    395524494\.29 |      19887.7976 |
| MP     | NC_MP    |        248614.00 |        498.6121 |
| MS     | NC_MS    |         43277.14 |        208.0316 |
| WS     | NC_WS    |        218108.47 |        467.0208 |
| OZ     | NC_OZ    |         33948.94 |        184.2524 |

``` r

mod137b <- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6100 - 1, family = binomial(),
               data = LeeCancer)
cat("Model 3 (Huang, Jou and Cho), without intercept:\n")
#> Model 3 (Huang, Jou and Cho), without intercept:
knitr::kable(calc_table(mod137b), digits = 4)
```

| method | nc_label | condition_number | condition_index |
|:-------|:---------|-----------------:|----------------:|
| RAW    | NC_RAW   |     2.900940e+08 |      17032.1463 |
| MP     | NC_MP    |     1.431963e+05 |        378.4129 |
| MS     | NC_MS    |     5.158097e+02 |         22.7114 |
| WS     | NC_WS    |     1.393080e+05 |        373.2398 |
| OZ     | NC_OZ    |     3.372989e+04 |        183.6570 |

## Comparison

For Model 2 (Ozkale), the classical, unweighted condition number
reported by Ozkale and Arican (2016), $`CN = 295.7026`$, does not match
any of
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)’s
five methods, because these always weight by the IRLS weights at
convergence – exactly the transformation that, per this package’s whole
premise, should be applied and that most reviewed works omit.

For Model 3 (Huang, Jou and Cho), the reported value ($`418.13`$) is
closest, in order of magnitude, to the condition index (square-root
scale) of `"MP"`/`"WS"` **without** an intercept (373.24 and 378.41
respectively) than to any of the with-intercept variants, but it does
not match exactly.

|  | Literature | multiCollglm |
|----|----|----|
| Model 2, stated transformation | Centered and scaled, **not** IRLS-weighted | `"OZ"` (centered+scaled, IRLS-weighted, no intercept): 128.2 |
| Model 3, transformation | Not specified | `"WS"` without intercept: 373.24 (closest to the published 418.13) |

## Conclusion

The three published models fitted to this dataset each illustrate, in
their own way, why the transformation matters: Lesaffre and Marx (1993)
work with only three already-known-to-be-collinear covariates; Ozkale
and Arican (2016) compute a classical diagnostic without IRLS weighting
that, by construction, cannot match any of this package’s methods (all
of which weight by the IRLS weights); and the number Ozkale (2021)
reports from her own eigenvalues shows an internal arithmetic
inconsistency we could not resolve. In addition, for the 5-covariate
no-intercept model, the fit Ozkale (2019) reports does not match the
unique global maximum-likelihood estimate for these data (verified by
convergence from four different starting points), which points to an
error in her paper rather than in the data, already independently
verified against SAS’s official documentation.

## References

- Lee, E.T. (1974). A computer program for linear logistic regression
  analysis. *Computer Programs in Biomedicine*, 4(2), 80-92.
  <https://doi.org/10.1016/0010-468X(74)90011-7>
- Lesaffre, E. and Marx, B.D. (1993). Collinearity in Generalized Linear
  Regression. *Communications in Statistics - Theory and Methods*,
  22(7), 1933-1952. <https://doi.org/10.1080/03610929308831126>
- Ozkale, M.R. and Arican, E. (2016). A new biased estimator in logistic
  regression model. *Statistics*, 50(2), 233-253.
  <https://doi.org/10.1080/02331888.2015.1123711>
- Ozkale, M.R. (2021). The red indicator and corrected VIFs in
  generalized linear models. *Communications in Statistics - Simulation
  and Computation*, 50(12), 4144-4170.
  <https://doi.org/10.1080/03610918.2019.1639740>
- Huang, C.-C.L., Jou, Y.-J. and Cho, H.-J. (2015). A New
  Multicollinearity Diagnostic for Generalized Linear Models. *Journal
  of Applied Statistics*.
  <https://doi.org/10.1080/02664763.2015.1126239>
