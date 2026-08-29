# Reproduction: Chatterjee and Hadi (1988) -- Nitrogen data

## Context

The `Nitrogen` dataset records weather conditions and nitrogen dioxide
concentration (`y`, p.p.h.m.) measured in September 1984 at a monitoring
station in the San Francisco Bay Area, with four regressors: wind speed
(`x1`), maximum temperature (`x2`), insolation (`x3`) and a stability
factor (`x4`). The dataset was introduced by Chatterjee and Hadi (1988),
and has been reused by at least nine later works to illustrate
collinearity diagnostics in a Gamma GLM with the inverse link.

Chatterjee and Hadi describe the dataset as having 26 observations, but
the original observation numbering jumps from observation 8 to
observation 10 (there is no observation 9): there are in fact only 25
observations, which is what `data(Nitrogen)` includes in this package.
As far as we know, this discrepancy had not previously been flagged in
the literature that reuses this dataset.

## Data

``` r

library(multiCollglm)
data(Nitrogen)
str(Nitrogen)
#> 'data.frame':    25 obs. of  5 variables:
#>  $ y : int  6 5 5 3 7 9 6 2 10 7 ...
#>  $ x1: num  11.1 12.1 12 17.8 9.5 7.2 11.5 13.4 10.8 13.8 ...
#>  $ x2: int  90 86 80 70 90 100 92 74 87 78 ...
#>  $ x3: int  382 380 372 352 358 362 302 316 339 328 ...
#>  $ x4: int  12 20 19 16 10 12 15 15 14 14 ...

mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"),
           data = Nitrogen)
summary(mod)
#> 
#> Call:
#> glm(formula = y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"), 
#>     data = Nitrogen)
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)  
#> (Intercept)  2.992e-02  2.652e-01   0.113   0.9113  
#> x1           1.966e-02  7.798e-03   2.522   0.0203 *
#> x2          -1.024e-03  2.478e-03  -0.413   0.6837  
#> x3           5.271e-05  1.748e-04   0.301   0.7661  
#> x4          -1.056e-03  2.279e-03  -0.463   0.6481  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.1126254)
#> 
#>     Null deviance: 5.7610  on 24  degrees of freedom
#> Residual deviance: 2.2855  on 20  degrees of freedom
#> AIC: 110.57
#> 
#> Number of Fisher Scoring iterations: 5
```

## Original diagnostic (a tour of the literature)

The first work to use this dataset is Ozkale and Arican (2016)
\[`review3`\], which computes the condition number on $`X'\hat{W}X`$ as
the ratio between the largest and smallest eigenvalue (without a square
root), with a reciprocal link and an intercept, **standardizing the four
regressors before adding the column of ones**. They obtain $`CN =
213.8097`$, which they interpret as evidence of severe collinearity.
This value is then repeated without being recalculated in Amin et
al. (2020) \[`review2`\], Lukman et al. (2022) \[`review16`\], Akram et
al. (2022) \[`review34`\] and Shewa and Ugwuowo (2023) \[`review94`\] –
in the last two cases the fitted model actually uses a logarithmic link,
not the reciprocal link that originally produced the value, so the
condition number they cite almost certainly does not correspond to their
own weight matrix.

Kurtoglu (2021) \[`review97`\] performs an analogous calculation but
with the IRLS weights of their own log-gamma model, obtaining a much
larger value: $`49,649,929`$. Bulut (2023) \[`review47`\] applies the
same dataset under an Inverse Gaussian distribution (not Gamma) with a
log link, and obtains $`10,263.65`$.

Only two works document a full transformation of the data: Ozkale and
Abbasi (2022) \[`review66`\] weight by $`\hat{W}^{1/2}`$, center, and
rescale each column to unit length, obtaining $`CN = 2.49\times10^{8}`$
(and computing the VIF via a Moore-Penrose pseudoinverse because the
resulting matrix is nearly singular, which they themselves note is why
their VIF values come out artificially low, 0.19-5.49, despite the
astronomically large condition number). Ozkale (2021) \[`review42`\],
with the same weighting/centering/unit-length-scaling procedure, reports
both an untransformed scenario ($`CN = 1.4820\times10^{4}`$) and a
transformed one ($`CN = 1.5508\times10^{16}`$, driven by a numerically
zero eigenvalue).

We were unable to reproduce the value from Ozkale and Arican (2016)
($`\hat\varphi^2 = 0.07572852`$, $`CN = 213.8097`$): refitting the model
under several standardizations consistently yields
$`\hat\varphi^2 \approx
0.1126`$ (a quantity invariant to any rescaling of the design matrix),
which suggests some unreported detail in their estimation procedure
rather than a data or scaling issue.

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
| RAW    | NC_RAW   |     7.727554e+07 |       8790.6504 |
| MP     | NC_MP    |     5.762046e+03 |         75.9081 |
| MS     | NC_MS    |     2.433344e+02 |         15.5992 |
| WS     | NC_WS    |     6.080451e+03 |         77.9772 |
| OZ     | NC_OZ    |     5.791925e+03 |         76.1047 |

`method = "MP"` requires the `multiColl` package; if it is not
installed, that row shows `NA` instead of failing.

## Comparison

|  | Literature | multiCollglm |
|----|----|----|
| Most-cited transformation | Standardize before fitting, reciprocal link (`review3`) | `"WS"`: no centering, rescaled to unit length after the IRLS fit |
| “Standard” value (`CN`, no square root) | 213.8097, repeated without recalculation in 5 works | `NC_WS` = 6080.4514 |
| With full centering + scaling (`review42`, `review66`) | $`1.55\times10^{16}`$ / $`2.49\times10^{8}`$: one eigenvalue numerically zero | `NC_OZ` = 5791.9254 (centered, without the exact-rank problem because the intercept is dropped before centering) |
| Condition index (square root) with `"WS"` | – | $`\approx`$ 78, consistent with the order of magnitude most works that keep the intercept report |

The pattern is the same as in the other datasets covered by this
package: the “213.8097” figure that dominates this literature is, for
the most part, an inherited citation rather than an independent
calculation, and the two only works that do fully transform the data
(`review42`, `review66`) obtain condition numbers several orders of
magnitude larger because centering the matrix weighted by
$`\hat{W}^{1/2}`$ pushes one eigenvalue to floating-point precision – a
numerical instability rather than a stable measure of collinearity.

This dataset also illustrates the structural reason already documented
in
[`?condition_number`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md):
for the Gamma family with its canonical (inverse) link, `sqrt(w) * eta`
equals exactly 1 for every observation, so centering the IRLS-weighted
matrix (as `"OZ"` and `"MS"` do) can leave it exactly rank-deficient –
with or without an intercept. That is exactly what explains the
near-zero eigenvalues in `review42` and `review66`.

## Conclusion

Keeping the intercept and applying only the unit-length rescaling to the
already IRLS-weighted matrix (`"WS"`, this package’s default method)
gives a moderate but well-defined condition number: no eigenvalue
collapses to zero, so the diagnostic is numerically interpretable and
comparable across datasets or across values of a ridge parameter.
Centering before scaling (`"OZ"`, `"MS"`) is still valid here, because
the intercept is dropped before centering, but it’s worth remembering
that, for this dataset and this family/link, that centering is exactly
what triggers the astronomical condition numbers reported by two of the
nine reviewed works.

## References

- Chatterjee, S. and Hadi, A.S. (1988). *Sensitivity Analysis in Linear
  Regression*. John Wiley and Sons.
- Ozkale, M.R. and Arican, E. (2016). A new biased estimator in logistic
  regression model. *Statistics*, 50(2), 233-253.
  <https://doi.org/10.1080/02331888.2015.1123711>
- Amin, M., Qasim, M., Amanullah, M. and Afzal, S. (2020). Performance
  of some ridge estimators for the gamma regression model. *Statistical
  Papers*, 61(3), 997-1026. <https://doi.org/10.1007/s00362-017-0971-z>
- Lukman, A.F., Ayinde, K., Kibria, B.M.G. and Adewuyi, E.T. (2022).
  Modified ridge-type estimator for the gamma regression model.
  *Communications in Statistics - Simulation and Computation*, 51(9),
  5009-5023. <https://doi.org/10.1080/03610918.2020.1752720>
- Akram, M.N., Kibria, B.M.G., Abonazel, M.R. and Afzal, N. (2022). On
  the performance of some biased estimators in the gamma regression
  model: simulation and applications. *Journal of Statistical
  Computation and Simulation*, 92(12), 2425-2447.
  <https://doi.org/10.1080/00949655.2022.2032059>
- Shewa, G.A. and Ugwuowo, F.I. (2023). Kibria-Lukman type estimator for
  gamma regression model. *Concurrency and Computation: Practice and
  Experience*, 35(1). <https://doi.org/10.1002/cpe.7441>
- Kurtoglu, F. (2021). Modified ridge parameter estimators for log-gamma
  model: Monte Carlo evidence with a graphical investigation.
  *Communications in Statistics - Simulation and Computation*, 50(7),
  2168-2183. <https://doi.org/10.1080/03610918.2019.1650181>
- Bulut, Y.M. (2023). Inverse Gaussian Liu-type estimator.
  *Communications in Statistics - Simulation and Computation*, 52(10),
  4864-4879. <https://doi.org/10.1080/03610918.2021.1971243>
- Ozkale, M.R. and Abbasi, A. (2022). Iterative restricted OK estimator
  in generalized linear models and the selection of tuning parameters
  via MSE and genetic algorithm. *Statistical Papers*, 63(6), 1979-2040.
  <https://doi.org/10.1007/s00362-022-01304-0>
- Ozkale, M.R. (2021). The red indicator and corrected VIFs in
  generalized linear models. *Communications in Statistics - Simulation
  and Computation*, 50(12), 4144-4170.
  <https://doi.org/10.1080/03610918.2019.1639740>
