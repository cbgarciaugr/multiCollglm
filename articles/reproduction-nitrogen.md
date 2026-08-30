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
```

## Diagnostic with multiCollglm and comparison with literature

The first work (within our selection) to use this dataset is Ozkale and
Arican (2016), which computes $`CN`$ on $`X'\hat{W}X`$ as the ratio
between the largest and smallest eigenvalue (without a square root),
with a reciprocal link and an intercept, **standardizing the four
regressors before adding the column of ones**. They obtain $`CN =
213.8097`$, which they interpret as evidence of severe collinearity.
This value is then repeated without being recalculated in Amin et al.
(2020), Lukman et al. (2022), Akram et al. (2022) and Shewa and Ugwuowo
(2023). In the last two cases the fitted model actually uses a
logarithmic link, not the reciprocal link that originally produced the
value, so the $`CN`$ they cite almost certainly does not correspond to
their own weight matrix.

Kurtoglu (2021) performs an analogous calculation but with the IRLS
weights of their own log-gamma model, obtaining a much larger value:
$`49,649,929`$. Bulut (2023) applies the same dataset under an Inverse
Gaussian distribution (not Gamma) with a log link, and obtains
$`10,263.65`$.

Only two works document a full transformation of the data: Ozkale and
Abbasi (2022) weight by $`\hat{W}^{1/2}`$, center, and rescale each
column to unit length, obtaining $`CN = 2.49\times10^{8}`$ (and
computing the VIF via a Moore-Penrose pseudoinverse because the
resulting matrix is nearly singular, which they themselves note is why
their VIF values come out artificially low, 0.19-5.49, despite the
astronomically large $`CN`$). Ozkale (2021), with the same
weighting/centering/unit-length-scaling procedure, reports both an
untransformed scenario ($`CN = 1.4820\times10^{4}`$) and a transformed
one ($`CN = 1.5508\times10^{16}`$, driven by a numerically zero
eigenvalue).

We were unable to reproduce the value from Ozkale and Arican (2016)
($`\hat\varphi^2 = 0.07572852`$, $`CN = 213.8097`$): refitting the model
under several standardizations consistently yields
$`\hat\varphi^2 \approx
0.1126`$ (a quantity invariant to any rescaling of the design matrix),
which suggests some unreported detail in their estimation procedure
rather than a data or scaling issue.

``` r

knitr::kable(tbl, digits = 4)
```

| method | nc_label |           CN |  sqrt(CN) |
|:-------|:---------|-------------:|----------:|
| RAW    | NC_RAW   | 7.727554e+07 | 8790.6504 |
| MP     | NC_MP    | 5.762046e+03 |   75.9081 |
| MS     | NC_MS    | 2.433344e+02 |   15.5992 |
| WS     | NC_WS    | 6.080451e+03 |   77.9772 |
| OZ     | NC_OZ    | 5.791925e+03 |   76.1047 |

This dataset also illustrates the structural reason already documented
in
[`?condition_number`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md):
for the Gamma family with its canonical (inverse) link, `sqrt(w) * eta`
equals exactly 1 for every observation, so centering the IRLS-weighted
matrix (as `"OZ"` and `"MS"` do) can leave it exactly rank-deficient –
with or without an intercept. That is exactly what explains the
near-zero eigenvalues in Ozkale (2021) and Ozkale and Abbasi (2022).

## References

Akram, M. N., B. M. G. Kibria, M. R. Abonazel, and N. Afzal. 2022. “On
the Performance of Some Biased Estimators in the Gamma Regression Model:
Simulation and Applications.” *Journal of Statistical Computation and
Simulation* 92 (12): 2425–47.
<https://doi.org/10.1080/00949655.2022.2032059>.

Amin, M., M. Qasim, M. Amanullah, and S. Afzal. 2020. “Performance of
Some Ridge Estimators for the Gamma Regression Model.” *Statistical
Papers* 61 (3): 997–1026. <https://doi.org/10.1007/s00362-017-0971-z>.

Bulut, Y. M. 2023. “Inverse Gaussian Liu-Type Estimator.”
*Communications in Statistics - Simulation and Computation* 52 (10):
4864–79. <https://doi.org/10.1080/03610918.2021.1971243>.

Chatterjee, S., and A. S. Hadi. 1988. *Sensitivity Analysis in Linear
Regression*. John Wiley; Sons.

Kurtoglu, F. 2021. “Modified Ridge Parameter Estimators for Log-Gamma
Model: Monte Carlo Evidence with a Graphical Investigation.”
*Communications in Statistics - Simulation and Computation* 50 (7):
2168–83. <https://doi.org/10.1080/03610918.2019.1650181>.

Lukman, A. F., K. Ayinde, B. M. G. Kibria, and E. T. Adewuyi. 2022.
“Modified Ridge-Type Estimator for the Gamma Regression Model.”
*Communications in Statistics - Simulation and Computation* 51 (9):
5009–23. <https://doi.org/10.1080/03610918.2020.1752720>.

Ozkale, M. R. 2021. “The Red Indicator and Corrected VIFs in Generalized
Linear Models.” *Communications in Statistics - Simulation and
Computation* 50 (12): 4144–70.
<https://doi.org/10.1080/03610918.2019.1639740>.

Ozkale, M. R., and A. Abbasi. 2022. “Iterative Restricted OK Estimator
in Generalized Linear Models and the Selection of Tuning Parameters via
MSE and Genetic Algorithm.” *Statistical Papers* 63 (6): 1979–2040.
<https://doi.org/10.1007/s00362-022-01304-0>.

Ozkale, M. R., and E. Arican. 2016. “A New Biased Estimator in Logistic
Regression Model.” *Statistics* 50 (2): 233–53.
<https://doi.org/10.1080/02331888.2015.1123711>.

Shewa, G. A., and F. I. Ugwuowo. 2023. “Kibria-Lukman Type Estimator for
Gamma Regression Model.” *Concurrency and Computation: Practice and
Experience* 35 (1). <https://doi.org/10.1002/cpe.7441>.
