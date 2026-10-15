# Reproduction: plastic plywood manufacturing data (Plastic)

## Context

The `Plastic` dataset comes from a statistical-process-control case
study at a medium-sized timber industry that manufactures laminated
plastic plywood (2016). It has 100 reference observations; the response
`Y` is the number of defects per unit area (count data), and the four
regressors are volumetric shrinkage (`X1`, %), assembly time (`X2`,
minutes), wood density (`X3`, g/cm3) and drying temperature (`X4`, deg
C). It is a widely cited dataset in the collinearity literature because
`X1` and `X2` are strongly correlated with each other, as are `X3` and
`X4`.

It has been reused by Sami, Amin and Butt (2022), Tanis and Asar (2026),
Mammadova and Ozkale (2021), Dawoud (2025), Tanis and Asar (2024) and
Ashraf et al. (2025), each with a different count-data model (Bell,
Conway-Maxwell-Poisson, Poisson, negative binomial).

## Data

``` r

library(multiCollglm)
data(Plastic)
str(Plastic)
#> 'data.frame':    100 obs. of  5 variables:
#>  $ Y : int  18 1 41 53 3 24 26 0 4 4 ...
#>  $ X1: num  12.39 8.49 9.83 12.39 8.58 ...
#>  $ X2: num  17.4 14.2 15.1 17.9 14 ...
#>  $ X3: num  0.52 0.52 0.57 0.54 0.53 0.55 0.54 0.53 0.54 0.52 ...
#>  $ X4: num  118 106 150 124 118 ...

mod <- glm(Y ~ X1 + X2 + X3 + X4, family = poisson(link = "log"),
           data = Plastic)
```

## Diagnostic with multiCollglm and comparison with literature

Tanis and Asar (2026) and Tanis and Asar (2024) fit a Bell regression
without an intercept, with the data centered and standardized so that
\\X'X\\ is in correlation form (Tanis and Asar (2024) state this
explicitly: *“The design matrix is centered and standardized so that
\\X'X\\ is in correlation form”*), obtaining \\CN = 74.5281\\.

Sami, Amin and Butt (2022) fit a Conway-Maxwell-Poisson model and
compute their \\CN\\ as the square root of the eigenvalue ratio of
\\X'\hat{W}X\\ **excluding the intercept**, with no additional
transformation reported, obtaining \\CN = 8634.73\\. Ashraf et al.
(2025) report **exactly the same value**, \\8634.73\\, despite fitting
three different distributions (Conway-Maxwell-Poisson, Poisson and
negative binomial) – that three different models, necessarily with
different IRLS weight matrices, produce the same \\CN\\ to four decimal
places is itself a red flag.

Dawoud (2025) uses only three of the four regressors (assembly time,
density, temperature) with Conway-Maxwell-Poisson and a log link,
obtaining \\CN = 13120.14\\, with no further detail on the
transformation. Mammadova and Ozkale (2021) report \\CN =
2.9558\times10^{5}\\ as the ratio between the largest and smallest
eigenvalue of \\X'\hat{W}X\\, without specifying a distribution, link,
or transformation.

``` r

cat("All four variables, with intercept:\n")
#> All four variables, with intercept:
knitr::kable(calc_table(mod), digits = 4)
```

| method | nc_label |           CN |   sqrt(CN) |
|:-------|:---------|-------------:|-----------:|
| WS     | NC_WS    |     768847.2 |   876.8393 |
| MP     | NC_MP    |     808281.3 |   899.0447 |
| X      | NC_X     |     354909.6 |   595.7429 |
| RAW    | NC_RAW   | 5411312565.6 | 73561.6243 |

``` r


mod1 <- glm(Y ~ X1 + X2 + X3 + X4 - 1, family = poisson(link = "log"),
            data = Plastic)
cat("All four variables, without intercept:\n")
#> All four variables, without intercept:
knitr::kable(calc_table(mod1), digits = 4)
```

| method | nc_label |           CN |   sqrt(CN) |
|:-------|:---------|-------------:|-----------:|
| WS     | NC_WS    |     104612.3 |   323.4382 |
| MP     | NC_MP    |     106011.9 |   325.5947 |
| X      | NC_X     |     117479.0 |   342.7521 |
| RAW    | NC_RAW   | 306087439\.4 | 17495.3548 |

``` r


mod2 <- glm(Y ~ X2 + X3 + X4 - 1, family = poisson(link = "log"),
             data = Plastic)
cat("Dawoud's (2025) subset (X2, X3, X4), without intercept:\n")
#> Dawoud's (2025) subset (X2, X3, X4), without intercept:
knitr::kable(calc_table(mod2), digits = 4)
```

| method | nc_label |           CN |  sqrt(CN) |
|:-------|:---------|-------------:|----------:|
| WS     | NC_WS    |     1846.908 |   42.9757 |
| MP     | NC_MP    |     1964.449 |   44.3221 |
| X      | NC_X     |     1808.590 |   42.5275 |
| RAW    | NC_RAW   | 27574987.094 | 5251.1891 |

No combination of family/link/method available in
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)
reproduces the three most-cited values in this literature (\\8634.73\\;
\\13120.14\\; \\2.9558\times10^{5}\\), whether with the four variables
or Dawoud’s (2025) subset of three, and with or without an intercept.
The values that are obtained span a wide range (from a few tens to
billions depending on the method), but none match the published figures.

That, together with the fact that Ashraf et al. (2025) report the same
\\8634.73\\ for three different distributions (when the IRLS weight
matrix, and hence \\CN\\, should differ between them) suggests that
\\8634.73\\ has been **cited from one work to the next without being
independently recalculated**, rather than being a verified computation
in each case. We were unable to establish exactly which transformation
first produced that figure.

Tanis and Asar’s (2026) value (74.5281, Bell without an intercept, data
centered and standardized into correlation form) could not be
independently verified because Bell regression is not part of R base’s
[`glm()`](https://rdrr.io/r/stats/glm.html) families; it is left as
reported by the original source.

## References

Ashraf, B., M. Amin, W. Emam, Y. Tashkandy, and M. Faisal. 2025.
“Negative Binomial Regression Model Estimation Using Stein Approach:
Methods, Simulation, and Applications.” *Journal of Mathematics* 2025
(1): 9134821. <https://doi.org/10.1155/jom/9134821>.

Dawoud, I. 2025. “New Biased Estimators for the Conway-Maxwell-Poisson
Model.” *Journal of Statistical Computation and Simulation* 95 (1):
117–36. <https://doi.org/10.1080/00949655.2024.2421317>.

Mammadova, U., and M. R. Ozkale. 2021. “Profile Monitoring for Count
Data Using Poisson and Conway-Maxwell-Poisson Regression-Based Control
Charts Under Multicollinearity Problem.” *Journal of Computational and
Applied Mathematics* 388: 113275.
<https://doi.org/10.1016/j.cam.2020.113275>.

Marcondes Filho, D., and A. M. P. Sant’Anna. 2016. *Statistical Process
Control Applied to a Plastic Plywood Manufacturing Process*.

Sami, F., M. Amin, and M. M. Butt. 2022. “On the Ridge Estimation of the
Conway-Maxwell Poisson Regression Model with Multicollinearity: Methods
and Applications.” *Concurrency and Computation: Practice and
Experience* 34 (1): e6477. <https://doi.org/10.1002/cpe.6477>.

Tanis, C., and Y. Asar. 2024. “Almost Unbiased Ridge Estimator in Bell
Regression Model: Theory and Application to Plastic Polywood Data.”
*Statistics* 58 (5): 1031–45.
<https://doi.org/10.1080/02331888.2024.2389416>.

Tanis, C., and Y. Asar. 2026. “Almost Unbiased Liu Estimator in Bell
Regression Model.” *Communications in Statistics - Simulation and
Computation*, ahead of print.
<https://doi.org/10.1080/03610918.2026.2648869>.
