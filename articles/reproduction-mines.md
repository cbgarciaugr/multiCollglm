# Reproduction: coal mine roof-fracture data (Mine)

## Context

The `Mine` dataset records 44 observations of roof fractures in coal
mines in the Appalachian region of West Virginia, with four continuous
predictors: inner burden thickness, percentage extraction of the
previously mined lower seam, lower seam height, and time since the mine
opened. It is the real-data application shared by five reviewed works,
each fitting a different count-data GLM: Kurtoglu and Ozkale (2019b) and
Kurtoglu and Ozkale (2019a) (same authors, a restricted ridge and a
restricted Liu estimator respectively, both with Poisson and a log
link); Sami, Amin and Butt (2022) (Conway-Maxwell-Poisson, a “COMPRR”
ridge estimator); Hussein and Algamal (2026) (Bell regression, a
Kibria-Lukman-type estimator); and Sami et al. (2026) (discrete Weibull,
a “DWRE” ridge estimator).

Sami, Amin and Butt (2022) state that they analyze a different
distribution but provide no diagnostic measures of their own (which
should differ across distributions), instead citing a previous paper to
claim that the dataset exhibits collinearity. Sami et al. (2026) do the
same, citing Sami, Amin and Butt (2022) despite also using a different
distribution; moreover, based on the code in Sami et al.’s (2026)
appendix, the diagnostic appears to be computed on the unweighted
$`X'X`$, and although the code includes a standardization step, it is
never actually applied.

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
```

## Diagnostic with multiCollglm and comparison with literature

`sqrt(CN)` (square-root scale) of `method = "RAW"` – no transformation
at all, directly on $`X'\hat{W}X`$ – is **2093.3158**, which matches, up
to rounding, the **2093.32** reported by Kurtoglu and Ozkale (2019b,
2019a).

Hussein and Algamal (2026), with a Bell distribution, obtain **2139.89**
(same order of magnitude, consistent with a different but Poisson-like
count distribution for these data).

``` r

knitr::kable(tbl, digits = 4)
```

| method | nc_label |           CN |  sqrt(CN) |
|:-------|:---------|-------------:|----------:|
| RAW    | NC_RAW   | 4381970.9542 | 2093.3158 |
| MP     | NC_MP    |     699.8016 |   26.4538 |
| MS     | NC_MS    |     101.5642 |   10.0779 |
| WS     | NC_WS    |     816.9398 |   28.5822 |
| OZ     | NC_OZ    |       4.3874 |    2.0946 |

This dataset is a good example of how the transformation of data can
affect the results in `CN`. Note that the diagnostic’s magnitude changes
substantially.

## References

Hussein, S. M. J. A., and Z. Y. Algamal. 2026. “A New Shrinkage
Estimator for the Bell Regression Model.” *Thailand Statistician* 24
(2): 462–73.

Kurtoglu, F., and M. R. Ozkale. 2019a. “Restricted Liu Estimator in
Generalized Linear Models: Monte Carlo Simulation Studies on Gamma and
Poisson Distributed Responses.” *Hacettepe Journal of Mathematics and
Statistics* 48 (4): 1250–76. <https://doi.org/10.15672/HJMS.2018.624>.

Kurtoglu, F., and M. R. Ozkale. 2019b. “Restricted Ridge Estimator in
Generalized Linear Models: Monte Carlo Simulation Studies on Poisson and
Binomial Distributed Responses.” *Communications in Statistics -
Simulation and Computation* 48 (4): 1191–218.
<https://doi.org/10.1080/03610918.2017.1408822>.

Sami, F., H. M. Aljohani, K. Iqbal, S. Ghorashi, and M. Abu Bakar
Siddique. 2026. “A New Estimator for Discrete Weibull Regression Model
with Correlated Variables: An Illustrative Example.” *Journal of
Statistical Computation and Simulation* 96 (11): 2737–57.
<https://doi.org/10.1080/00949655.2026.2636232>.

Sami, F., M. Amin, and M. M. Butt. 2022. “On the Ridge Estimation of the
Conway-Maxwell Poisson Regression Model with Multicollinearity: Methods
and Applications.” *Concurrency and Computation: Practice and
Experience* 34 (1): e6477. <https://doi.org/10.1002/cpe.6477>.
