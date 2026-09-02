# Reproduction: cancer remission data (Lee 1974)

## Context

Lee’s (1974) cancer remission dataset records clinical observations from
27 patients, of whom 9 had complete remission (`y = 1`) and 18 had
incomplete remission (`y = 0`), together with six quantitative
covariates: cellularity index (CELL), smear index (SMEAR), infiltrate
index (INFIL), labeling index (LI), percentage of blast cells (BLAST)
and body temperature (TEMP). It is a classic example for collinearity
diagnostics in binomial logistic regression, used with different subsets
of covariates by Lesaffre and Marx (1993), Ozkale and Arican (2016)
together with Ozkale (2021), and Huang, Jou and Cho (2015).

Lesaffre and Marx (1993) use this dataset only with three covariates
(LI, TEMP, CELL) that show severe collinearity. The remaining reviewed
works use all covariates except, in some cases, LI.

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

## Diagnostic with multiCollglm and comparison with literature

### Model 1: Lesaffre and Marx (1993) – CELL, LI, TEMP

``` r

mod1_leecancer <- glm(y ~ x1 + x4 + x6100, family = binomial(), data = LeeCancer)
cat("Model 1 (Lesaffre and Marx), with intercept:\n")
#> Model 1 (Lesaffre and Marx), with intercept:
knitr::kable(calc_table(mod1_leecancer), digits = 4)
```

| method | nc_label |           CN |   sqrt(CN) |
|:-------|:---------|-------------:|-----------:|
| RAW    | NC_RAW   | 1.168157e+08 | 10808.1330 |
| MP     | NC_MP    | 1.088698e+05 |   329.9542 |
| MS     | NC_MS    | 2.567380e+02 |    16.0230 |
| WS     | NC_WS    | 9.931110e+04 |   315.1366 |
| OZ     | NC_OZ    | 2.788063e+02 |    16.6975 |

Lesaffre and Marx (1993) used this example to compare the $`CN`$
proposed by Mackinnon and Puterman (1989) (329.95) with the one proposed
by Weissfeld and Sereika (1991) (315.14) and also with the one obtained
directly from design matrix $`X`$ with unit length transformation
(190.78). Results obtained with `multiCollglm` match those reported in
the original paper.

This is in fact the very example Lesaffre and Marx (1993) use (Sec. 6.1)
to introduce their own
[`ml_collinearity()`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
criterion, comparing `sqrt(CN)` of the original `X` against `sqrt(CN)`
of the Mackinnon-Puterman information matrix:

``` r

ml_collinearity(mod1_leecancer)
#> ML-collinearity diagnostic (Lesaffre and Marx, 1993)
#> 
#> sqrt(CN) of X (original design matrix, multiColl::CNs(X)[1]): 11.7696
#> sqrt(CN) of W (MacKinnon-Puterman information matrix, method = "MP"): 329.9542
#> ratio (sqrt(CN)_W / sqrt(CN)_X): 28.0345
#> 
#> Collinearity in X (sqrt(CN)_X > 30): no
#> ML-collinearity (ratio_wx > 5 and sqrt(CN)_W > 30): YES
```

`multiCollglm` reproduces both of their reported values exactly
(`sqrt(CN)_X = 190.78`, `sqrt(CN)_W = 329.95`), and their own conclusion
for this model: the ratio (28.03) stays well below the
`ratio_threshold = 5` cutoff, so this is a case of ordinary collinearity
in `X` (`sqrt(CN)_X = 190.78 > 30`), not ML-collinearity.

### Model 2: Ozkale and Arican (2016) / Ozkale (2021) – CELL, SMEAR, INFIL, BLAST, TEMP (no LI)

``` r

mod2_leecancer <- glm(y ~ x1 + x2 + x3 + x5 + x6100, family = binomial(),
                data = LeeCancer)
cat("Model 2 (Ozkale), with intercept:\n")
#> Model 2 (Ozkale), with intercept:
knitr::kable(calc_table(mod2_leecancer), digits = 4)
```

| method | nc_label |            CN |   sqrt(CN) |
|:-------|:---------|--------------:|-----------:|
| RAW    | NC_RAW   | 219838806\.31 | 14826.9621 |
| MP     | NC_MP    |     163327.20 |   404.1376 |
| MS     | NC_MS    |      23296.31 |   152.6313 |
| WS     | NC_WS    |     143044.11 |   378.2117 |
| OZ     | NC_OZ    |      14371.56 |   119.8814 |

``` r

mod2b_leecancer <- glm(y ~ x1 + x2 + x3 + x5 + x6100 - 1, family = binomial(),
                 data = LeeCancer)
cat("Model 2 (Ozkale), without intercept (their own specification):\n")
#> Model 2 (Ozkale), without intercept (their own specification):
knitr::kable(calc_table(mod2b_leecancer), digits = 4)
```

| method | nc_label |           CN |   sqrt(CN) |
|:-------|:---------|-------------:|-----------:|
| RAW    | NC_RAW   | 1.942571e+08 | 13937.6155 |
| MP     | NC_MP    | 7.769953e+04 |   278.7464 |
| MS     | NC_MS    | 4.202298e+02 |    20.4995 |
| WS     | NC_WS    | 7.476822e+04 |   273.4378 |
| OZ     | NC_OZ    | 1.643415e+04 |   128.1957 |

Ozkale and Arican (2016) standardize and center the explanatory
variables before fitting, explicitly stating that “standardization
removes any non-essential ill-conditioning resulting from the
intercept”. They obtain $`CN = 295.7026`$ (eigenvalue ratio, no square
root), which on the square-root scale would be $`17.196`$. This value
does not match any of
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)’s
five methods

Ozkale (2021), with the same model, reports the untransformed
eigenvalues $`15.6562, 0.4601, 0.1557, 0.0145, 0.0001`$ and from those a
condition number of $`1.3032\times10^{5}`$, a value that is obtained
exactly neither as the ratio ($`15.6562/0.0001 = 156,562`$) nor as its
square root ($`\approx 395.7`$) from those very eigenvalues she reports.
The same happens with her centered-and-scaled scenario: eigenvalues
$`4.4133,0.3657, 0.2009, 0.0200, 0.0002`$, from which the ratio would be
$`22,066.5`$, not the $`2.6099\times10^{4}`$ she reports. This value
does not match any of
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)’s
five methods

### Model 3: Huang, Jou and Cho (2015) (all covariates)

``` r

mod3_leecancer <- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6100, family = binomial(),
              data = LeeCancer)
cat("Model 3 (Huang, Jou and Cho), with intercept:\n")
#> Model 3 (Huang, Jou and Cho), with intercept:
knitr::kable(calc_table(mod3_leecancer), digits = 4)
```

| method | nc_label |            CN |   sqrt(CN) |
|:-------|:---------|--------------:|-----------:|
| RAW    | NC_RAW   | 395524494\.29 | 19887.7976 |
| MP     | NC_MP    |     248614.00 |   498.6121 |
| MS     | NC_MS    |      43277.14 |   208.0316 |
| WS     | NC_WS    |     218108.47 |   467.0208 |
| OZ     | NC_OZ    |      33948.94 |   184.2524 |

``` r

mod3b_leecancer <- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6100 - 1, family = binomial(),
               data = LeeCancer)
cat("Model 3 (Huang, Jou and Cho), without intercept:\n")
#> Model 3 (Huang, Jou and Cho), without intercept:
knitr::kable(calc_table(mod3b_leecancer), digits = 4)
```

| method | nc_label |           CN |   sqrt(CN) |
|:-------|:---------|-------------:|-----------:|
| RAW    | NC_RAW   | 2.900940e+08 | 17032.1463 |
| MP     | NC_MP    | 1.431963e+05 |   378.4129 |
| MS     | NC_MS    | 5.158097e+02 |    22.7114 |
| WS     | NC_WS    | 1.393080e+05 |   373.2398 |
| OZ     | NC_OZ    | 3.372989e+04 |   183.6570 |

Huang, Jou and Cho (2015) use all six covariates and report
$`CN = 418.13`$. This value does not match any of
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)’s
five methods

## References

Huang, C.-C. L., Y.-J. Jou, and H.-J. Cho. 2015. “A New
Multicollinearity Diagnostic for Generalized Linear Models.” *Journal of
Applied Statistics*, ahead of print.
<https://doi.org/10.1080/02664763.2015.1126239>.

Lee, E. T. 1974. “A Computer Program for Linear Logistic Regression
Analysis.” *Computer Programs in Biomedicine* 4 (2): 80–92.
<https://doi.org/10.1016/0010-468X(74)90011-7>.

Lesaffre, E., and B. D. Marx. 1993. “Collinearity in Generalized Linear
Regression.” *Communications in Statistics - Theory and Methods* 22 (7):
1933–52. <https://doi.org/10.1080/03610929308831126>.

Mackinnon, M. J., and M. L. Puterman. 1989. “Collinearity in Generalized
Linear Models.” *Communications in Statistics - Theory and Methods* 18
(9): 3463–72. <https://doi.org/10.1080/03610928908830102>.

Ozkale, M. R. 2021. “The Red Indicator and Corrected VIFs in Generalized
Linear Models.” *Communications in Statistics - Simulation and
Computation* 50 (12): 4144–70.
<https://doi.org/10.1080/03610918.2019.1639740>.

Ozkale, M. R., and E. Arican. 2016. “A New Biased Estimator in Logistic
Regression Model.” *Statistics* 50 (2): 233–53.
<https://doi.org/10.1080/02331888.2015.1123711>.

Weissfeld, L. A., and S. M. Sereika. 1991. “A Multicollinearity
Diagnostic for Generalized Linear Models.” *Communications in
Statistics - Theory and Methods* 20 (4): 1183–98.
<https://doi.org/10.1080/03610929108830558>.
