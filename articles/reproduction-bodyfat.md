# Reproduction: Bodyfat data (TH.data)

## Context

The `bodyfat` dataset from the `TH.data` package records
body-composition and anthropometric measurements from 71 healthy German
women. The response, `DEXfat`, is the percentage body fat measured by
dual-energy X-ray absorptiometry (DXA); the regressors include age,
waist and hip circumference, elbow and knee breadth, and several derived
anthropometric variables (`anthro3a`, `anthro3b`, `anthro3c`,
`anthro4`). The original goal is to predict body fat from simpler
anthropometric measurements.

Five of the reviewed works use this dataset under a Gamma distribution:
Almuqrin and AbaOud (2025) \[`review23`\], Al-Ghamdi et al. (2025)
\[`review44`\], Asar and Korkmaz (2022) \[`review51`\], Seifollahi et
al. (2024) \[`review76`\] and El-Masry et al. (2025) \[`review127`\].
All of them state that they compute the condition number as the square
root of the ratio between the largest and smallest eigenvalue of
$`X'\hat{W}X`$.

## Data

``` r

library(multiCollglm)
if (!requireNamespace("TH.data", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
  message("The 'TH.data' package is not installed; skipping this article.")
}
data("bodyfat", package = "TH.data")
str(bodyfat)
#> 'data.frame':    71 obs. of  10 variables:
#>  $ age         : num  57 65 59 58 60 61 56 60 58 62 ...
#>  $ DEXfat      : num  41.7 43.3 35.4 22.8 36.4 ...
#>  $ waistcirc   : num  100 99.5 96 72 89.5 83.5 81 89 80 79 ...
#>  $ hipcirc     : num  112 116.5 108.5 96.5 100.5 ...
#>  $ elbowbreadth: num  7.1 6.5 6.2 6.1 7.1 6.5 6.9 6.2 6.4 7 ...
#>  $ kneebreadth : num  9.4 8.9 8.9 9.2 10 8.8 8.9 8.5 8.8 8.8 ...
#>  $ anthro3a    : num  4.42 4.63 4.12 4.03 4.24 3.55 4.14 4.04 3.91 3.66 ...
#>  $ anthro3b    : num  4.95 5.01 4.74 4.48 4.68 4.06 4.52 4.7 4.32 4.21 ...
#>  $ anthro3c    : num  4.5 4.48 4.6 3.91 4.15 3.64 4.31 4.47 3.47 3.6 ...
#>  $ anthro4     : num  6.13 6.37 5.82 5.66 5.91 5.14 5.69 5.7 5.49 5.25 ...

mod <- glm(DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
             anthro3a + anthro3b + anthro3c + anthro4,
           family = Gamma(link = "log"), data = bodyfat)
summary(mod)
#> 
#> Call:
#> glm(formula = DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + 
#>     kneebreadth + anthro3a + anthro3b + anthro3c + anthro4, family = Gamma(link = "log"), 
#>     data = bodyfat)
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  -0.219467   0.233499  -0.940   0.3510    
#> age           0.001606   0.001001   1.605   0.1138    
#> waistcirc     0.004196   0.002086   2.012   0.0487 *  
#> hipcirc       0.010772   0.002497   4.315 5.95e-05 ***
#> elbowbreadth  0.014118   0.031775   0.444   0.6584    
#> kneebreadth   0.042557   0.022519   1.890   0.0635 .  
#> anthro3a     -0.122724   0.161763  -0.759   0.4510    
#> anthro3b      0.140431   0.175752   0.799   0.4274    
#> anthro3c      0.129183   0.064844   1.992   0.0508 .  
#> anthro4       0.163679   0.201576   0.812   0.4199    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.01038517)
#> 
#>     Null deviance: 9.55963  on 70  degrees of freedom
#> Residual deviance: 0.60848  on 61  degrees of freedom
#> AIC: 362.74
#> 
#> Number of Fisher Scoring iterations: 4
```

## Original diagnostic (a tour of the literature)

Almuqrin and AbaOud (2025), Al-Ghamdi et al. (2025) and El-Masry et al.
(2025) report an identical value, $`CN = 3394.30`$; none of them
explicitly states which link function was used, so we could not
reproduce that specific value. Asar and Korkmaz (2022) do specify a log
link and untransformed data, and obtain $`CN = 3410.63`$. Seifollahi et
al. (2024) also claims to use a log link but reports a different value,
$`CN =
4026.23`$ – an inconsistency within the literature itself, since two
works that claim to follow exactly the same procedure (Gamma, log,
untransformed) do not agree with each other.

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
| RAW    | NC_RAW   |     1.163241e+07 |       3410.6317 |
| MP     | NC_MP    |     1.239331e+05 |        352.0413 |
| MS     | NC_MS    |     1.196367e+04 |        109.3785 |
| WS     | NC_WS    |     1.239331e+05 |        352.0413 |
| OZ     | NC_OZ    |     9.811259e+02 |         31.3229 |

A notable feature of this table: `"MP"` and `"WS"` agree **exactly**
(not just approximately). The reason is structural, specific to this
family/link combination: for a Gamma GLM with a **log** link, the IRLS
weight of every observation at convergence is

``` math
w_i = \frac{1}{V(\mu_i)}\left(\frac{d\mu_i}{d\eta_i}\right)^2 = \frac{1}{\mu_i^2}\cdot\mu_i^2 = 1,
```

i.e. **constant and equal to 1 for every observation**, regardless of
the fitted values. This can be checked directly:

``` r

range(mod$weights)
#> [1] 1 1
```

With IRLS weights identically equal to 1, weighting by $`\hat{W}^{1/2}`$
does nothing (it’s the same as the unweighted matrix), so rescaling to
unit length *before* fitting (`"MP"`) or *after* fitting (`"WS"`)
produces exactly the same matrix, and hence the same condition number.
This is the log-link analogue, for Gamma, of the identity already
documented for Gamma with the canonical inverse link
(`sqrt(w) * eta ≡ 1`, see
[`?condition_number`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)):
every family/link combination has its own structural identity in the
IRLS weights, and it is worth keeping in mind when interpreting any
weighted collinearity diagnostic.

## Comparison

|  | Literature | multiCollglm |
|----|----|----|
| Untransformed (`"RAW"`) | $`CN = 3410.63`$ (Asar and Korkmaz 2022, square root) | 3410.6317 |
| Rescaled to unit length, no centering | $`CN \approx 352.04`$ | 352.0413 |
| Centered + scaled (`"OZ"`) | – | 31.3229 |

Asar and Korkmaz’s (2022) value reproduces with `method = "RAW"` to four
decimal places: **3410.6317** versus their **3410.63**. The “consensus”
value from the other three works (3394.30) comes very close but not
exactly, which suggests they used a slightly different link (or some
minor sampling variant) that they did not document; and Seifollahi et
al.’s (2024) value, despite claiming the same procedure as Asar and
Korkmaz, is close to neither.

## Conclusion

Once the already IRLS-weighted matrix is rescaled to unit length without
centering (`"WS"`), the condition number drops from the thousands to a
few hundred: the conclusion of severe collinearity holds, but is no
longer dominated by the disparity in scale between variables as
different as age (in years) and elbow breadth (in centimeters). The more
interesting finding for this dataset, however, is the structural
identity of the IRLS weights under Gamma with a log link: because they
are always equal to 1, two of the package’s five methods coincide
exactly, which does not happen in general for other family/link
combinations (see the Nitrogen article, where none of the five methods
agree with each other).

## References

- Penrose, K.W., Nelson, A.G. and Fisher, A.G. (1985). Generalized body
  composition prediction equation for men using simple measurement
  techniques. *Medicine & Science in Sports & Exercise*, 17(2), 189.
- Almuqrin, M.A. and AbaOud, M. (2025). Developing a new modified
  two-parameter Liu estimator for the gamma regression model: Method,
  simulation and application to health data. *Alexandria Engineering
  Journal*, 129, 1212-1222. <https://doi.org/10.1016/j.aej.2025.08.033>
- Al-Ghamdi, M.N., Abonazel, M.R., Dawoud, I., Algamal, Z.Y. and Azazy,
  A.R. (2025). A new estimator of the gamma regression model: theory,
  simulation, and application to body fat data. *Communications in
  Mathematical Biology and Neuroscience*, 2025, Article 53.
  <https://doi.org/10.28919/cmbn/9149>
- Asar, Y. and Korkmaz, M. (2022). Almost unbiased Liu-type estimators
  in gamma regression model. *Journal of Computational and Applied
  Mathematics*, 403, 113819. <https://doi.org/10.1016/j.cam.2021.113819>
- Seifollahi, S., Bevrani, H. and Kamary, K. (2024). Inequality
  restricted estimator for gamma regression: Bayesian approach as a
  solution to the multicollinearity. *Communications in Statistics -
  Theory and Methods*, 53(23), 8297-8311.
  <https://doi.org/10.1080/03610926.2023.2281267>
- El-Masry, A.M., Abu-Hamed, A.A., Abonazel, M.R., Omara, T.M. and
  Khattab, I.G. (2025). A new gamma regression estimate of female body
  fat percentage: advanced statistical modeling. *Communications in
  Mathematical Biology and Neuroscience*, 2025, Article 133.
  <https://doi.org/10.28919/cmbn/9493>
