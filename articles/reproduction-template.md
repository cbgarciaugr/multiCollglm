# Reproduction: \[Author(s), Year\] -- \[short title of the work\]

## Context

\[2-3 sentence summary of the original work: what GLM they fitted, which
variables, and what conclusion they drew about collinearity. Full
citation in the References section at the end.\]

## Data

``` r

# TODO: replace with the real data (or a faithful simulation) from the
# original work. Illustrative example with simulated data:
set.seed(1)
n <- 200
x1 <- rnorm(n)
x2 <- x1 + rnorm(n, sd = 0.05) # almost collinear with x1
x3 <- rnorm(n)
mu <- exp(1 + 0.3 * x1 + 0.3 * x2 - 0.2 * x3)
y <- rgamma(n, shape = 5, rate = 5 / mu)

mod <- glm(y ~ x1 + x2 + x3, family = Gamma(link = "inverse"))
summary(mod)
#> 
#> Call:
#> glm(formula = y ~ x1 + x2 + x3, family = Gamma(link = "inverse"))
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.399412   0.014927  26.758  < 2e-16 ***
#> x1          -0.138631   0.189221  -0.733    0.465    
#> x2          -0.003724   0.188467  -0.020    0.984    
#> x3           0.046692   0.009657   4.835 2.69e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.2649019)
#> 
#>     Null deviance: 125.402  on 199  degrees of freedom
#> Residual deviance:  54.882  on 196  degrees of freedom
#> AIC: 683.06
#> 
#> Number of Fisher Scoring iterations: 5
```

## Original diagnostic (without the correct transformation)

\[Explain what transformation (or lack of one) the original work
applied: for example, classical VIF on the unweighted matrix, or data
centered before computing the condition number.\]

``` r

# TODO: reproduce the calculation exactly as the original work did it.
# Illustrative example: classical least-squares VIF (car::vif), which
# ignores the GLM's IRLS weights and is therefore not the correct
# transformation for a non-Gaussian GLM.
if (requireNamespace("car", quietly = TRUE)) {
  car::vif(mod)
} else {
  message("The 'car' package is not installed; skipping this reference VIF.")
}
#>         x1         x2         x3 
#> 613.707129 612.280115   1.059993
```

\[Interpret the result exactly as the original work interpreted it: did
they conclude there was, or was not, problematic collinearity?\]

## Diagnostic with multiCollglm (correct transformation)

``` r

library(multiCollglm)

condition_number(mod)
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 3.2960 0.5203 0.1835 0.0002
#> 
#> Condition number (eigenvalue scale): 15364.1028
#> Condition index (classical, sqrt): 123.9520
bkw_diagnostics(mod)
#> Collinearity diagnostics (Belsley-Kuh-Welsch)
#> 
#>      Eigenvalue Sing.value Condition_index
#> dim1      3.296      1.815           1.000
#> dim2      0.520      0.721           2.517
#> dim3      0.184      0.428           4.238
#> dim4      0.000      0.015         123.952
#> 
#> Variance-decomposition proportions (row = coefficient, column = component):
#>              dim1  dim2  dim3  dim4
#> (Intercept) 0.020 0.011 0.969 0.001
#> x1          0.000 0.000 0.000 1.000
#> x2          0.000 0.000 0.000 1.000
#> x3          0.031 0.920 0.027 0.022
#> 
#> >> Possible collinearity problems (condition index >= 10 and proportion >= 0.5 on >= 2 variables):
#>   - dim4 (index = 123.95): x1, x2
rvif_diagnostics(mod)
#> Redefined Variance Inflation Factor (RVIF)
#> 
#>                 RVIF      %
#> (Intercept)    4.036 75.220
#> x1          2335.173 99.957
#> x2          2328.223 99.957
#> x3             1.708 41.462
```

## Comparison

\[Table or text comparing both diagnostics: does the conclusion about
the existence/severity of collinearity change? Which variables are
implicated according to each method?\]

|  | Original diagnostic | multiCollglm |
|----|----|----|
| Transformation applied | \[e.g. none / centered\] | Weighted by IRLS, scaled to unit length, without centering |
| Conclusion about collinearity | \[…\] | \[…\] |
| Variables implicated | \[…\] | \[…\] |

## Conclusion

\[2-4 sentences: why the transformation matters in this specific case,
and what would have changed in the original work’s conclusions had it
been applied correctly.\]

## References

- \[Full citation of the original work, in APA/the format used by the
  package.\]
- Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). *Regression
  Diagnostics: Identifying Influential Data and Sources of
  Collinearity*. Wiley.
- Salmeron, R., Garcia, C.B. and Garcia, J. (2025). A redefined Variance
  Inflation Factor: overcoming the limitations of the Variance Inflation
  Factor. *Computational Economics*, 65, 337-363.
  <https://doi.org/10.1007/s10614-024-10575-8>
