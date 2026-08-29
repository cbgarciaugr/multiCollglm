# Nitrogen dioxide air-pollution data (Chatterjee and Hadi 1988)

Weather conditions and nitrogen dioxide concentrations recorded over
September 1984 at a monitoring station in the San Francisco Bay Area,
originally reported by Chatterjee and Hadi (1988). Chatterjee and Hadi
describe the dataset as containing 26 observations, but the original
observation numbering jumps from observation 8 to observation 10, with
no observation 9 reported: the dataset actually has only 25
observations, which is what is included here.

## Usage

``` r
Nitrogen
```

## Format

An object of class `data.frame` with 25 rows and 5 columns.

## Source

Chatterjee, S. and Hadi, A.S. (1988). *Sensitivity Analysis in Linear
Regression*. John Wiley and Sons.

## Details

- y:

  Nitrogen dioxide concentration, in parts per hundred million
  (p.p.h.m.).

- x1:

  Mean wind speed, in miles per hour (m.p.h.).

- x2:

  Maximum temperature, in degrees Fahrenheit.

- x3:

  Insolation, in langleys per day.

- x4:

  Stability factor, in degrees Fahrenheit.

Used in
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md)'s
reproduction article for a GLM with a Gamma family and canonical
(inverse) link, where it also illustrates a structural singularity: for
that family/link combination, `sqrt(w) * eta` is exactly 1 for every
observation, so centering the IRLS-weighted design matrix (as in
`method = "OZ"` or `"MS"`) makes it exactly rank-deficient regardless of
whether an intercept is present.

## Examples

``` r
data(Nitrogen)
mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"),
           data = Nitrogen)
condition_number(mod, method = "WS")
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 4.7741 0.1578 0.0443 0.0230 0.0008
#> 
#> Condition number (NC_WS) (eigenvalue scale): 6080.4514
#> Condition index (NC_WS) (classical, sqrt): 77.9772
```
