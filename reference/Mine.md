# Coal mine roof-fracture data

Counts of roof fractures recorded at coal mines in the Appalachian
region of West Virginia, together with four continuous predictors
describing the geometry and history of each mine. Used by several
published count-data GLM collinearity illustrations (Poisson,
Conway-Maxwell-Poisson, Bell and discrete-Weibull regressions).

## Usage

``` r
Mine
```

## Format

A data frame with 44 rows and 5 columns:

- y:

  Number of roof fractures (count response).

- x1:

  Inner burden thickness.

- x2:

  Percentage extraction of the previously mined lower seam.

- x3:

  Lower seam height.

- x4:

  Time since the mine was opened.

## Source

Myers, R.H. (1990). *Classical and Modern Regression with Applications*
(2nd ed.). Duxbury Press. Reproduced (correctly) by Marx, B.D. (1992).

## Details

Two transcription errors present in some circulating copies of this
dataset – row 30 with `y = 22` instead of `y = 2`, and row 36 with
`x1 = 0` instead of `x1 = 80` – have been corrected here after
cross-checking against a second paper's appendix data table. With the
correction, ordinary-least-squares coefficients, eigenvalues, condition
indices and the ridge parameter `k` all reproduce Marx (1992) to high
precision; the uncorrected version does not.

## Examples

``` r
data(Mine)
mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(link = "log"),
           data = Mine)
condition_number(mod, method = "WS")
#> Eigenvalues of (X'WX)_ul (decreasing order):
#> [1] 4.1168 0.4495 0.3411 0.0874 0.0050
#> 
#> CN (NC_WS) (eigenvalue scale): 816.9398
#> sqrt(CN) (NC_WS) (classical, singular-value scale): 28.5822
```
