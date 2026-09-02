# Plastic plywood manufacturing process data

100 reference observations from a medium-sized timber industry
manufacturing laminated plastic plywood, reported by Marcondes Filho and
Sant'Anna (2016) as a statistical-process-control case study and
subsequently reused as a count-data GLM collinearity example (Poisson,
Conway-Maxwell-Poisson, Bell and negative-binomial regressions). The
input variables show strong, well-documented pairwise correlations
(shrinkage with assembly time; density with drying temperature).
Verified byte-identical (up to formatting) against the values reported
in the 2016 primary source.

## Usage

``` r
Plastic
```

## Format

An object of class `data.frame` with 100 rows and 5 columns.

## Source

Marcondes Filho, D. and Sant'Anna, A.M.P. (2016). Statistical process
control applied to a plastic plywood manufacturing process.

## Details

- Y:

  Number of surface defects per laminated-plywood area (count response).

- X1:

  Volumetric shrinkage, in percent.

- X2:

  Assembly time, in minutes.

- X3:

  Wood density, in g/cm^3.

- X4:

  Drying temperature, in degrees Celsius.

## Examples

``` r
data(Plastic)
mod <- glm(Y ~ X1 + X2 + X3 + X4, family = poisson(link = "log"),
           data = Plastic)
condition_number(mod, method = "WS")
#> Eigenvalues of (X'WX)_ul (decreasing order):
#> [1] 4.9872 0.0080 0.0046 0.0001 0.0000
#> 
#> CN (NC_WS) (eigenvalue scale): 768847.1927
#> sqrt(CN) (NC_WS) (classical, singular-value scale): 876.8393
```
