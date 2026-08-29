# Lee (1974) cancer remission data

Clinical observations from 27 leukemia patients, of whom 9 experienced
complete cancer remission (`y = 1`) and 18 did not (`y = 0`), together
with six quantitative covariates. A classic teaching example for
collinearity diagnostics in binomial logistic regression, reproduced
(with varying subsets of covariates) by Lesaffre and Marx (1993), Ozkale
(2019) and others; also distributed as the `Remission` dataset in
SAS/STAT example documentation. Cross-checked against the SAS
documentation page for this example and against the original SAS
`datalines`: exact match.

## Usage

``` r
LeeCancer
```

## Format

An object of class `data.frame` with 27 rows and 7 columns.

## Source

Lee, E.T. (1974). A computer program for linear logistic regression
analysis. *Computer Programs in Biomedicine*, 4(2), 80-92.
[doi:10.1016/0010-468X(74)90011-7](https://doi.org/10.1016/0010-468X%2874%2990011-7)

## Details

- y:

  Complete remission indicator (1 = complete remission, 0 = incomplete
  remission).

- x1:

  CELL, cellularity index.

- x2:

  SMEAR, smear index.

- x3:

  INFIL, infiltrate index.

- x4:

  LI, labeling index.

- x5:

  BLAST, percentage of blast cells.

- x6:

  TEMP, highest body temperature (recorded on a /100 scale, as in the
  original source; some reproduction articles multiply it back by 100 to
  match a paper's own units).

## Examples

``` r
data(LeeCancer)
mod <- glm(y ~ x1 + x4 + x6, family = binomial(), data = LeeCancer)
condition_number(mod, method = "WS")
#> Eigenvalues of X_unit'X_unit (decreasing order):
#> [1] 3.9034 0.0927 0.0038 0.0000
#> 
#> Condition number (NC_WS) (eigenvalue scale): 99311.0952
#> Condition index (NC_WS) (classical, sqrt): 315.1366
```
