# multiCollglm 0.4.0

## Breaking changes

* The default definition is now Weissfeld and Sereika (1991) (`method = "WS"`):
  the IRLS-weighted design matrix scaled to unit length by column, **without
  centering**. (This is what `center = FALSE, scale = "unit"`, the previous
  default, computed.)
* `method` now takes four values: `"WS"` (default), `"MP"` (MacKinnon and
  Puterman, 1989), `"X"` (new) and `"RAW"`.
* **Removed** `method = "OZ"` (Ozkale, 2021) and `method = "MS"` (Marx and
  Smith, 1990). Both centered the data, which drops the intercept from the
  analysis (so non-essential collinearity cannot be detected) and, for
  some families such as Gamma with the inverse link, produces an exactly
  singular matrix. Calling them stops with an explanatory message.
* **Removed** the arguments `center` and `scale` of `condition_number()`,
  `bkw_diagnostics()` and `rvif_diagnostics()`. Only `method` selects the
  matrix now. Passing `center` or `scale` stops with an explanatory message.

## New features

* `method = "X"`: eigenvalues of X_ul'X_ul, the *original* design matrix
  scaled to unit length, without the IRLS weights (kappa_X of Lesaffre
  and Marx, 1993). For a model with intercept it is exactly
  `VisCollin::colldiag(mod, add.intercept = TRUE)`.
* `method` is available in `condition_number()`, `bkw_diagnostics()` and
  `rvif_diagnostics()` (and their `glmnet` / `cv.glmnet` methods), and is
  reported in the printed output and in the returned object.
* New article *Comparison with VisCollin*, run on the package's real datasets.
* pkgdown website: "Examples" menu with links to every example, an examples
  overview article, and a GitHub Actions workflow that builds and publishes the
  site to the `gh-pages` branch.

## Other changes

* The `multiColl` package is no longer used: `method = "MP"` and
  `ml_collinearity()` are computed natively. `multiColl` has been removed
  from `Suggests`; `VisCollin` has been added (used only by the comparison
  article and one test).
* A zero eigenvalue is now detected with a numerical-rank tolerance
  (`max(eigenvalue) * p * .Machine$double.eps`), so a perfectly collinear
  design gives `Inf` and a warning rather than an arbitrary huge number.
* The reproduction articles use the four methods; the Nitrogen article
  shows, with base R code, why centering makes the weighted matrix singular
  for the Gamma family with inverse link.
* Tests rewritten: each method against its manual definition, Table 1 of the
  manuscript on the packaged datasets, equivalence with `VisCollin::colldiag()`.

# multiCollglm 0.3.0

* Initial documented version: `condition_number()` with five definitions
  (`"RAW"`, `"MP"`, `"MS"`, `"WS"`, `"OZ"`), `bkw_diagnostics()`,
  `rvif_diagnostics()`, `ml_collinearity()`, `glmnet` / `cv.glmnet` methods,
  vignette and pkgdown articles.
