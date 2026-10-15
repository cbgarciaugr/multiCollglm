# Articles

### Get started

- [Introduction to
  multiCollglm](https://cbgarciaugr.github.io/multiCollglm/articles/multiCollglm.md):

### Examples

Real-data applications. Start with the overview, which links to each
example and summarizes the four condition-number definitions on all the
data sets.

- [Examples: real-data
  applications](https://cbgarciaugr.github.io/multiCollglm/articles/examples.md):

  Index of the worked examples of multiCollglm: five reproductions of
  published collinearity diagnostics on real data, and a comparison with
  the VisCollin package, with a summary of the four condition-number
  definitions on every data set.

- [Reproduction: cancer remission data (Lee
  1974)](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-lee-cancer.md):

  Reproduces the three logistic-model variants fitted to Lee’s (1974)
  cancer remission dataset – Lesaffre and Marx (1993), Ozkale and Arican
  (2016)/Ozkale (2021), and Huang, Jou and Cho (2015) – and contrasts
  them against multiCollglm’s four methods (WS, MP, X, RAW).

- [Reproduction: Chatterjee and Hadi (1988) -- Nitrogen
  data](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-nitrogen.md):

  Reproduces the collinearity diagnostic of the Nitrogen dataset (Gamma,
  inverse link), contrasts it against nine published
  reproductions/variants, and shows how multiCollglm’s four
  condition-number methods (WS, MP, X, RAW) behave, and why centering
  the weighted matrix is not possible for this model.

- [Reproduction: coal mine roof-fracture data
  (Mine)](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-mines.md):

  Reproduces the collinearity diagnostic of the coal mine roof-fracture
  dataset, used by five works fitting different count-data GLMs, and
  shows that the “RAW” (untransformed) method exactly reproduces the
  most-cited value in the literature.

- [Reproduction: Bodyfat data
  (TH.data)](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-bodyfat.md):

  Reproduces the collinearity diagnostic of the bodyfat dataset (Gamma,
  log link) used by five works, and documents a structural fact about
  the log link for the Gamma family that explains why three of the
  package’s methods coincide exactly in this case.

- [Reproduction: plastic plywood manufacturing data
  (Plastic)](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-plastic.md):

  Reproduces the collinearity diagnostic of the plastic plywood
  manufacturing dataset, used by six works fitting different count-data
  GLMs, and documents a citation error (CN = 8634.73) that propagates
  unverified across at least three publications.

### Comparison with other software

How multiCollglm relates to the VisCollin package, on the same data.

- [Comparison with
  VisCollin](https://cbgarciaugr.github.io/multiCollglm/articles/comparison-viscollin.md):

  Compares multiCollglm with the VisCollin package (Friendly), on the
  real datasets of the package: what VisCollin::colldiag() computes when
  it is given a GLM, when it coincides exactly with multiCollglm (method
  = “X”), and where the two diverge (IRLS weights, factors, weights and
  offsets, visualization).
