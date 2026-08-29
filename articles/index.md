# Articles

### Collinearity reproductions

Published cases where collinearity was diagnosed without applying the
correct data transformation, reproduced and compared against
multiCollglm’s results.

- [Reproduction: cancer remission data (Lee
  1974)](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-lee-cancer.md):

  Reproduces the three logistic-model variants fitted to Lee’s (1974)
  cancer remission dataset – Lesaffre and Marx (1993), Ozkale and Arican
  (2016)/Ozkale (2021), and Huang, Jou and Cho (2015) – and contrasts
  them against multiCollglm’s five methods.

- [Reproduction: Chatterjee and Hadi (1988) -- Nitrogen
  data](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-nitrogen.md):

  Reproduces the collinearity diagnostic of the Nitrogen dataset (Gamma,
  inverse link), contrasts it against nine published
  reproductions/variants, and shows how each of multiCollglm’s five
  condition-number methods (RAW, MP, MS, WS, OZ) responds to this
  model’s structural singularity.

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
  the log link for the Gamma family that explains why two of the
  package’s methods coincide exactly in this case.

- [Reproduction: plastic plywood manufacturing data
  (Plastic)](https://cbgarciaugr.github.io/multiCollglm/articles/reproduction-plastic.md):

  Reproduces the collinearity diagnostic of the plastic plywood
  manufacturing dataset, used by six works fitting different count-data
  GLMs, and documents a citation error (CN = 8634.73) that propagates
  unverified across at least three publications.
