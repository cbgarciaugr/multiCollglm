# multiCollglm: collinearity diagnostics for GLMs and glmnet

Provides
[`condition_number()`](https://cbgarciaugr.github.io/multiCollglm/reference/condition_number.md),
[`bkw_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/bkw_diagnostics.md)
and
[`rvif_diagnostics()`](https://cbgarciaugr.github.io/multiCollglm/reference/rvif_diagnostics.md)
to assess collinearity in generalized linear models, all computed over
the design matrix weighted by the IRLS weights at convergence and scaled
to unit length (without centering). Works directly on `glm` objects, and
on `glmnet`/`cv.glmnet` objects by internally refitting a `glm` on the
active set of variables at the given lambda.
[`ml_collinearity()`](https://cbgarciaugr.github.io/multiCollglm/reference/ml_collinearity.md)
adds the complementary diagnostic of Lesaffre and Marx (1993), comparing
the condition number of the original design matrix against that of the
IRLS-weighted information matrix to distinguish ordinary collinearity
from ML-collinearity.

## See also

Useful links:

- <https://github.com/cbgarciaugr/multiCollglm>

- <https://cbgarciaugr.github.io/multiCollglm/>

- Report bugs at <https://github.com/cbgarciaugr/multiCollglm/issues>

## Author

**Maintainer**: Catalina Garcia Garcia <cbgarcia@go.ugr.es>

Authors:

- Catalina Garcia Garcia <cbgarcia@go.ugr.es>

- Roman Salmeron Gomez
