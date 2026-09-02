#' multiCollglm: collinearity diagnostics for GLMs and glmnet
#'
#' Provides [condition_number()], [bkw_diagnostics()] and
#' [rvif_diagnostics()] to assess collinearity in generalized linear
#' models, all computed over the design matrix weighted by the IRLS
#' weights at convergence and scaled to unit length (without centering).
#' Works directly on `glm` objects, and on `glmnet`/`cv.glmnet` objects by
#' internally refitting a `glm` on the active set of variables at the
#' given lambda. [ml_collinearity()] adds the complementary diagnostic of
#' Lesaffre and Marx (1993), comparing the condition number of the
#' original design matrix against that of the IRLS-weighted information
#' matrix to distinguish ordinary collinearity from ML-collinearity.
#'
#' @keywords internal
"_PACKAGE"
