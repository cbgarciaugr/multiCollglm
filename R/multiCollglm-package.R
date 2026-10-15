#' multiCollglm: collinearity diagnostics for GLMs and glmnet
#'
#' Provides [condition_number()], [bkw_diagnostics()] and
#' [rvif_diagnostics()] to assess collinearity in generalized linear
#' models. By default they are computed, following Weissfeld and Sereika
#' (1991), over the design matrix weighted by the IRLS weights at
#' convergence and scaled to unit length (without centering); the `method`
#' argument also offers the definition of MacKinnon and Puterman (1989),
#' the one of the original (unweighted) design matrix, and the raw
#' information matrix.
#' Works directly on `glm` objects, and on `glmnet`/`cv.glmnet` objects by
#' internally refitting a `glm` on the active set of variables at the
#' given lambda. [ml_collinearity()] adds the complementary diagnostic of
#' Lesaffre and Marx (1993), comparing the condition number of the
#' original design matrix against that of the IRLS-weighted information
#' matrix to distinguish ordinary collinearity from ML-collinearity.
#'
#' @keywords internal
"_PACKAGE"
