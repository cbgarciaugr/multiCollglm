#' multiCollglm: collinearity diagnostics for GLMs and glmnet
#'
#' Provides [condition_number()], [bkw_diagnostics()] and
#' [rvif_diagnostics()] to assess collinearity in generalized linear
#' models, all computed over the design matrix weighted by the IRLS
#' weights at convergence and scaled to unit length (without centering).
#' Works directly on `glm` objects, and on `glmnet`/`cv.glmnet` objects by
#' internally refitting a `glm` on the active set of variables at the
#' given lambda.
#'
#' @keywords internal
"_PACKAGE"
