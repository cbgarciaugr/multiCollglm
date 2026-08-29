#' Nitrogen dioxide air-pollution data (Chatterjee and Hadi 1988)
#'
#' Weather conditions and nitrogen dioxide concentrations recorded over
#' September 1984 at a monitoring station in the San Francisco Bay Area,
#' originally reported by Chatterjee and Hadi (1988). Chatterjee and Hadi
#' describe the dataset as containing 26 observations, but the original
#' observation numbering jumps from observation 8 to observation 10, with
#' no observation 9 reported: the dataset actually has only 25
#' observations, which is what is included here.
#'
#' Used in [condition_number()]'s reproduction article for a GLM with a
#' Gamma family and canonical (inverse) link, where it also illustrates a
#' structural singularity: for that family/link combination,
#' `sqrt(w) * eta` is exactly 1 for every observation, so centering the
#' IRLS-weighted design matrix (as in `method = "OZ"` or `"MS"`) makes it
#' exactly rank-deficient regardless of whether an intercept is present.
#'
#' @format A data frame with 25 rows and 5 columns:
#' \describe{
#'   \item{y}{Nitrogen dioxide concentration, in parts per hundred million
#'     (p.p.h.m.).}
#'   \item{x1}{Mean wind speed, in miles per hour (m.p.h.).}
#'   \item{x2}{Maximum temperature, in degrees Fahrenheit.}
#'   \item{x3}{Insolation, in langleys per day.}
#'   \item{x4}{Stability factor, in degrees Fahrenheit.}
#' }
#' @source Chatterjee, S. and Hadi, A.S. (1988). \emph{Sensitivity Analysis
#'   in Linear Regression}. John Wiley and Sons.
#' @examples
#' data(Nitrogen)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = Gamma(link = "inverse"),
#'            data = Nitrogen)
#' condition_number(mod, method = "WS")
"Nitrogen"

#' Coal mine roof-fracture data
#'
#' Counts of roof fractures recorded at coal mines in the Appalachian
#' region of West Virginia, together with four continuous predictors
#' describing the geometry and history of each mine. Used by several
#' published count-data GLM collinearity illustrations (Poisson,
#' Conway-Maxwell-Poisson, Bell and discrete-Weibull regressions).
#'
#' Two transcription errors present in some circulating copies of this
#' dataset -- row 30 with `y = 22` instead of `y = 2`, and row 36 with
#' `x1 = 0` instead of `x1 = 80` -- have been corrected here after
#' cross-checking against a second paper's appendix data table. With the
#' correction, ordinary-least-squares coefficients, eigenvalues, condition
#' indices and the ridge parameter `k` all reproduce Marx (1992) to high
#' precision; the uncorrected version does not.
#'
#' @format A data frame with 44 rows and 5 columns:
#' \describe{
#'   \item{y}{Number of roof fractures (count response).}
#'   \item{x1}{Inner burden thickness.}
#'   \item{x2}{Percentage extraction of the previously mined lower seam.}
#'   \item{x3}{Lower seam height.}
#'   \item{x4}{Time since the mine was opened.}
#' }
#' @source Myers, R.H. (1990). \emph{Classical and Modern Regression with
#'   Applications} (2nd ed.). Duxbury Press. Reproduced (correctly) by
#'   Marx, B.D. (1992).
#' @examples
#' data(Mine)
#' mod <- glm(y ~ x1 + x2 + x3 + x4, family = poisson(link = "log"),
#'            data = Mine)
#' condition_number(mod, method = "WS")
"Mine"

#' Lee (1974) cancer remission data
#'
#' Clinical observations from 27 leukemia patients, of whom 9 experienced
#' complete cancer remission (`y = 1`) and 18 did not (`y = 0`), together
#' with six quantitative covariates. A classic teaching example for
#' collinearity diagnostics in binomial logistic regression, reproduced
#' (with varying subsets of covariates) by Lesaffre and Marx (1993),
#' Ozkale (2019) and others; also distributed as the `Remission` dataset
#' in SAS/STAT example documentation. Cross-checked against the SAS
#' documentation page for this example and against the original SAS
#' `datalines`: exact match.
#'
#' @format A data frame with 27 rows and 7 columns:
#' \describe{
#'   \item{y}{Complete remission indicator (1 = complete remission,
#'     0 = incomplete remission).}
#'   \item{x1}{CELL, cellularity index.}
#'   \item{x2}{SMEAR, smear index.}
#'   \item{x3}{INFIL, infiltrate index.}
#'   \item{x4}{LI, labeling index.}
#'   \item{x5}{BLAST, percentage of blast cells.}
#'   \item{x6}{TEMP, highest body temperature (recorded on a /100 scale,
#'     as in the original source; some reproduction articles multiply it
#'     back by 100 to match a paper's own units).}
#' }
#' @source Lee, E.T. (1974). A computer program for linear logistic
#'   regression analysis. \emph{Computer Programs in Biomedicine}, 4(2),
#'   80-92. \doi{10.1016/0010-468X(74)90011-7}
#' @examples
#' data(LeeCancer)
#' mod <- glm(y ~ x1 + x4 + x6, family = binomial(), data = LeeCancer)
#' condition_number(mod, method = "WS")
"LeeCancer"

#' Plastic plywood manufacturing process data
#'
#' 100 reference observations from a medium-sized timber industry
#' manufacturing laminated plastic plywood, reported by Marcondes Filho
#' and Sant'Anna (2016) as a statistical-process-control case study and
#' subsequently reused as a count-data GLM collinearity example (Poisson,
#' Conway-Maxwell-Poisson, Bell and negative-binomial regressions). The
#' input variables show strong, well-documented pairwise correlations
#' (shrinkage with assembly time; density with drying temperature).
#' Verified byte-identical (up to formatting) against the values reported
#' in the 2016 primary source.
#'
#' @format A data frame with 100 rows and 5 columns:
#' \describe{
#'   \item{Y}{Number of surface defects per laminated-plywood area
#'     (count response).}
#'   \item{X1}{Volumetric shrinkage, in percent.}
#'   \item{X2}{Assembly time, in minutes.}
#'   \item{X3}{Wood density, in g/cm^3.}
#'   \item{X4}{Drying temperature, in degrees Celsius.}
#' }
#' @source Marcondes Filho, D. and Sant'Anna, A.M.P. (2016). Statistical
#'   process control applied to a plastic plywood manufacturing process.
#' @examples
#' data(Plastic)
#' mod <- glm(Y ~ X1 + X2 + X3 + X4, family = poisson(link = "log"),
#'            data = Plastic)
#' condition_number(mod, method = "WS")
"Plastic"
