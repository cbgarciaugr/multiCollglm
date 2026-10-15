# Prepares the four bundled example datasets used throughout the package's
# reproduction articles (vignettes/articles/reproduccion-*.Rmd). Not run at
# build/check/install time (data-raw/ is listed in .Rbuildignore); this
# script only documents how data/*.rda were produced and lets them be
# regenerated from the original CSV snapshots kept in this same folder.
#
# The fifth dataset used in the articles, Bodyfat (Penrose, Nelson and
# Fisher 1985), is NOT bundled here: it already ships in the 'TH.data'
# package (data("bodyfat", package = "TH.data")), so re-distributing it
# here would be redundant and raise an unnecessary licensing question.

Nitrogen <- read.csv("data-raw/Nitrogen.csv")
# Chatterjee and Hadi (1988): San Francisco Bay Area, September 1984.
# y = NO2 concentration (p.p.h.m.); x1 = mean wind speed (m.p.h.);
# x2 = maximum temperature (F); x3 = insolation (langleys/day);
# x4 = stability factor (F). Chatterjee and Hadi describe the dataset as
# having 26 observations, but the original observation numbering jumps
# from 8 to 10 (no observation 9): only 25 observations actually exist,
# which is what this data frame contains.
usethis_or_save <- function(obj, name) {
  assign(name, obj)
  save(list = name, file = file.path("data", paste0(name, ".rda")),
       compress = "bzip2")
}
usethis_or_save(Nitrogen, "Nitrogen")

Mine <- read.csv("data-raw/Mine.csv")
# Coal mine roof-fracture data (Myers 1990; used by Marx 1992 and others
# for GLM collinearity/ridge illustrations). y = number of fractures;
# x1 = inner burden thickness; x2 = percentage extraction of the
# previously mined lower seam; x3 = lower seam height; x4 = time since
# the mine opened. Two transcription errors present in some circulating
# copies of this dataset (row 30: y = 22 should be y = 2; row 36: x1 = 0
# should be x1 = 80) have been corrected here after cross-checking against
# a second paper's appendix table; with the correction, OLS coefficients,
# eigenvalues, condition indices and ridge-k all reproduce Marx (1992) to
# high precision.
usethis_or_save(Mine, "Mine")

LeeCancer <- read.csv("data-raw/LeeCancer.csv")
# Lee (1974) cancer remission data, as compiled by Lee and widely
# reproduced since (e.g. SAS/STAT example 'Remission', Lesaffre and Marx
# 1993, Ozkale 2019). 27 patients; y = complete remission indicator
# (1 = yes, 9 patients; 0 = no, 18 patients). x1 = CELL (cellularity
# index), x2 = SMEAR (smear index), x3 = INFIL (infiltrate index),
# x4 = LI (labeling index), x5 = BLAST (percentage of blast cells),
# x6 = TEMP (highest body temperature, originally recorded /100 -- the
# reproduction articles multiply it back by 100 where a paper's own
# variable is on that scale). Cross-checked against the SAS documentation
# page for this example and against user-supplied SAS 'datalines': exact
# match.
usethis_or_save(LeeCancer, "LeeCancer")

Plastic <- read.csv("data-raw/Plastic.csv")
# Plastic plywood manufacturing process data (Marcondes Filho and
# Sant'Anna 2016), 100 reference observations. Y = number of surface
# defects per unit area (count); X1 = volumetric shrinkage (%);
# X2 = assembly time (minutes); X3 = wood density (g/cm^3); X4 = drying
# temperature (deg C). Verified byte-identical (up to formatting) against
# the values reported in the original 2016 source.
usethis_or_save(Plastic, "Plastic")
