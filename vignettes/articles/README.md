# Reproduction articles (pkgdown)

This folder does not contain package vignettes: these are **pkgdown
articles**. They are rendered only for the package website
(`pkgdown::build_site()` / `pkgdown::build_articles()`) and are **not**
compiled during `R CMD build`/`R CMD check` or installed with the package
— that's why `vignettes/articles` is listed in `.Rbuildignore`. This lets
each article use external datasets or take longer to run without
affecting the package's build/check or CRAN submission.

## How to add a new reproduction

1. Copy `reproduction-template.Rmd` with a descriptive name, e.g.:

   ```
   reproduction-perez-2019-energy-consumption.Rmd
   ```

2. Fill in the sections following the template (context, data, original
   diagnostic exactly as done in the published work, diagnostic with
   `multiCollglm`, comparison, conclusion, references).

3. Before including third-party data, check its license. If it cannot be
   redistributed, use a simulation faithful to its statistical properties
   and say so explicitly in the "Data" section.

4. Add the article to `_pkgdown.yml` (`articles` section, "Reproductions"
   group) so it appears in the website's menu.

5. Build the website locally to review the result before publishing:

   ```r
   pkgdown::build_site()
   # or just the articles:
   pkgdown::build_articles()
   ```
