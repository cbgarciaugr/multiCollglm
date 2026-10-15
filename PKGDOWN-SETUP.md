# Publishing the pkgdown website

The website is built by `.github/workflows/pkgdown.yaml` and published
to the `gh-pages` branch. One-time setup:

1.  **Source on `main`.** The package source (R/, man/, vignettes/,
    DESCRIPTION, \_pkgdown.yml, .github/…) must live on the `main`
    branch. `gh-pages` must only contain the generated site (the
    workflow creates it).
2.  **Workflow.** Keep `.github/workflows/pkgdown.yaml` on `main`.
3.  **GitHub Pages.** Settings \> Pages \> Build and deployment \>
    Source: *Deploy from a branch*, Branch: `gh-pages`, folder
    `/ (root)`.
4.  **Permissions.** Settings \> Actions \> General \> Workflow
    permissions: *Read and write permissions*.
5.  Push to `main` (or run the workflow from the *Actions* tab). After
    it ends in green, the site is at
    <https://cbgarciaugr.github.io/multiCollglm/>.

Build it locally with:

``` r

install.packages(c("pkgdown", "VisCollin", "TH.data", "rvif", "glmnet"))
pkgdown::build_site()      # writes docs/ (ignored by git and by R CMD build)
```

To add an example: copy `vignettes/articles/reproduction-template.Rmd`,
then list it in `_pkgdown.yml` (the `articles:` section and the
`Examples` menu in `navbar:`) and in the table of
`vignettes/articles/examples.Rmd`.
