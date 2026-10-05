# Equivalence tests

Tests that compare latents with another implementation: `numDeriv` derivatives,
`mclust`, `MASS::glm.nb()` and `MASS::polr()`, `glm()` fits, Latent GOLD and
Mplus results, the reference implementations in
`tests/testthat/helper-*-reference.R`, and the numbers stated in the workflow
articles. They are for development only: this directory is excluded from the
package build, so they never run in `R CMD check` or on CRAN.

Run them from the package root:

```r
devtools::load_all()
testthat::test_dir("tests/equivalence", load_package = "none")
```

They need `numDeriv`, `mclust`, `MASS` and `knitr` installed. The workflow
claims read `vignettes/articles/`.
