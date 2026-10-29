# Equivalence tests

Tests that compare latents with another implementation: `numDeriv` derivatives,
`mclust`, `MASS::glm.nb()` and `MASS::polr()`, `lm()`/`glm()` fits,
`mice::pool.scalar()`, `stats::optim()`/`optimize()`/`integrate()` used as
independent references, Latent GOLD, Mplus and glca results, published values,
reference implementations written in the tests (dense, brute-force and
path-enumeration likelihoods, the log-space recursion in
`helper-lta-log-reference.R`), and the numbers stated in the workflow articles.
They are for development only: this directory is excluded from the package
build, so they never run in `R CMD check` or on CRAN. CI runs them in the
`equivalence` workflow.

A file here has the name of the `tests/testthat/` file whose model it checks;
`tests/testthat/` keeps the package's own contracts, invariants and refusals.

Run them from the package root:

```r
Sys.setenv(NOT_CRAN = "true")   # some blocks carry skip_on_cran()
devtools::load_all()
testthat::test_dir("tests/equivalence", load_package = "none")
```

They need `numDeriv`, `mclust`, `MASS`, `mice`, `withr` and `knitr` installed.
The workflow claims read `vignettes/articles/`; the Latent GOLD step-three
comparison reads `equivalence/latentgold/`.
