library(latents)
env <- new.env(parent = asNamespace('latents'))
expr <- parse('tests/testthat/test-astra-diagnostics.R')
eval(expr[[1L]], env)
eval(expr[[2L]], env)
fit <- env$astra_diagnostic_fits()$cross
expected <- env$astra_hand_quality(fit$group_posteriors, fit$sampling_weights)$modal
actual <- get_results(fit, 'group_classes')$n_assigned
cat('actual counts:', actual, '; expected positive-support counts:', expected,
    '; nobs:', nobs(fit), '\n')
stopifnot(identical(actual, expected), sum(actual) == nobs(fit))
