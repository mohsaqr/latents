pkgload::load_all('.',quiet=TRUE)
original <- .lta_coefficient_table
# Reproduce the original defect in memory: the table discards the caller's
# already-computed information and recomputes at its own default step.
mutant <- function(fit, block, level, vcov_type, inference = NULL) {
  original(fit, block, level, vcov_type, inference = NULL)
}
assignInNamespace('.lta_coefficient_table',mutant,'latents')
testthat::test_file('tests/testthat/test-engine-fixes.R')
