pkgload::load_all('.', quiet = TRUE)
original <- new.env(parent = asNamespace('latents'))
sys.source('validation/release-fixes-2026-10-28/regression/before-inference.R.txt', original)
invisible(lapply(c('.mixture_pack', '.mixture_inference'), function(name) {
  assignInNamespace(name, get(name, original), ns = 'latents')
}))
result <- testthat::test_file('tests/testthat/test-release-audit-regression.R',
  reporter = 'summary', stop_on_failure = FALSE)
table <- as.data.frame(result)
print(table[c('test', 'failed', 'warning', 'passed')])
stopifnot(sum(table$failed) > 0L)
