pkgload::load_all('.', quiet = TRUE)
source('equivalence/harness/registry.R')
original <- latents:::.mixture_pack
fixed <- function(spec, params) {
  params$sigma2 <- unname(params$sigma2)
  original(spec, params)
}
assignInNamespace('.mixture_pack', fixed, ns = 'latents')
run <- run_equivalence('flexmix', check_api = FALSE)
print(run)
print(as.data.frame(run, what = 'summary'))
write.csv(run$comparisons, 'validation/release-audit-2026-10-28/equivalence/flexmix-patched-comparisons.csv', row.names = FALSE)
write.csv(run$suites, 'validation/release-audit-2026-10-28/equivalence/flexmix-patched-suites.csv', row.names = FALSE)
stopifnot(.equivalence_passes(run, require_all = TRUE))
