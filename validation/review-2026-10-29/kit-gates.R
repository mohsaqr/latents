# Review reproduction: operate only on copies of saved tables in memory.
# Run from the latents repository root. No source or saved evidence is changed.
library(testthat)
root <- getwd()
exercise <- function(kit, test, compare_name, mutation) {
  env <- new.env(parent = globalenv())
  env$comparison <- read.csv(file.path(root, 'equivalence', kit, 'comparison.csv'))
  env$comparison <- mutation(env$comparison)
  env[[compare_name]] <- function(...) env$comparison
  env$mplus_growth_dir <- env$stata_dir <- env$pool_packages_dir <- file.path(root, 'equivalence', kit)
  expression <- tail(parse(file.path(root, 'equivalence', test)), 1)[[1]]
  cat('\n', kit, ': rows=', nrow(env$comparison), '; NA agreements=', sum(is.na(env$comparison$agrees)), '\n', sep='')
  # Evaluate the original wrapper's test_that expression with its expensive
  # comparison function replaced by the mutated saved comparison table.
  eval(expression, env)
}
exercise('mplus-growth', 'test-mplus-growth.R', 'compare_all_mplus_growth', function(d) {
  ix <- startsWith(d$quantity, 's.e.')
  d$latents[ix] <- d$difference[ix] <- NA_real_
  d$agrees[ix] <- NA
  d
})
exercise('stata-published', 'test-stata-published.R', 'compare_all_stata', function(d) {
  ix <- startsWith(d$quantity, 's.e.')
  d$latents[ix] <- d$difference[ix] <- NA_real_
  d$agrees[ix] <- NA
  d
})
exercise('pool-packages', 'test-pool-packages.R', 'compare_all_pool', function(d) {
  d$latents <- d$difference <- NA_real_
  d$agrees <- NA
  d
})
source('equivalence/package-kit-runner.R')
d <- read.csv('equivalence/mplus-growth/comparison.csv')
k <- read.csv('equivalence/mplus-growth/known-differences.csv')
ix <- d$case == 'runner_gmm_k3' & d$quantity == 'number of parameters'
stopifnot(sum(ix) == 1L)
d$latents[ix] <- 9999
d$difference[ix] <- abs(d$latents[ix] - d$reference_value[ix])
d$agrees[ix] <- FALSE
z <- apply_register(d, k)
cat('\nParameter-count mutation:\n')
print(z$comparison[ix, c('case', 'quantity', 'reference_value', 'latents', 'status')])
cat('unexplained=', sum(z$comparison$status == 'disagrees'), ' stale=', nrow(z$stale), '\n')
stopifnot(z$comparison$status[ix] == 'known difference',
          sum(z$comparison$status == 'disagrees') == 0L,
          nrow(z$stale) == 0L)
