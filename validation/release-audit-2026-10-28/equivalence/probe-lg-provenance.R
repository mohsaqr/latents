source('equivalence/latentgold-wider/compare-functions.R')
directory <- 'equivalence/latentgold-wider'
cases <- readRDS(file.path(directory, 'kit/cases.rds'))
checks <- lapply(c('INPUTS', 'OUTPUTS'), function(which) {
  recorded <- read.table(file.path(directory, 'returned', paste0(which, '.md5')), col.names = c('md5', 'file'))
  prefix <- if (which == 'INPUTS') 'kit' else 'returned'
  matched <- unname(tools::md5sum(file.path(directory, prefix, recorded$file))) == recorded$md5
  data.frame(type = which, checked = length(matched), matched = sum(matched))
})
print(do.call(rbind, checks))
data_matches <- vapply(names(cases), function(name) {
  dat <- read.table(file.path(directory, 'kit', paste0(name, '.dat')), header = TRUE)
  isTRUE(all.equal(dat, cases[[name]]$data, check.attributes = FALSE, tolerance = 1e-14))
}, logical(1))
print(data_matches)
stopifnot(all(data_matches))
listing <- read_listing(file.path(directory, 'returned/r01.lst'))
theta <- setNames(listing$parameters$estimate, listing$parameters$term)
variance <- diag(listing$parameters$standard_error^2)
maps <- list(identity)
baseline <- compare_case(listing, maps, theta, variance, listing$log_likelihood, listing$n_parameters, 'self-test')
stopifnot(all(baseline$agrees))
perturbed <- theta
perturbed[1L] <- perturbed[1L] + 0.01
negative <- compare_case(listing, maps, perturbed, variance, listing$log_likelihood, listing$n_parameters, 'negative-control')
print(subset(negative, !agrees))
stopifnot(sum(!negative$agrees) == 1L)
cat('Output hashes and .rds/.dat parity checked here; production test checks input hashes only.\n')
