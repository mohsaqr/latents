# dentistry_k3: lca() and glca agree on the log likelihood, the shares and the
# response probabilities, and on Linux their posteriors differed by more than
# 1e-4. Which side is further from the maximum? Each solution is scored by a
# plain EM written here (binary items, three classes): its log likelihood, the
# size of one EM step taken from it, and its distance to the point that EM
# reaches from it after a long run.
#
# From the project root: Rscript equivalence/lca-packages/dentistry-posteriors.R
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("equivalence", "lca-packages", "compare-functions.R"))
case <- Filter(function(case) identical(case$name, "dentistry_k3"), lca_cases())[[1L]]
items <- as.matrix(case$data[case$items])
stopifnot("dentistry items are coded 0/1" = all(items %in% 0:1))

# Parameters in one shape: shares (K) and P(item = 1 | class) (K x items).
as_parameters <- function(solution) {
  list(shares = solution$shares,
       p_one = vapply(case$items, function(v) solution$probabilities[[v]][, 2L],
                      numeric(length(solution$shares))))
}
joint <- function(parameters) {
  log_density <- items %*% t(log(parameters$p_one)) + (1 - items) %*% t(log1p(-parameters$p_one))
  sweep(log_density, 2L, log(parameters$shares), `+`)
}
log_likelihood <- function(parameters) {
  values <- joint(parameters)
  top <- apply(values, 1L, max)
  sum(top + log(rowSums(exp(values - top))))
}
posteriors <- function(parameters) {
  values <- joint(parameters)
  weights <- exp(values - apply(values, 1L, max))
  weights / rowSums(weights)
}
em_step <- function(parameters) {
  responsibility <- posteriors(parameters)
  list(shares = colMeans(responsibility),
       p_one = crossprod(responsibility, items) / colSums(responsibility))
}
gap <- function(a, b) max(abs(a$shares - b$shares), abs(a$p_one - b$p_one))
# EM has no closed form: the iteration is the algorithm, bounded and reported.
long_run <- function(parameters, iterations = 400000L) {
  step <- NA_real_
  for (i in seq_len(iterations)) {
    updated <- em_step(parameters)
    step <- gap(updated, parameters)
    parameters <- updated
    if (step < 1e-15) break
  }
  list(parameters = parameters, iterations = i, last_step = step)
}

ours <- fit_latents_lca(case)
theirs <- fit_glca(case)
perm <- .match_classes(ours, theirs)
theirs_matched <- list(shares = theirs$shares[perm],
                       p_one = as_parameters(theirs)$p_one[perm, , drop = FALSE])
solutions <- list(latents = as_parameters(ours), glca = theirs_matched)
limits <- lapply(solutions, long_run)
options(digits = 12, width = 160)
report <- do.call(rbind, lapply(names(solutions), function(name) {
  start <- solutions[[name]]
  limit <- limits[[name]]
  data.frame(
    fit = name, log_likelihood = log_likelihood(start),
    one_em_step = gap(em_step(start), start),
    limit_log_likelihood = log_likelihood(limit$parameters),
    em_iterations = limit$iterations, last_step = limit$last_step,
    parameters_to_limit = gap(start, limit$parameters),
    posteriors_to_limit = max(abs(posteriors(start) - posteriors(limit$parameters))))
}))
print(report, row.names = FALSE)
cat(sprintf("\nlatents against glca: parameters %.3g, posteriors %.3g (reported: %.3g)\n",
            gap(solutions$latents, solutions$glca),
            max(abs(posteriors(solutions$latents) - posteriors(solutions$glca))),
            max(abs(ours$posteriors - theirs$posteriors[, perm, drop = FALSE]))))
cat(sprintf("the two limits: parameters %.3g, posteriors %.3g, log likelihood %.3g\n",
            gap(limits$latents$parameters, limits$glca$parameters),
            max(abs(posteriors(limits$latents$parameters) - posteriors(limits$glca$parameters))),
            log_likelihood(limits$latents$parameters) - log_likelihood(limits$glca$parameters)))
print(utils::packageVersion("glca"))
