# LMest 3.2.8 time-heterogeneous latent Markov model (modBasic = 0): one
# transition matrix per occasion; binary items. Three scenarios. Reference:
# JStats tests/fixtures/latent/lta-time-varying-ref.json.
source(file.path("equivalence", "lta-extensions", "common.R"))
ref <- jsonlite::fromJSON(file.path(jstats, "tests", "fixtures", "latent",
                                    "lta-time-varying-ref.json"))
rows <- lapply(names(ref$scenarios), function(name) {
  s <- ref$scenarios[[name]]
  long <- to_long(s$data)
  items <- paste0("y", seq_len(s$M))
  long[items] <- lapply(long[items], factor)
  fit <- lta(long, items, "subject", n_profiles = s$K, time = "occasion",
             categorical = items, transitions = "occasion", n_starts = 10,
             seed = 1, tol = 1e-10, max_iter = 3000)
  ours_rho <- vapply(fit$measurement[[1L]]$response_probabilities,
                     function(p) p[, 2L], numeric(s$K))
  # LMest rho is K x M x categories; category 2 is "1".
  perm <- align(ours_rho, matrix(s$expected$rho[, , 2L], s$K))
  table <- get_results(fit, "transitions")
  tau_gap <- max(vapply(seq_len(s$Ti - 1L), function(t) {
    slice <- subset(table, occasion == t + 1L)
    ours <- matrix(0, s$K, s$K)
    ours[cbind(match(slice$from, paste0("profile_", 1:s$K)),
               match(slice$to, paste0("profile_", 1:s$K)))] <- slice$probability
    max(abs(ours[perm, perm] - matrix(s$expected$tau[t, , ], s$K)))
  }, numeric(1)))
  initial <- get_results(fit, "initial")$probability[perm]
  rbind(
    row_result(paste0("lmest_", name), "log likelihood", fit$log_likelihood,
               s$expected$logLik, 1e-4),
    row_result(paste0("lmest_", name), "parameter count", fit$n_parameters,
               s$expected$npar, 0),
    row_result(paste0("lmest_", name), "max |transition diff| per occasion", tau_gap, 0, 1e-4),
    row_result(paste0("lmest_", name), "max |initial diff|",
               max(abs(initial - s$expected$piv)), 0, 1e-4))
})
result <- do.call(rbind, rows)
print(result, row.names = FALSE, digits = 8)
saveRDS(result, file.path("equivalence", "lta-extensions", "result-lmest-occasion.rds"))
