pkgload::load_all(".", quiet = TRUE)
set.seed(110508)
n <- 240L
T <- 4L
p <- 5L
state <- matrix(1L, n, T)
state[, 1] <- sample(1:2, n, TRUE)
invisible(lapply(2:3, function(t) {
  state[, t] <<- ifelse(runif(n) < .2, 3L - state[, t - 1], state[, t - 1])
}))
state[, T] <- ifelse(runif(n) < .5, 3L, state[, T - 1L])
means <- rbind(rep(0, p), c(2.5, 2.5, 0, 0, 0), rep(5, p))
vars <- paste0("y", seq_len(p))
dd <- do.call(rbind, lapply(seq_len(T), function(t) {
  values <- means[state[, t], , drop = FALSE] + matrix(rnorm(n * p), n, p)
  colnames(values) <- vars
  data.frame(id = seq_len(n), time = t, values)
}))
warnings <- character()
fit <- withCallingHandlers(lta(dd, vars, "id", 3, "time", n_starts = 4, seed = 32,
  tol = 1e-10, cores = 1), warning = function(w) {
    warnings <<- c(warnings, paste(class(w)[1], conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
print(warnings)
top <- which.max(rowSums(fit$means))
tab <- as.data.frame(get_results(fit, "transitions"))
print(tab[tab$from == top, c("from", "to", "probability", "expected_count", "estimated")], digits = 12)
cat("origin total:", sum(tab$expected_count[tab$from == top]), "n_parameters:", fit$n_parameters,
    "boundary:", fit$boundary, "\n")
general_warnings <- character()
general <- withCallingHandlers(lta(dd, vars, "id", 3, "time", transitions = "occasion",
  n_starts = 3, seed = 32, tol = 1e-10, cores = 1), warning = function(w) {
    general_warnings <<- c(general_warnings, paste(class(w)[1], conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
print(general_warnings)
top_general <- which.max(rowSums(general$measurement[[1]]$means))
tab_general <- as.data.frame(get_results(general, "transitions"))
print(tab_general[tab_general$from == paste0("profile_", top_general), ], digits = 12)
totals <- vapply(general$expectation$moments[[1]]$pairs, function(moves) sum(moves[, top_general, ]), 0)
cat("General origin totals per move:", totals, "n_parameters:", general$n_parameters, "\n")
saveRDS(list(fit = fit, table = tab, general = general, general_table = tab_general,
  totals = totals, warnings = warnings, general_warnings = general_warnings),
  "validation/review-2026-11-05/engines/near-empty.rds")
