# A state reached only at the final wave, moderate separation, K = 3, T = 4.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
n <- 300; T <- 4; vars <- paste0("y", 1:5)
mu <- rbind(c(0,0,0,0,0), c(2.5,2.5,0,0,0), c(5,5,5,5,5))
set.seed(3)
state <- matrix(1L, n, T); state[, 1] <- sample(1:2, n, TRUE)
for (t in 2:(T-1)) state[, t] <- ifelse(runif(n) < .2, 3L - state[, t-1], state[, t-1])
state[, T] <- ifelse(runif(n) < .5, 3L, state[, T-1])   # state 3 appears only at wave 4
d <- do.call(rbind, lapply(seq_len(T), function(t) {
  y <- mu[state[, t], ] + matrix(rnorm(n * 5), n, 5); colnames(y) <- vars
  data.frame(id = seq_len(n), time = t, y)
}))
for (seed in 1:3) {
  w <- character()
  f <- withCallingHandlers(lta(d, vars, "id", 3, "time", n_starts = 5, seed = seed),
    warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
  top <- which.max(rowSums(f$means))
  cat(sprintf("\nhomogeneous seed %d: logLik %.4f n_par %d; warnings: %s\n", seed, f$log_likelihood, f$n_parameters, paste(unique(w), collapse = ", ")))
  print(subset(get_results(f, "transitions"), from == top)[, c("from","to","probability","expected_count","estimated")])
}
for (seed in 1:3) {
  w <- character()
  g <- withCallingHandlers(lta(d, vars, "id", 3, "time", weights = NULL, initial_covariates = character(), mover_stayer = FALSE, order = 1, transitions = "occasion", n_starts = 5, seed = seed),
    warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
  top <- which.max(rowSums(g$measurement[[1]]$means))
  cat(sprintf("\ngeneral (occasion) seed %d: logLik %.4f n_par %d; warnings: %s\n", seed, g$log_likelihood, g$n_parameters, paste(unique(w), collapse = ", ")))
  print(subset(get_results(g, "transitions"), from == paste0("profile_", top)))
  tot <- vapply(g$expectation$moments[[1]]$pairs, function(p) sum(p[, top, ]), 0)
  cat("expected moves out of the top profile, per move:", signif(tot, 3), "\n")
}
