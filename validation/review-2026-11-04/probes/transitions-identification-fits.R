# What lta() does on specifications the path distribution cannot identify.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(n, T, seed) {
  set.seed(seed)
  P <- matrix(c(.8,.2,.3,.7), 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    s <- integer(T); s[1] <- sample(2, 1, prob = c(.6,.4))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 3)[s]), y2 = rnorm(T, c(1, -2)[s], .8), y3 = rnorm(T, c(0, 2)[s]))
  }))
}
vars <- c("y1","y2","y3")
run <- function(label, data, ...) {
  w <- character()
  f <- withCallingHandlers(lta(data, vars, "id", n_profiles = 2, time = "time", n_starts = 6, seed = 1, ...),
    warning = function(cnd) { w <<- c(w, class(cnd)[1]); invokeRestart("muffleWarning") })
  inf <- tryCatch({ pi <- parameter_inference(f); sprintf("SEs given; median SE %.3g, max %.3g", median(pi$standard_error, na.rm = TRUE), max(pi$standard_error, na.rm = TRUE)) },
                  error = function(e) paste("refused:", class(e)[1]))
  cat(sprintf("\n== %s\n logLik %.5f  n_parameters %d  BIC %.2f  conv %s  replicated %d/%d\n warnings: %s\n inference: %s\n",
     label, f$log_likelihood, f$n_parameters, f$bic, f$converged, f$n_best_replicated, nrow(f$starts),
     paste(unique(w), collapse = ", "), inf))
  invisible(f)
}
d2 <- sim(300, 2, 1)
f1 <- run("T=2, H=1 (saturated in the latent paths)", d2)
f2 <- run("T=2, H=2", d2, n_group_classes = 2)
print(round(f2$group_probabilities,3)); print(round(f2$transition_probabilities,3))
f2b <- suppressWarnings(lta(d2, vars, "id", 2, "time", n_group_classes = 2, n_starts = 6, seed = 2))
cat("other seed H=2: logLik", f2b$log_likelihood, "; shares", round(f2b$group_probabilities, 3), "\n")
print(round(f2b$transition_probabilities, 3))
f3 <- run("T=2, stayer", d2, mover_stayer = TRUE)
print(get_results(f3, "group_classes")); print(subset(get_results(f3, "transitions"), from != to))
f3b <- suppressWarnings(lta(d2, vars, "id", 2, "time", mover_stayer = TRUE, n_starts = 6, seed = 2))
cat("other seed stayer: logLik", f3b$log_likelihood, "; shares", round(f3b$group_probabilities, 3), "\n")
print(subset(get_results(f3b, "transitions"), from != to))
d3 <- sim(300, 3, 2)
f4 <- run("T=3, H=1, saturated order 2", d3, order = 2)
f5 <- run("T=3, H=2, order 2", d3, order = 2, n_group_classes = 2)
f6 <- run("T=3, order 2, stayer", d3, order = 2, mover_stayer = TRUE)
f7 <- run("T=3, H=2, occasion transitions", d3, transitions = "occasion", n_group_classes = 2)
f8 <- run("T=3, H=1, occasion transitions", d3, transitions = "occasion")
f9 <- run("T=3, H=2 homogeneous, identified", d3, n_group_classes = 2)
g <- enumerate_classes(d2, vars, "id", time = "time", n_profiles = 2, n_group_classes = 1:2, n_starts = 3, seed = 1)
print(get_results(g))
