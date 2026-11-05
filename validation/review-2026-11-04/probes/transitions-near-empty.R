# A profile occupied only at the last occasion, with well separated indicators.
# The docs promise such a row is uniform, estimated = FALSE and not counted.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
n <- 200
for (sep in c(6, 12, 30, 60)) {
  set.seed(1)
  s1 <- rep(1L, n)
  s2 <- ifelse(runif(n) < .4, 2L, 1L)
  d <- rbind(data.frame(id = 1:n, time = 1, y1 = rnorm(n, c(0, sep)[s1]), y2 = rnorm(n, c(0, sep)[s1])),
             data.frame(id = 1:n, time = 2, y1 = rnorm(n, c(0, sep)[s2]), y2 = rnorm(n, c(0, sep)[s2])))
  w <- character()
  f <- withCallingHandlers(lta(d, c("y1","y2"), "id", 2, "time", n_starts = 2, seed = 1),
    warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
  high <- which.max(f$means[, "y1"])
  cat(sprintf("\nsep %g: high profile = %d; n_par %d; warnings: %s\n", sep, high, f$n_parameters, paste(unique(w), collapse=", ")))
  tb <- get_results(f, "transitions")
  print(subset(tb, from == high))
  print(round(f$initial_probabilities, 12))
  r <- tryCatch({p <- parameter_inference(f); subset(p, parameter == "transition_probability")[, c("outcome","term","estimate","standard_error","conf_low","conf_high")]}, error = function(e) paste("refused", class(e)[1]))
  print(r)
}
# same with the general engine (occasion transitions, T=2 equals homogeneous)
set.seed(1); sep <- 12
s2 <- ifelse(runif(n) < .4, 2L, 1L)
d <- rbind(data.frame(id = 1:n, time = 1, y1 = rnorm(n, 0), y2 = rnorm(n, 0)),
           data.frame(id = 1:n, time = 2, y1 = rnorm(n, c(0, sep)[s2]), y2 = rnorm(n, c(0, sep)[s2])))
w <- character()
g <- withCallingHandlers(lta(d, c("y1","y2"), "id", 2, "time", transitions = "occasion", n_starts = 2, seed = 1),
  warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
cat("\nGeneral engine sep 12: n_par", g$n_parameters, " warnings:", paste(unique(w), collapse=", "), "\n")
print(get_results(g, "transitions"))
print(g$empty_transition_rows)
print(summary(colSums(g$expectation$moments[[1]]$pairs[[1]][,,1])))
