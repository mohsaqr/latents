suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
gen <- function(n, slope, seed, z_fun = rnorm) {
  set.seed(seed); T <- 3
  pid <- rep(seq_len(n), each = T); occ <- rep(seq_len(T), n)
  z <- z_fun(n * T)
  s <- matrix(0L, n, T); s[, 1] <- sample(2, n, TRUE)
  for (t in 2:T) { p <- plogis(-1.5 + slope * z[occ == t]); s[, t] <- ifelse(runif(n) < p, 3L - s[, t-1], s[, t-1]) }
  data.frame(id = pid, time = occ, z = z, y1 = rnorm(n * T, c(0, 3)[as.vector(t(s))]), y2 = rnorm(n * T, c(0, 3)[as.vector(t(s))]))
}
report <- function(label, d) {
  w <- character()
  f <- withCallingHandlers(lta(d, c("y1","y2"), "id", 2, "time", transition_covariates = "z", n_starts = 3, seed = 1),
    warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
  lt <- lapply(f$designs$transition, function(D) exp(.lta_log_transition(D, array(.lta_parameters(f)$transition[,,,1], dim(f$transition_coefficients)[1:3]))))
  smallest <- min(unlist(lt))
  ps <- round(c(f$transition_coefficients), 3)
  cat(sprintf("\n%s: warnings [%s]; boundary %s; smallest row-level move probability %.2e; coefficients %s\n", label, paste(unique(w), collapse=","), f$boundary, smallest, paste(ps, collapse = " ")))
  print(get_results(f, "transitions")[, c("occasion","from","to","probability")][1:4, ])
  r <- tryCatch(subset(parameter_inference(f), term == "z")[, c("from","to","estimate","standard_error")], error = function(e) paste("refused:", class(e)[1]))
  print(r); invisible(f)
}
d <- gen(300, 2, 11)
report("no extreme row, slope 2 per SD", d)
dx <- d; dest <- which(dx$time == 2)[3]; dx$z[dest] <- -6.5
report("one destination row at -6.5 SD", dx)
report("slope 3 per SD, n = 300", gen(300, 3, 12))
report("skewed covariate, exp(z), slope 1.2", gen(300, 1.2, 13, function(k) scale(rexp(k))[, 1]))

cat("\n--- trailing-dot term\n")
sim <- gen(150, 1, 3)
sim$title <- factor(ifelse(sim$id %% 2 == 0, "Aa", "Zed."))
fd <- suppressWarnings(lta(sim, c("y1","y2"), "id", 2, "time", transition_covariates = "title", n_starts = 2, seed = 1))
print(dimnames(fd$transition_coefficients)[[1]])
print(.lta_coefficient_table(fd, "transition", .95, "observed")[, c("from","to","term","estimate","standard_error","estimated")])
cat("identified coordinates", sum(.lta_identified_coordinates(fd)), "; n_parameters", fd$n_parameters, "\n")
report("slope 3 per SD, n = 1000, no outlier", gen(1000, 3, 21))
