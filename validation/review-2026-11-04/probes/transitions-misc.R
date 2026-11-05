suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(n, T, seed) {
  set.seed(seed)
  P <- matrix(c(.8,.2,.3,.7), 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    s <- integer(T); s[1] <- sample(2, 1, prob = c(.6,.4))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 3)[s]), y2 = rnorm(T, c(1, -2)[s], .8), y3 = rnorm(T, c(0, 2)[s]), z = rnorm(T), w = rep(rnorm(1), T))
  }))
}
vars <- c("y1","y2","y3")
d <- sim(250, 4, 1)
cat("--- order 2 + occasion: second_order_transitions has one row per pair x destination?\n")
f <- suppressWarnings(lta(d, vars, "id", 2, "time", order = 2, transitions = "occasion", n_starts = 3, seed = 1))
print(names(get_results(f, "second_order_transitions"))); print(nrow(get_results(f, "second_order_transitions")))
print(get_results(f, "second_order_transitions"))
print(subset(get_results(f, "second_order_coefficients"), previous == "profile_1" & from == "profile_1")[, c("to","term","estimate")])
print(get_results(f, "transitions"))
f0 <- suppressWarnings(lta(d, vars, "id", 2, "time", order = 2, n_starts = 3, seed = 1))
print(get_results(f0, "transitions"))

cat("\n--- Score at the reported optimum, EM-only paths (structures) and finish paths\n")
opt <- function(label, ...) {
  dd <- d; dd$u <- rep(c(1, 2, 3), length.out = 250)[dd$id]
  dm <- dd; dm$y1[seq(1, nrow(dm), 9)] <- NA; dm <- dm[!(dm$id %% 6 == 0 & dm$time == 4), ]
  args <- list(...)
  data <- if (isTRUE(args$fiml_gap)) dm else dd
  args$fiml_gap <- NULL
  if (identical(args$missing, "fiml") && !identical(data, dm)) data <- dm
  f <- suppressWarnings(do.call(lta, c(list(data, vars, "id", n_profiles = 2, time = "time", n_starts = 2, seed = 1), args)))
  s <- colSums(.lta_group_scores(.lta_theta(f), f))
  id <- .lta_identified_coordinates(f)
  cat(sprintf("%-40s conv %s  LL %.5f  max|score| %.2e  iter %s\n", label, f$converged, f$log_likelihood, max(abs(s[id])), f$iterations))
  invisible(f)
}
opt("diagonal occasion measurement", measurement = "occasion")
opt("equal variance, occasion measurement", measurement = "occasion", variance_model = "equal")
opt("EEE occasion measurement", measurement = "occasion", model = "EEE")
opt("EEI occasion measurement fiml", measurement = "occasion", model = "EEI", missing = "fiml")
opt("VEI occasion transitions", transitions = "occasion", model = "VEI")
opt("VEE fiml", transitions = "occasion", model = "VEE", missing = "fiml")
opt("EEE weights", weights = "u", model = "EEE")
opt("EVV weights occasion measurement", weights = "u", model = "EVV", measurement = "occasion")
opt("EEI weights fiml", weights = "u", model = "EEI", missing = "fiml")
opt("VVI stayer", model = "VVI", mover_stayer = TRUE)
cat("\n--- homogeneous vs general on structures (same data, general via weights of 1)\n")
d$one <- 1
invisible(lapply(c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","VVE","EEV","VEV","EVV","VVV"), function(m) {
  a <- suppressWarnings(lta(d, vars, "id", 2, "time", model = m, n_starts = 3, seed = 1))
  b <- suppressWarnings(lta(d, vars, "id", 2, "time", model = m, n_starts = 3, seed = 1, weights = "one"))
  cat(sprintf("%s LL hom %.6f gen %.6f diff %.2e npar %d %d\n", m, a$log_likelihood, b$log_likelihood, a$log_likelihood - b$log_likelihood, a$n_parameters, b$n_parameters))
}))
cat("\n--- second order + occasion + covariates + H=2: SEs vs numeric Hessian of the likelihood on raw coordinates\n")
fg <- suppressWarnings(lta(sim(300, 4, 3), vars, "id", 2, "time", order = 2, transitions = "occasion", transition_covariates = "z", initial_covariates = "w", n_starts = 2, seed = 1))
theta <- .lta_theta(fg); ll <- function(v) .lta_expectation(fg$x, fg$codes, fg$layout, fg$designs, .lta_unpack(v, fg), fg$occasion_of_row, NULL, NULL)$log_likelihood
g <- function(v) colSums(.lta_group_scores(v, fg))
cat("boundary", fg$boundary, " conv", fg$converged, " np", length(theta), "\n")
if (!fg$boundary) {
  hess <- optimHess(theta, function(v) -ll(v), function(v) -g(v))
  se_num <- sqrt(diag(solve(hess)))
  se_pkg <- sqrt(diag(vcov(fg)))
  cat("max rel SE diff", max(abs(se_pkg / se_num - 1)), "\n")
} else print("boundary fit: see below")
fi <- suppressWarnings(lta(sim(300, 4, 3), vars, "id", 2, "time", transitions = "occasion", transition_covariates = "z", initial_covariates = "w", n_starts = 2, seed = 1))
theta <- .lta_theta(fi); g <- function(v) colSums(.lta_group_scores(v, fi))
hess <- optimHess(theta, function(v) -.lta_expectation(fi$x, fi$codes, fi$layout, fi$designs, .lta_unpack(v, fi), fi$occasion_of_row, NULL, NULL)$log_likelihood, function(v) -g(v))
cat("first order occasion + cov: max rel SE diff", max(abs(sqrt(diag(vcov(fi))) / sqrt(diag(solve(hess))) - 1)), " boundary", fi$boundary, "\n")
pi <- parameter_inference(fi, level = .9)
cat("CI check: (conf_high-estimate)/SE =", unique(round((pi$conf_high - pi$estimate)/pi$standard_error, 6)), " qnorm(.95) =", round(qnorm(.95), 6), "\n")
cat("two-sided p check:", max(abs(pi$p_value - 2*pnorm(-abs(pi$estimate/pi$standard_error)))), "\n")

cat("\n--- bootstrap crash reproduction\n")
try_boot <- function(label, f, data) {
  r <- tryCatch(withCallingHandlers({ b <- parameter_inference(f, data = data, method = "bootstrap", iter = 15, n_starts = 1, seed = 1)
     sprintf("completed, n_valid %d of 15", attr(b, "n_valid")) }, warning = function(c) invokeRestart("muffleWarning")),
     error = function(e) sprintf("ERROR class [%s]: %s", paste(class(e), collapse = ","), conditionMessage(e)))
  cat(label, "->", r, "\n")
}
sim2 <- function(n, T, seed) { x <- sim(n, T, seed); x }
d3 <- sim2(60, 3, 6); extra <- transform(subset(sim2(1, 5, 7), time > 3), id = 1L); d3 <- rbind(d3, extra)
o <- suppressWarnings(lta(d3, vars, "id", 2, "time", transitions = "occasion", n_starts = 2, seed = 1)); try_boot("occasion transitions, one longer person", o, d3)
m <- suppressWarnings(lta(d3, vars, "id", 2, "time", measurement = "occasion", n_starts = 2, seed = 1)); try_boot("occasion measurement, one longer person", m, d3)
hc <- suppressWarnings(lta(d3, vars, "id", 2, "time", weights = NULL, n_starts = 2, seed = 1)); try_boot("homogeneous, one longer person", hc, d3)
rare <- sim2(60, 3, 8); rare$f <- factor(ifelse(rare$id == 1, "rare", ifelse(rare$id %% 2 == 0, "a", "b")))
c1 <- suppressWarnings(lta(rare, vars, "id", 2, "time", transition_covariates = "f", n_starts = 2, seed = 1)); try_boot("factor covariate with a rare level", c1, rare)
