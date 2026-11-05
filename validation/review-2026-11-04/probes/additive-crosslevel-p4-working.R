suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
vars <- c("y1", "y2")
# No group structure at all: every member draws her profile iid, whatever the group.
draw <- function(seed, J = 100, n = 6) { set.seed(seed)
  k <- sample(1:2, J * n, TRUE, c(.5, .5)); pm <- rbind(c(-1.5, -1), c(1.5, 1))
  data.frame(team = rep(sprintf("t%03d", seq_len(J)), each = n), y1 = rnorm(J * n, pm[k, 1], .6), y2 = rnorm(J * n, pm[k, 2], .6)) }
crit <- function(f) c(LL = f$log_likelihood, q = f$n_parameters, bic = f$bic, bic_ind = f$bic_individual)
rows <- lapply(1:6, function(seed) {
  D <- draw(seed)
  one <- function(fam, H) suppressWarnings(multilpa(D, vars, "team", n_profiles = 2, n_group_classes = H, family = fam, seed = 1, n_starts = 4, variance_model = "equal", between_variance = "equal"))
  f1 <- one("full_cross_level", 1); f2 <- one("full_cross_level", 2)
  r2 <- one("restricted_cross_level", 2)
  p1 <- suppressWarnings(multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 1, variance_model = "equal", seed = 1, n_starts = 4))
  p2 <- suppressWarnings(multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 2, variance_model = "equal", seed = 1, n_starts = 4))
  data.frame(seed = seed,
    full_gain = f2$log_likelihood - f1$log_likelihood, rest_gain = r2$log_likelihood - f1$log_likelihood,
    prof_gain = p2$log_likelihood - p1$log_likelihood,
    full_dBIC_H2_minus_H1 = f2$bic - f1$bic, rest_dBIC_H2_minus_H1 = r2$bic - f1$bic, prof_dBIC_H2_minus_H1 = p2$bic - p1$bic,
    full_dBICind = f2$bic_individual - f1$bic_individual,
    full_comp_c1 = min(f2$composition[, 1]), full_comp_c2 = max(f2$composition[, 1]),
    prof_comp_c1 = min(p2$profile_probabilities[, 1]), prof_comp_c2 = max(p2$profile_probabilities[, 1]))
})
print(do.call(rbind, rows), digits = 4)
# Decomposing one case: how much of the full gain is the rating part and how much the group-mean part?
D <- draw(1)
f2 <- suppressWarnings(multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1, n_starts = 4, variance_model = "equal", between_variance = "equal"))
print(f2$composition); print(f2$group_means); print(f2$group_mean_variances); print(f2$group_probabilities)
# Terminal loss of a sensitivity seed: what the caller sees
G <- draw(2)
sfit <- suppressMessages(lpa(G, vars, n_profiles = 2, n_starts = 1, seed = 1))
lost <- function(indices, collected, cores, retry = FALSE) lapply(seq_along(indices), function(i) if (i == 2L || retry) NULL else collected(indices[[i]]))
print(tryCatch(testthat::with_mocked_bindings(
  suppressWarnings(sensitivity(sfit, seeds = 1:3, n_starts = 1, cores = 2)),
  .latents_fork_pass = lost, .latents_socket_pass = lost), error = function(e) paste("ERROR", class(e)[1], ":", conditionMessage(e))))
# doc claim: dispersion fit raises latents_weak_class
set.seed(9); J <- 60; X <- do.call(rbind, lapply(seq_len(J), function(j) { b <- rnorm(2, 0, .5); data.frame(team = paste0("t", j), y1 = rnorm(5, b[1], 1), y2 = rnorm(5, b[2], 1)) }))
w <- character(); invisible(withCallingHandlers(tryCatch(summary(multilpa(X, vars, "team", n_group_classes = 2, family = "dispersion", seed = 1, n_starts = 2)), error = function(e) cat("err", class(e)[1], "\n")), warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") }))
print(unique(w))
# Effect of group size on the prevalence column, moderate imbalance
pm <- rbind(c(-1.5, -1), c(1.5, 1))
set.seed(12); type <- rep(1:2, 40); M <- do.call(rbind, lapply(seq_along(type), function(j) { n <- c(5L, 15L)[type[j]]; k <- sample(1:2, n, TRUE, if (type[j] == 1) c(.2, .8) else c(.75, .25)); data.frame(team = sprintf("t%02d", j), y1 = rnorm(n, pm[k, 1], .6), y2 = rnorm(n, pm[k, 2], .6)) }))
f <- multilpa(M, vars, "team", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1, n_starts = 3)
r <- multilpa(M, vars, "team", n_profiles = 2, n_group_classes = 2, family = "restricted_cross_level", seed = 1, n_starts = 3)
print(subset(get_results(f, "profiles"), level == "individual" & indicator == "y1")[, c("profile", "mean", "prevalence")])
print(subset(get_results(r, "profiles"), level == "individual" & indicator == "y1")[, c("profile", "mean", "prevalence")])
cat("full: share of individuals (posterior)", colMeans(f$subject_posteriors), "\n"); print(get_results(f, "group_classes"))
