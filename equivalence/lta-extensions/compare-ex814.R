# Mplus User's Guide example 8.14 (Mplus 8.8): LTA with a covariate on the
# initial classes (c1 ON x) and origin-specific covariate effects on the
# transitions (c2 ON x within each c1 class); binary items, thresholds equal
# over time. Reference: JStats validation/r-reference/mplus-ug-ex8.14-full.out.txt
source(file.path("equivalence", "lta-extensions", "common.R"))
raw <- read.table(file.path(jstats, "validation", "r-reference", "mplus-ug-ex8.14.dat"))
names(raw) <- c(paste0("u1", 1:5), paste0("u2", 1:5), "x", "true1", "true2")
long <- do.call(rbind, lapply(1:2, function(t) {
  items <- raw[, paste0("u", t, 1:5)]
  names(items) <- paste0("u", 1:5)
  cbind(person = seq_len(nrow(raw)), occasion = t, x = raw$x, items)
}))
long[paste0("u", 1:5)] <- lapply(long[paste0("u", 1:5)], factor)
fit <- lta(long, paste0("u", 1:5), "person", n_profiles = 3, time = "occasion",
           categorical = paste0("u", 1:5), transition_covariates = "x",
           initial_covariates = "x", n_starts = 10, seed = 1, tol = 1e-10,
           max_iter = 2000)

mplus_ll <- -12902.140
mplus_thresholds <- rbind(c(-1.107, -1.442, -1.148, -0.932, -1.196),
                          c(0.953, 0.934, -0.865, -1.088, -0.794),
                          c(0.858, 1.082, 1.125, 1.012, 1.155))
mplus_initial_x <- c(-1.552, -2.061)          # C1#1, C1#2 ON X (reference C1#3)
mplus_odds_x <- rbind(c(1, 0.361, 1.991), c(2.003, 1, 10.520), c(1.685, 3.728, 1))

response <- fit$measurement[[1L]]$response_probabilities
ours_thresholds <- vapply(response, function(p) log(p[, 1L] / p[, 2L]), numeric(3))
perm <- align(ours_thresholds, mplus_thresholds)
# Initial logits on x against the profile aligned to Mplus's reference class.
initial_x <- c(fit$initial_coefficients["x", , 1L], 0)[perm]
initial_x <- initial_x[1:2] - initial_x[3]
# Transition odds ratios for x against staying, aligned.
odds <- diag(3)
invisible(lapply(1:3, function(k) {
  others <- setdiff(1:3, k)
  odds[k, others] <<- exp(fit$transition_coefficients["x", , k, 1L])
}))
odds <- odds[perm, perm]
result <- rbind(
  row_result("ex8.14", "log likelihood", fit$log_likelihood, mplus_ll, 5e-3),
  row_result("ex8.14", "parameter count", fit$n_parameters, 31, 0),
  row_result("ex8.14", "max |threshold diff|", max(abs(ours_thresholds[perm, ] - mplus_thresholds)), 0, 2e-3),
  row_result("ex8.14", "max |initial logit on x diff|", max(abs(initial_x - mplus_initial_x)), 0, 2e-3),
  row_result("ex8.14", "max |log odds ratio diff|", max(abs(log(odds) - log(mplus_odds_x))), 0, 5e-3))
print(result, row.names = FALSE, digits = 8)
saveRDS(result, file.path("equivalence", "lta-extensions", "result-ex814.rds"))
