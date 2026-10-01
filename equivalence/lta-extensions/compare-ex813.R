# Mplus User's Guide example 8.13 (Mplus 8.8): LTA with a known group cg
# (from g) affecting the initial classes (c1 ON cg) and with group-specific
# transitions (c2 ON c1 within each cg; c2 ON cg). Represented here by g as an
# initial and a transition covariate: per origin, an intercept and a g effect,
# which is exactly one free transition row per group. Mplus's likelihood also
# includes the known-class proportion (one parameter), added analytically.
# Reference: JStats validation/r-reference/lta-mplus-ug-ex8.13.out.txt
source(file.path("equivalence", "lta-extensions", "common.R"))
raw <- read.table(file.path(jstats, "validation", "r-reference", "mplus-ug-ex8.13.dat"))
names(raw) <- c(paste0("u1", 1:5), paste0("u2", 1:5), "g", "true_cg", "true1", "true2")
long <- do.call(rbind, lapply(1:2, function(t) {
  items <- raw[, paste0("u", t, 1:5)]
  names(items) <- paste0("u", 1:5)
  cbind(person = seq_len(nrow(raw)), occasion = t, g = raw$g, items)
}))
long[paste0("u", 1:5)] <- lapply(long[paste0("u", 1:5)], factor)
fit <- lta(long, paste0("u", 1:5), "person", n_profiles = 3, time = "occasion",
           categorical = paste0("u", 1:5), transition_covariates = "g",
           initial_covariates = "g", n_starts = 20, seed = 1, tol = 1e-10,
           max_iter = 2000)
share <- mean(raw$g)
known_class <- sum(raw$g) * log(share) + sum(1 - raw$g) * log(1 - share)
mplus_thresholds <- rbind(c(-1.253, -0.895, -0.974, -1.039, -0.838),
                          c(1.061, 0.959, 0.860, 0.829, 0.941),
                          c(0.960, 0.719, -1.078, -0.824, -1.286))
ours_thresholds <- vapply(fit$measurement[[1L]]$response_probabilities,
                          function(p) log(p[, 1L] / p[, 2L]), numeric(3))
perm <- align(ours_thresholds, mplus_thresholds)
# C1 ON CG#1 (cg#1 is g = 0) against C1#3: minus our g effect, re-referenced.
effect <- c(fit$initial_coefficients["g", , 1L], 0)[perm]
initial_cg <- -(effect[1:2] - effect[3])
result <- rbind(
  row_result("ex8.13", "log likelihood (+ known-class term)",
             fit$log_likelihood + known_class, -14596.240, 5e-3),
  row_result("ex8.13", "parameter count (+ known-class proportion)",
             fit$n_parameters + 1L, 32, 0),
  row_result("ex8.13", "max |threshold diff|",
             max(abs(ours_thresholds[perm, ] - mplus_thresholds)), 0, 2e-3),
  row_result("ex8.13", "max |C1 ON CG#1 diff|",
             max(abs(initial_cg - c(0.908, 1.948))), 0, 2e-3))
print(result, row.names = FALSE, digits = 8)
saveRDS(result, file.path("equivalence", "lta-extensions", "result-ex813.rds"))
