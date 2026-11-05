suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
vars <- c("y1", "y2")
draw <- function(J = 60, seed = 5, big = FALSE) {
  set.seed(seed)
  type <- rep(1:2, length.out = J)
  pm <- rbind(c(-1.5, -1), c(1.5, 1))
  do.call(rbind, lapply(seq_along(type), function(j) {
    n <- if (big) c(4L, 40L)[type[j]] else sample(3:8, 1)
    s <- if (type[j] == 1L) c(.8, .2) else c(.25, .75)
    k <- sample(1:2, n, TRUE, s)
    data.frame(team = sprintf("t%02d", j), type = type[j],
               y1 = rnorm(n, pm[k, 1], .6), y2 = rnorm(n, pm[k, 2], .6))
  }))
}
D <- draw()
# Independent likelihood: L_j = sum_h w_h f(m_j|h) prod_i sum_k pi_{k|h} f(x_ij|k)
hand <- function(fit, D) {
  ids <- rownames(fit$group_posteriors)
  H <- length(fit$group_probabilities); K <- nrow(fit$profile_means)
  # restricted's stored composition is descriptive; pass prevalences explicitly
  pi_mat <- if (fit$family == "full_cross_level") fit$composition else
    matrix(fit$profile_prevalences, H, K, byrow = TRUE)
  vapply(ids, function(id) {
    g <- D[D$team == id, ]
    x <- as.matrix(g[, vars])
    m <- colMeans(x)
    comps <- vapply(seq_len(H), function(h) {
      lm <- sum(dnorm(m, fit$group_means[h, ], sqrt(fit$group_mean_variances[h, ]), log = TRUE))
      li <- sum(apply(x, 1, function(row) {
        v <- vapply(seq_len(K), function(k) log(pi_mat[h, k]) + sum(dnorm(row, fit$profile_means[k, ], sqrt(fit$profile_variances[k, ]), log = TRUE)), 0)
        mv <- max(v); mv + log(sum(exp(v - mv)))
      }))
      log(fit$group_probabilities[h]) + lm + li
    }, 0)
    mx <- max(comps); mx + log(sum(exp(comps - mx)))
  }, 0)
}
full <- multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1, n_starts = 3, tol = 1e-10)
rest <- multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 2, family = "restricted_cross_level", seed = 1, n_starts = 3, tol = 1e-10)
cat("full reported", full$log_likelihood, "hand", sum(hand(full, D)), "\n")
cat("rest reported", rest$log_likelihood, "hand", sum(hand(rest, D)), "\n")
print(full); print(rest)
# Counts by hand: K=2,H=2,d=2 varying/varying
cat("n_parameters full", full$n_parameters, "hand", 2*2+2*2 + 2*(2-1) + 2*2+2*2 + 1, "; rest", rest$n_parameters, "hand", 2*2+2*2 + 1 + 2*2+2*2 + 1, "\n")
cat("BIC: full", full$bic, "hand", -2*full$log_likelihood + full$n_parameters*log(60), " bic_ind", full$bic_individual, "hand", -2*full$log_likelihood + full$n_parameters*log(nrow(D)), "\n")
# Finite-difference stationarity of the hand likelihood at the full fit (free parameters)
hand_theta <- function(theta, fit0, D) {
  f <- fit0
  f$profile_means[] <- theta[1:4]; f$profile_variances[] <- exp(theta[5:8])
  f$group_means[] <- theta[9:12]; f$group_mean_variances[] <- exp(theta[13:16])
  w <- plogis(theta[17]); f$group_probabilities[] <- c(w, 1 - w)
  p <- plogis(theta[18:19]); f$composition[] <- rbind(c(p[1], 1 - p[1]), c(p[2], 1 - p[2]))
  sum(hand(f, D))
}
th0 <- c(as.vector(full$profile_means), log(as.vector(full$profile_variances)), as.vector(full$group_means), log(as.vector(full$group_mean_variances)), qlogis(full$group_probabilities[1]), qlogis(full$composition[, 1]))
grad <- vapply(seq_along(th0), function(i) { e <- replace(numeric(19), i, 1e-4); (hand_theta(th0 + e, full, D) - hand_theta(th0 - e, full, D)) / 2e-4 }, 0)
cat("full: max |d working LL / d theta| at the estimate:", max(abs(grad)), "\n")
# Composition
print(get_results(full, "composition")); print(get_results(rest, "composition"))
w <- rest$group_posteriors[rest$sufficient_statistics$index, ]; r <- rest$subject_posteriors
byhand <- t(w) %*% r / colSums(w)
cat("restricted descriptive composition diff by hand", max(abs(byhand - rest$composition)), "\n")
cat("restricted estimated prevalences (model-implied in every class):", rest$profile_prevalences, "\n")
cat("full sums", rowSums(full$composition), "rest sums", rowSums(rest$composition), "\n")
cat("full: descriptive composition by the same formula\n")
wf <- full$group_posteriors[full$sufficient_statistics$index, ]; print(t(wf) %*% full$subject_posteriors / colSums(wf)); print(full$composition)
# `prevalence` column: "a profile's share of the individuals" under group-size imbalance
B <- draw(J = 60, seed = 7, big = TRUE)
fb <- multilpa(B, vars, "team", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1, n_starts = 3)
rb <- multilpa(B, vars, "team", n_profiles = 2, n_group_classes = 2, family = "restricted_cross_level", seed = 1, n_starts = 3)
cat("unequal sizes (class 1 groups n=4, class 2 groups n=40). N =", nrow(B), "\n")
cat("full table prevalence:", unique(subset(get_results(fb, "profiles"), level == "individual")$prevalence), "\n")
cat("full: posterior share of individuals:", colMeans(fb$subject_posteriors), "\n")
print(fb$composition); print(fb$group_probabilities)
cat("restricted table prevalence:", unique(subset(get_results(rb, "profiles"), level == "individual")$prevalence), "; posterior share", colMeans(rb$subject_posteriors), "\n")
print(fb$profile_means)
# One group class: both coincide
o1 <- multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 1, family = "full_cross_level", seed = 1, n_starts = 2, tol = 1e-11)
o2 <- multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 1, family = "restricted_cross_level", seed = 1, n_starts = 2, tol = 1e-11)
cat("H=1 LL full", o1$log_likelihood, "restricted", o2$log_likelihood, "; params", o1$n_parameters, o2$n_parameters, "\n")
# Unit rescale of the working likelihood: shift (N + J) d log s
Ds <- D; Ds$y1 <- D$y1 * 10; Ds$y2 <- D$y2 * 10
fs <- multilpa(Ds, vars, "team", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1, n_starts = 3, tol = 1e-10)
cat("rescale shift full", full$log_likelihood - fs$log_likelihood, "theory (N+J)d log 10", (nrow(D) + 60) * 2 * log(10), "\n")
# Does the working-likelihood restriction hold? Which verbs accept a mix?
add <- multilpa(D, vars, "team", n_group_classes = 2, family = "additive", seed = 1, n_starts = 2)
prof <- multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 2, seed = 1, n_starts = 2)
show <- function(expr) tryCatch({ v <- expr; if (is.data.frame(v)) print(v) else print(v); invisible() }, error = function(e) cat("ERROR", class(e)[1], ":", substr(conditionMessage(e), 1, 110), "\n"))
cat("--- compare_models(full, rest)\n"); show(compare_models(full, rest))
cat("--- compare_models(full, prof)\n"); show(compare_models(full, prof))
cat("--- bootstrap_lrt(rest, full)\n"); show(bootstrap_lrt(rest, full, iter = 2, n_starts = 1, seed = 1))
cat("--- bootstrap_lrt(prof, full)\n"); show(bootstrap_lrt(prof, full, iter = 2, n_starts = 1, seed = 1))
cat("--- AIC(full, add, prof)\n"); show(AIC(full, add, prof)); show(BIC(full, add, prof))
cat("--- enumerate family mix\n"); show(enumerate_classes(D, vars, "team", family = c("additive", "full_cross_level"), n_group_classes = 2, n_starts = 1))
cat("--- plots/tables on fits\n"); show(get_best(full)); show(get_checks(full)); show(diagnostics(full))
cat("LL: full", full$log_likelihood, "additive", add$log_likelihood, "profiles", prof$log_likelihood, "\n")

# Weights vs duplicates in both cross-level forms
dupl <- D[D$team %in% sprintf("t%02d", 1:20), ]; dupl$team <- paste0(dupl$team, "_d")
D2 <- rbind(D, dupl); Dw <- D; Dw$w <- ifelse(Dw$team %in% sprintf("t%02d", 1:20), 2, 1)
for_family <- function(fam) {
  a <- multilpa(D2, vars, "team", n_profiles = 2, n_group_classes = 2, family = fam, seed = 1, n_starts = 3, tol = 1e-10)
  b <- multilpa(Dw, vars, "team", n_profiles = 2, n_group_classes = 2, family = fam, seed = 1, n_starts = 3, tol = 1e-10, weights = "w")
  align <- function(m1, m2) min(max(abs(m1 - m2)), max(abs(m1 - m2[nrow(m2):1, ])))
  cat(fam, ": profile means diff", align(a$profile_means, b$profile_means), "group means diff", align(a$group_means, b$group_means), "prevalence diff", min(max(abs(a$profile_prevalences - b$profile_prevalences)), max(abs(a$profile_prevalences - rev(b$profile_prevalences)))), "group prob diff", min(max(abs(a$group_probabilities - b$group_probabilities)), max(abs(a$group_probabilities - rev(b$group_probabilities)))), "\n")
  cat(fam, ": LL dup", a$log_likelihood, "LL weighted", b$log_likelihood, "weighted * 80/60", b$log_likelihood * 80 / 60, "\n")
  cat(fam, ": composition\n"); print(a$composition); print(b$composition)
}
for_family("full_cross_level"); for_family("restricted_cross_level")
# Zero-weight groups: the counts in fit table / print
Dz <- D; Dz$w <- 1; Dz$w[Dz$team %in% sprintf("t%02d", 1:5)] <- 0
nz <- sum(Dz$w == 0)
z <- multilpa(Dz, vars, "team", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1, n_starts = 2, weights = "w")
sub <- multilpa(subset(D, !team %in% sprintf("t%02d", 1:5)), vars, "team", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1, n_starts = 2)
cat("zero weights: LL", z$log_likelihood, "subset LL", sub$log_likelihood, "; bic", z$bic, sub$bic, "\n")
print(get_results(z, "fit")); print(get_results(sub, "fit")); print(z); cat("nobs", nobs(z), "\n")
zadd <- multilpa(Dz, vars, "team", n_group_classes = 2, family = "additive", seed = 1, n_starts = 2, weights = "w")
print(get_results(zadd, "fit")[, c("n_groups", "n_obs")])
zr <- multilpa(Dz, vars, "team", n_profiles = 2, n_group_classes = 2, family = "restricted_cross_level", seed = 1, n_starts = 2, weights = "w")
subr <- multilpa(subset(D, !team %in% sprintf("t%02d", 1:5)), vars, "team", n_profiles = 2, n_group_classes = 2, family = "restricted_cross_level", seed = 1, n_starts = 2)
cat("restricted zero weights: LL", zr$log_likelihood, "subset LL", subr$log_likelihood, "\n"); print(get_results(zr, "fit")[, c("n_groups", "n_obs", "bic")]); print(get_results(subr, "fit")[, c("n_groups", "n_obs", "bic")])
