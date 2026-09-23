# Run from the package root: Rscript equivalence/houle-2026/check-likelihood.R
# Synthetic structural check, NOT a reproduction of the published data analysis.
pkgload::load_all(quiet = TRUE)

# Measurement values printed in the supplement's dispersion-heterogeneity
# second stage. Mixing weights and observations below are independently chosen.
measurement <- list(
  means = matrix(c(-2.49658, 2.49721, 1.92801, 1.54718, -.01351, .00624,
                   -1.72604, -2.05703, 2.48742, -2.50084), 5, 2, byrow = TRUE),
  variances = matrix(c(.08409, .08564, .64117, .10539, .37173, .36234,
                       .28702, .62374, .33606, .33277), 5, 2, byrow = TRUE),
  group_probabilities = c(.4, .6),
  profile_probabilities = rbind(c(.35, .25, .15, .15, .10),
                                c(.10, .15, .20, .25, .30)))
observations <- data.frame(
  group = rep(seq_len(3), each = 3),
  burn = c(-2.4, 1.8, 0.1, -1.5, 2.3, .2, -2.2, 1.6, 2.7),
  eng = c(2.6, 1.4, -.2, -2.2, -2.4, .3, 2.4, 1.7, -2.6))
str(observations)
stopifnot(!anyNA(observations), !anyDuplicated(observations))

scored <- multilpa(observations, c("burn", "eng"), "group",
                  n_profiles = 5, n_group_classes = 2, start = measurement,
                  n_starts = 1, max_iter = 0)

# Enumerate all 2 * 5^3 latent assignments within each group. This explicitly
# sums the joint density instead of using the package's nested E-step.
assignments <- as.matrix(expand.grid(h = 1:2, k1 = 1:5, k2 = 1:5, k3 = 1:5))
oracle <- lapply(split(seq_len(nrow(observations)), observations$group), function(rows) {
  stopifnot(length(rows) == 3L)
  weights <- apply(assignments, 1L, function(a) {
    stopifnot(length(a) == 4L)
    k <- a[2:4]
    measurement$group_probabilities[a[1]] *
      prod(measurement$profile_probabilities[cbind(a[1], k)]) *
      prod(dnorm(observations$burn[rows], measurement$means[k, 1],
                 sqrt(measurement$variances[k, 1]))) *
      prod(dnorm(observations$eng[rows], measurement$means[k, 2],
                 sqrt(measurement$variances[k, 2])))
  })
  posterior <- weights / sum(weights)
  list(log_likelihood = log(sum(weights)),
       group = vapply(1:2, function(h) sum(posterior[assignments[, 1] == h]), numeric(1)),
       subject = t(vapply(2:4, function(i) {
         vapply(1:5, function(k) sum(posterior[assignments[, i] == k]), numeric(1))
       }, numeric(5))))
})
expected_ll <- sum(vapply(oracle, `[[`, numeric(1), "log_likelihood"))
expected_group <- do.call(rbind, lapply(oracle, `[[`, "group"))
expected_subject <- do.call(rbind, lapply(oracle, `[[`, "subject"))
# The package still counts freely estimated measurement at max_iter = 0:
# 20 measurement + 8 conditional profile weights + 1 group weight = 29.
# A genuinely fixed-measurement second stage would count only 9 free weights.
errors <- c(log_likelihood = abs(scored$log_likelihood - expected_ll),
            group_posterior = max(abs(scored$group_posteriors - expected_group)),
            subject_posterior = max(abs(scored$subject_posteriors - expected_subject)))
stopifnot(all(errors < 1e-12),
          max(abs(scored$means - measurement$means)) < 1e-12,
          max(abs(scored$variances - measurement$variances)) < 1e-12,
          scored$n_parameters == 29)

# A start is not a fixed-measurement constraint: even one EM iteration updates
# the supplied measurement. These nine observations do not establish a fit.
updated <- multilpa(observations, c("burn", "eng"), "group",
                   n_profiles = 5, n_group_classes = 2, start = measurement,
                   n_starts = 1, max_iter = 1)
mean_change <- max(abs(updated$means - measurement$means))
stopifnot(mean_change > 1e-3,
          updated$log_likelihood >= scored$log_likelihood - 1e-10)
print(errors)
cat(sprintf("Maximum mean change after one EM iteration: %.8f\n", mean_change))
cat("PASS: exhaustive likelihood/posteriors; supplied starts update during EM.\n")
