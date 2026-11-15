# Standard errors of the spherical structures (EII, VII) in the homogeneous
# lta() against depmixS4 (Visser & Speekenbrink, 2010), the one external
# transition-model implementation that can fit them: a Gaussian response per
# indicator with the standard deviations tied by linear equality constraints
# (all of them for EII, within each state for VII). The other eight
# constrained structures tie volume, shape or orientation nonlinearly and no
# installed package fits them in a transition model; they are checked against
# the independent reference and the extended engine in
# tests/equivalence/test-lta-structure-inference.R.
#
# depmixS4 is set to latents' estimates, where its likelihood must equal
# latents' (1e-8). Its likelihood is then differentiated by numDeriv's
# Richardson extrapolation in a free chart of the constrained model, and the
# standard errors carried to latents' natural parameters must agree within
# 1e-4 relative (criteria fixed before the comparison was run).
#
# depmixS4's own standardError() is compared too, within the 5% that
# equivalence/test-depmixs4-standard-errors.R allows it: it differentiates by
# simple finite differences, which is why it is not the primary reference. On
# an exploratory dataset (not retained) it was 8.7% high on one EII mean while
# numDeriv's Hessian of depmixS4's own likelihood reproduced latents to 1e-7;
# on the retained dataset it agrees within 0.1%.

skip_if_not_installed("depmixS4")
skip_if_not_installed("numDeriv")

depmix_structure_fixture <- readRDS(file.path(
  "..", "tests", "equivalence", "fixtures", "lta-structure-references.rds"))$cases$one_class

depmix_structure_case <- function(structure) {
  frame <- depmix_structure_fixture$data
  fit <- suppressWarnings(lta(frame, c("a", "b"), "person", 2L, time = "wave",
                              model = structure, n_starts = 4L, seed = 1L,
                              tol = 1e-12, max_iter = 20000L))
  k <- fit$n_profiles
  model <- depmixS4::depmix(list(a ~ 1, b ~ 1), data = frame, nstates = k,
                            family = list(stats::gaussian(), stats::gaussian()),
                            ntimes = as.integer(table(frame$person)))
  responses <- unlist(lapply(seq_len(k), function(profile) {
    as.vector(rbind(fit$means[profile, ], sqrt(fit$variances[profile, ])))
  }), use.names = FALSE)
  model <- depmixS4::setpars(model, c(
    fit$initial_probabilities[1L, ],
    as.vector(t(matrix(fit$transition_probabilities[, , 1L], k, k))), responses))
  base <- depmixS4::getpars(model)
  mean_at <- k + k^2 + c(1L, 3L, 5L, 7L)
  sd_at <- k + k^2 + c(2L, 4L, 6L, 8L)
  # The standard deviations' equality pattern, as depmixS4's `equal` takes it.
  equal <- rep(1, length(base))
  equal[sd_at] <- if (structure == "EII") 2 else c(2, 2, 3, 3)
  sd_of <- if (structure == "EII") c(1L, 1L, 1L, 1L) else c(1L, 1L, 2L, 2L)
  softmax <- function(logits) exp(c(logits, 0) - max(c(logits, 0))) /
    sum(exp(c(logits, 0) - max(c(logits, 0))))
  # A free chart of the constrained model: logits against the last profile,
  # the four means, and the log standard deviations that remain free.
  expand <- function(z) {
    p <- base
    p[seq_len(k)] <- softmax(z[1L])
    p[k + 1:2] <- softmax(z[2L])
    p[k + 3:4] <- softmax(z[3L])
    p[mean_at] <- z[4:7]
    p[sd_at] <- exp(z[7L + sd_of])
    p
  }
  start <- c(log(base[1L] / base[2L]), log(base[3L] / base[4L]),
             log(base[5L] / base[6L]), base[mean_at],
             log(base[sd_at][!duplicated(sd_of)]))
  loglik <- function(z) as.numeric(depmixS4::logLik(depmixS4::setpars(model, expand(z))))
  # latents' natural order: means, the spread block (shared once for EII),
  # the group-class share, initial and transition probabilities.
  natural <- function(z) {
    p <- expand(z)
    variances <- p[sd_at]^2
    c(p[mean_at], if (structure == "EII") variances[1:2] else variances, 1,
      p[seq_len(k)], p[k + seq_len(k^2)])
  }
  covariance <- solve(-numDeriv::hessian(loglik, start))
  map <- numDeriv::jacobian(natural, start)
  their <- depmixS4::standardError(model, equal = equal)$se
  list(fit = fit, model = model, loglik = loglik(start),
       richardson = sqrt(pmax(diag(map %*% covariance %*% t(map)), 0)),
       estimate = unname(natural(start)),
       # depmixS4's own standard errors of the same natural parameters (the
       # variance from the standard deviation: SE(v) = 2 sd SE(sd)).
       standard_error = c(their[mean_at],
                          if (structure == "EII") 2 * base[sd_at][1:2] * their[sd_at][1:2] else
                            2 * base[sd_at] * their[sd_at],
                          0, their[seq_len(k)], their[k + seq_len(k^2)]))
}

test_that("spherical-structure standard errors of lta() agree with depmixS4's likelihood", {
  invisible(lapply(c("EII", "VII"), function(structure) {
    built <- depmix_structure_case(structure)
    ours <- parameter_inference(built$fit)
    expect_lt(abs(built$loglik - built$fit$log_likelihood), 1e-8, label = structure)
    expect_equal(ours$estimate, built$estimate, tolerance = 1e-10, label = structure)
    informative <- built$richardson > 1e-8
    relative <- abs(ours$standard_error - built$richardson)[informative] /
      built$richardson[informative]
    expect_lt(max(relative), 1e-4, label = structure)
  }))
})

test_that("depmixS4's standardError() agrees with lta() within its finite-difference error", {
  invisible(lapply(c("EII", "VII"), function(structure) {
    built <- depmix_structure_case(structure)
    ours <- parameter_inference(built$fit)$standard_error
    informative <- built$richardson > 1e-8
    relative <- abs(built$standard_error - ours)[informative] / ours[informative]
    expect_lt(max(relative), 0.05, label = structure)
  }))
})
