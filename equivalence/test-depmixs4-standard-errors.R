# Transition-model standard errors against depmixS4's standardError()
# (Visser & Speekenbrink, 2010). depmixS4 is set to latents' estimates, where
# its likelihood equals latents' to about 1e-13, and differentiates it
# numerically (nlme::fdHess). That finite-difference Hessian is itself only
# accurate to about 1-3%: against a Richardson-extrapolated Hessian of the same
# likelihood latents agrees to 1e-6 and depmixS4 differs by 2.7% on the
# three-profile case. The comparison is therefore at 5%, which a wrong score
# or a wrong delta method would exceed by far; the tight check against
# Richardson differentiation lives in tests/testthat/test-transition-inference.R.

skip_if_not_installed("depmixS4")

depmix_case <- function(profiles, categorical, seed) {
  set.seed(seed)
  stay <- 0.75
  transition <- matrix((1 - stay) / (profiles - 1), profiles, profiles)
  diag(transition) <- stay
  frame <- do.call(rbind, lapply(seq_len(80L), function(group) {
    path <- unlist(Reduce(function(previous, step) {
      sample.int(profiles, 1L, prob = transition[previous, ])
    }, seq_len(5L), init = sample.int(profiles, 1L), accumulate = TRUE))
    data.frame(g = group, t = seq_len(6L), state = path)
  }))
  centres <- seq(-2, 2, length.out = profiles)
  frame$y1 <- stats::rnorm(nrow(frame), centres[frame$state], 0.8)
  frame$y2 <- stats::rnorm(nrow(frame), 0.75 * centres[frame$state], 0.9)
  if (categorical) {
    frame$c1 <- factor(vapply(frame$state, function(state) {
      sample(c("a", "b", "c"), 1L,
             prob = c(0.7, 0.2, 0.1)[1L + (seq_len(3L) + state - 2L) %% 3L])
    }, character(1)), levels = c("a", "b", "c"))
  }
  indicators <- c("y1", "y2", if (categorical) "c1")
  fit <- lta(frame, indicators, "g", n_profiles = profiles, time = "t",
             categorical = if (categorical) "c1" else character(),
             n_starts = 8, seed = seed, max_iter = 5000, tol = 1e-12)
  families <- c(list(stats::gaussian(), stats::gaussian()),
                if (categorical) list(depmixS4::multinomial("identity")))
  model <- depmixS4::depmix(
    lapply(indicators, function(name) stats::as.formula(paste(name, "~ 1"))),
    data = frame, nstates = profiles, family = families,
    ntimes = as.integer(table(frame$g)))
  responses <- unlist(lapply(seq_len(profiles), function(profile) {
    c(as.vector(rbind(fit$means[profile, ], sqrt(fit$variances[profile, ]))),
      if (categorical) fit$response_probabilities[["c1"]][profile, ])
  }), use.names = FALSE)
  model <- depmixS4::setpars(model, c(
    fit$initial_probabilities[1L, ],
    as.vector(t(matrix(fit$transition_probabilities[, , 1L], profiles, profiles))),
    responses))
  list(fit = fit, model = model, profiles = profiles)
}

test_that("transition standard errors agree with depmixS4's numerical ones", {
  invisible(lapply(list(list(2L, FALSE, 101L), list(3L, FALSE, 102L),
                        list(2L, TRUE, 104L)), function(case) {
    built <- do.call(depmix_case, unname(case))
    fit <- built$fit
    k <- built$profiles
    expect_lt(abs(as.numeric(depmixS4::logLik(built$model)) - fit$log_likelihood),
              1e-8)
    theirs <- depmixS4::standardError(built$model)$se
    ours <- parameter_inference(fit)
    lookup <- function(parameter, outcome, term) {
      ours$standard_error[ours$parameter == parameter & ours$outcome == outcome &
                            ours$term == term]
    }
    profile <- function(index) sprintf("profile_%d", index)
    initial <- vapply(seq_len(k), function(state) {
      lookup("initial_probability", profile(state), "group_class_1")
    }, numeric(1))
    moves <- unlist(lapply(seq_len(k), function(from) {
      vapply(seq_len(k), function(to) {
        lookup("transition_probability", profile(to),
               sprintf("group_class_1:%s", profile(from)))
      }, numeric(1))
    }))
    block <- (length(theirs) - k - k^2) / k
    means <- vapply(seq_len(k), function(state) lookup("mean", profile(state), "y1"),
                    numeric(1))
    # depmixS4 reports the standard deviation; SE(variance) = 2 sd SE(sd).
    variances <- vapply(seq_len(k), function(state) {
      lookup("variance", profile(state), "y1")
    }, numeric(1))
    their_means <- theirs[k + k^2 + (seq_len(k) - 1L) * block + 1L]
    their_variances <- 2 * sqrt(fit$variances[, "y1"]) *
      theirs[k + k^2 + (seq_len(k) - 1L) * block + 2L]
    relative <- function(left, right) max(abs(left - right) / right)
    expect_lt(relative(initial, theirs[seq_len(k)]), 0.05)
    expect_lt(relative(moves, theirs[k + seq_len(k^2)]), 0.05)
    expect_lt(relative(means, their_means), 0.05)
    expect_lt(relative(variances, their_variances), 0.05)
  }))
})
