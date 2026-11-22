# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Wald inference for latent transition models. The analytic score is checked
# against numerical differentiation of the forward-backward likelihood, and the
# standard errors against a Hessian of that likelihood alone, so neither is
# compared with the code under test.

skip_on_cran()

.lta_inference_fixture <- function(n = 120L, waves = 5L, seed = 1L,
                                   stay = c(0.85, 0.6), categorical = FALSE,
                                   missing_share = 0, ragged = FALSE) {
  set.seed(seed)
  group_class <- stats::rbinom(n, 1L, 0.5) + 1L
  frame <- do.call(rbind, lapply(seq_len(n), function(person) {
    path <- Reduce(function(previous, wave) {
      if (stats::runif(1L) < stay[group_class[person]]) previous else 3L - previous
    }, seq_len(waves - 1L), init = sample(1:2, 1L), accumulate = TRUE)
    data.frame(person = person, wave = seq_len(waves), state = path)
  }))
  frame$a <- stats::rnorm(nrow(frame), c(0, 2)[frame$state])
  frame$b <- stats::rnorm(nrow(frame), c(0, 1.8)[frame$state]) + 0.3 * frame$a
  if (categorical) {
    frame$q <- ifelse(stats::runif(nrow(frame)) < c(0.8, 0.25)[frame$state],
                      "yes", "no")
  }
  if (missing_share > 0) {
    frame$a[stats::runif(nrow(frame)) < missing_share] <- NA
  }
  if (ragged) frame <- frame[!(frame$wave == waves & frame$person %% 3L == 0L), ]
  frame$state <- NULL
  frame
}

.lta_inference_fit <- function(data, vars = c("a", "b"), ...) {
  quietly(lta(data, vars, "person", 2L, time = "wave", n_starts = 3, seed = 1,
              tol = 1e-9, max_iter = 20000, ...))
}

# The likelihood and the analytic score, on the coordinates inference uses.
.lta_inference_pieces <- function(fit) {
  x <- matrix(as.numeric(fit$indicator_data), nrow(fit$indicator_data))
  centers <- colMeans(x, na.rm = TRUE)
  centred <- sweep(x, 2L, centers, "-")
  layout <- latents:::.multilpa_sequence_layout(fit$group_index, fit$time_values,
                                      fit$n_groups, fit$occasions)
  view <- latents:::.multilpa_transition_view(fit, centers)
  blocks <- latents:::.multilpa_transition_shared_blocks(view, "unconstrained")
  widths <- c(measurement = length(blocks$measurement), group = length(blocks$group))
  list(
    theta = unname(latents:::.multilpa_transition_encode(fit, "unconstrained", centers)),
    negative = function(value) {
      -latents:::.multilpa_transition_expectation(
        centred, layout,
        latents:::.multilpa_transition_decode(value, view, widths[["measurement"]],
                                    widths[["group"]]),
        fit$categorical_data)$log_likelihood
    },
    score = function(value) {
      colSums(latents:::.multilpa_transition_group_scores(value, centred, layout, view,
                                                widths, fit$categorical_data))
    })
}

# The reported standard errors (vcov() on the unconstrained scale) against
# the inverse of numDeriv's Hessian of the likelihood alone. 1e-4 relative,
# the original bound: Richardson-extrapolated second differences of a
# likelihood of several hundred terms are accurate to about that.
.lta_check_errors <- function(fit, setting) {
  pieces <- .lta_inference_pieces(fit)
  # Richardson-extrapolated differences of the likelihood alone.
  # optimHess() is not accurate enough to be the reference: its own error on
  # these likelihoods is about 1e-3 relative at its default step and grows
  # with roundoff at a finer one.
  reference <- sqrt(diag(solve(numDeriv::hessian(pieces$negative, pieces$theta))))
  reported <- sqrt(diag(vcov(fit, scale = "unconstrained")))
  expect_reference(reported, reference, reference = "numDeriv",
                   covers = c("lta()", "vcov.multilpa_transitions()"),
                   quantity = sprintf("standard errors, unconstrained scale, %s", setting),
                   tolerance = 1e-4, measure = "relative",
                   note = "Reference: square roots of the diagonal of the inverse of numDeriv::hessian() of the package's likelihood (.multilpa_transition_expectation), at the reported estimates.")
}

test_that("standard errors match a Hessian of the likelihood alone", {
  skip_without_reference("numDeriv")
  # The fixtures are drawn in the original order: the first case, then the second.
  first <- .lta_inference_fixture()
  second <- .lta_inference_fixture(n = 250L, stay = c(0.95, 0.4), ragged = TRUE)
  .lta_check_errors(.lta_inference_fit(first), "one group class")
  .lta_check_errors(.lta_inference_fit(second, n_group_classes = 2L),
                    "two group classes, ragged sequences")
})
