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
  layout <- .multilpa_sequence_layout(fit$group_index, fit$time_values,
                                      fit$n_groups, fit$occasions)
  view <- .multilpa_transition_view(fit, centers)
  blocks <- .multilpa_transition_shared_blocks(view, "unconstrained")
  widths <- c(measurement = length(blocks$measurement), group = length(blocks$group))
  list(
    theta = unname(.multilpa_transition_encode(fit, "unconstrained", centers)),
    negative = function(value) {
      -.multilpa_transition_expectation(
        centred, layout,
        .multilpa_transition_decode(value, view, widths[["measurement"]],
                                    widths[["group"]]),
        fit$categorical_data)$log_likelihood
    },
    score = function(value) {
      colSums(.multilpa_transition_group_scores(value, centred, layout, view,
                                                widths, fit$categorical_data))
    })
}

.lta_check_score <- function(fit, label) {
  pieces <- .lta_inference_pieces(fit)
  set.seed(9)
  away <- pieces$theta + 0.05 * stats::rnorm(length(pieces$theta))
  numerical <- vapply(seq_along(away), function(index) {
    step <- replace(numeric(length(away)), index, 1e-5)
    -(pieces$negative(away + step) - pieces$negative(away - step)) / 2e-5
  }, numeric(1))
  expect_lt(max(abs(pieces$score(away) - numerical)), 1e-5, label = label)
}

test_that("the analytic score is the derivative of the forward-backward likelihood", {
  # Every block and layout the score has to handle: one and two group classes,
  # full and equal covariance, a categorical indicator, missing values and
  # ragged sequences on both occasion scales.
  .lta_check_score(.lta_inference_fit(.lta_inference_fixture()), "diagonal")
  two_classes <- .lta_inference_fixture(n = 250L, stay = c(0.95, 0.4))
  .lta_check_score(.lta_inference_fit(two_classes, n_group_classes = 2L),
                   "two group classes")
  .lta_check_score(.lta_inference_fit(two_classes, n_group_classes = 2L,
                                      variance_model = "equal"), "equal variances")
  .lta_check_score(.lta_inference_fit(.lta_inference_fixture(),
                                      covariance_model = "full"), "full covariance")
  .lta_check_score(.lta_inference_fit(.lta_inference_fixture(categorical = TRUE),
                                      vars = c("a", "b", "q"), categorical = "q"),
                   "categorical")
  .lta_check_score(.lta_inference_fit(.lta_inference_fixture(missing_share = 0.15),
                                      missing = "fiml", covariance_model = "full"),
                   "missing, full covariance")
  .lta_check_score(.lta_inference_fit(.lta_inference_fixture(ragged = TRUE),
                                      occasions = "grid"), "ragged, grid")
})

test_that("transition score indexing also holds with three profiles", {
  fit <- quietly(lta(.lta_inference_fixture(), c("a", "b"), "person", 3L,
                      n_group_classes = 2L, time = "wave", n_starts = 1,
                      seed = 4, max_iter = 20))
  # A derivative identity holds away from an optimum too. Three profiles
  # exercise distinct transition rows and columns beyond a binary swap.
  .lta_check_score(fit, "three profiles, two group classes")
})

test_that("coef(), vcov() and parameter_inference() name every parameter alike", {
  invisible(lapply(c("diagonal", "full"), function(covariance_model) {
    fit <- .lta_inference_fit(.lta_inference_fixture(),
                              covariance_model = covariance_model)
    inference <- parameter_inference(fit)
    expect_identical(names(coef(fit)), rownames(vcov(fit)))
    expect_identical(.multilpa_parameter_names(inference), names(coef(fit)))
    expect_equal(inference$estimate, unname(coef(fit)))
    expect_identical(rownames(confint(fit)), names(coef(fit)))
  }))
})

test_that("the probability rows are internally consistent", {
  fit <- .lta_inference_fit(.lta_inference_fixture(n = 250L, stay = c(0.95, 0.4)),
                            n_group_classes = 2L)
  inference <- parameter_inference(fit)
  transitions <- subset(inference, parameter == "transition_probability")
  # Each row of a two-profile transition matrix is one free probability, so
  # its two entries share one standard error.
  by_row <- split(transitions$standard_error, transitions$term)
  expect_true(all(vapply(by_row, function(errors) {
    isTRUE(all.equal(errors[1L], errors[2L]))
  }, logical(1))))
  expect_true(all(is.na(transitions$p_value)))
  # Estimates are the fitted transition probabilities, row by row.
  expect_equal(sort(transitions$estimate),
               sort(as.vector(fit$transition_probabilities)))
  robust <- parameter_inference(fit, vcov_type = "robust")
  expect_true(all(is.finite(robust$standard_error)))
})

test_that("inference refuses what it cannot answer, by class", {
  data <- .lta_inference_fixture()
  loose <- quietly(lta(data, c("a", "b"), "person", 2L, time = "wave",
                       n_starts = 1, seed = 1, max_iter = 3))
  expect_error(parameter_inference(loose), class = "latents_no_converge")
  fit <- .lta_inference_fit(data)
  # A profile never occupied before the last occasion leaves its row empty.
  empty <- fit
  empty$empty_transition_rows[1L, 1L] <- TRUE
  expect_error(parameter_inference(empty), class = "latents_boundary_fit")
  altered <- data
  altered$a[1L] <- altered$a[1L] + 1
  expect_error(parameter_inference(fit, altered),
               class = "latents_bad_inference_data")
})

# The conditional likelihood's numerical Hessian over the free logits, with
# the held ones at their fitted values: independent of the analytic score and
# of the code that chooses and drops the held coordinates.
.lta_conditional_errors <- function(fit, held) {
  x <- matrix(as.numeric(fit$indicator_data), nrow(fit$indicator_data))
  centers <- colMeans(x, na.rm = TRUE)
  centred <- sweep(x, 2L, centers, "-")
  layout <- .multilpa_sequence_layout(fit$group_index, fit$time_values,
                                      fit$n_groups, fit$occasions)
  view <- .multilpa_transition_view(fit, centers)
  blocks <- .multilpa_transition_shared_blocks(view, "unconstrained")
  theta <- .multilpa_transition_encode(fit, "unconstrained", centers)
  free <- !names(theta) %in% held
  negative <- function(value) {
    -.multilpa_transition_expectation(
      centred, layout,
      .multilpa_transition_decode(replace(unname(theta), free, value), view,
                                  length(blocks$measurement), length(blocks$group)),
      fit$categorical_data)$log_likelihood
  }
  sqrt(diag(solve(numDeriv::hessian(negative, unname(theta[free])))))
}

test_that("`fix` holds a bound initial probability and infers the rest conditionally", {
  skip_if_not_installed("numDeriv")
  # The LTA article's model: one group class starts nobody in profile 1.
  activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  moves <- quietly(lta(course_engagement, activity, "student", n_profiles = 3,
                       n_group_classes = 2, time = "sequence", n_starts = 3,
                       seed = 1), c(.multilpa_expected_warnings, "latents_boundary"))
  expect_true(moves$boundary)
  expect_error(parameter_inference(moves), class = "latents_boundary_fit")
  fixed <- parameter_inference(moves, boundary = "fix")
  expect_identical(attr(fixed, "fixed_at_bound"),
                   "profile.initial_probability.profile_1.group_class_1")
  expect_identical(nrow(fixed), length(coef(moves)) - 1L)
  expect_true(all(is.finite(fixed$standard_error) & fixed$standard_error > 0))
  errors <- .lta_conditional_errors(moves, "initial_logit[1,1]")
  expect_equal(sqrt(diag(attr(fixed, "covariance_unconstrained"))), errors,
               tolerance = 1e-4, ignore_attr = TRUE)
  # vcov() and confint() answer from the same conditional covariance.
  expect_equal(sqrt(diag(vcov(moves, boundary = "fix"))), fixed$standard_error,
               ignore_attr = TRUE)
  expect_identical(nrow(confint(moves, boundary = "fix")), nrow(fixed))
  robust <- parameter_inference(moves, boundary = "fix", vcov_type = "robust")
  expect_true(all(is.finite(robust$standard_error)))
  # Tables hold the bound as profile fits' tables do (1.3.5; NA in 1.3.4).
  profiles <- get_results(moves)
  expect_true(all(is.finite(profiles$mean_standard_error)))
  expect_equal(profiles$mean_standard_error,
               subset(fixed, parameter == "mean")$standard_error)
  # The bootstrap is the other answer to the same fit.
  boot <- quietly(parameter_inference(moves, course_engagement,
                                      method = "bootstrap", iter = 5,
                                      n_starts = 1, seed = 3))
  initial <- subset(boot, parameter == "initial_probability")
  expect_identical(nrow(initial), 6L)
  expect_true(all(initial$conf_low >= 0 & initial$conf_high <= 1))
})

test_that("`fix` holds a bound response probability and its complement row", {
  skip_if_not_installed("numDeriv")
  set.seed(4)
  data <- data.frame(person = rep(seq_len(80), each = 3), wave = rep(1:3, 80))
  data$a <- stats::rnorm(240)
  data$q <- sample(c("lo", "mid", "hi"), 240, TRUE)
  fit <- quietly(lta(data, c("a", "q"), "person", 2L, time = "wave",
                     categorical = "q", n_starts = 2, seed = 1),
                 c(.multilpa_expected_warnings, "latents_boundary"))
  held <- .lta_transition_bound(fit)
  expect_true("response_logit[2,q,hi]" %in% held$unconstrained)
  fixed <- parameter_inference(fit, boundary = "fix")
  expect_true("measurement.response.profile_2.q:hi" %in% attr(fixed, "fixed_at_bound"))
  expect_equal(sqrt(diag(attr(fixed, "covariance_unconstrained"))),
               .lta_conditional_errors(fit, held$unconstrained),
               tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("the held coordinates follow the profile models' rule", {
  data <- .lta_inference_fixture()
  fit <- .lta_inference_fit(data)
  # Nothing on a bound: nothing held, and the two values agree.
  expect_identical(.lta_transition_bound(fit),
                   list(natural = character(), unconstrained = character()))
  expect_identical(parameter_inference(fit, boundary = "fix"),
                   parameter_inference(fit))
  expect_identical(attr(parameter_inference(fit), "fixed_at_bound"), character())
  expect_identical(vcov(fit, boundary = "fix"), vcov(fit))
  # A floored first entry holds its own logit.
  first <- fit
  first$transition_probabilities[1L, , 1L] <- c(1e-10, 1 - 1e-10)
  expect_identical(.lta_transition_bound(first)$unconstrained,
                   "transition_logit[1,1,1]")
  # Its complement is on the upper bound and leaves the table with it.
  expect_identical(.lta_transition_bound(first)$natural,
                   c("transition_probability[1,1,1]", "transition_probability[1,1,2]"))
  # A floored reference (last) entry holds the largest logit only.
  three <- fit
  three$n_profiles <- 3L
  three$transition_probabilities <- array(1 / 3, c(3L, 3L, 1L))
  three$transition_probabilities[2L, , 1L] <- c(0.7, 0.3 - 1e-10, 1e-10)
  three$initial_probabilities <- matrix(1 / 3, 1L, 3L)
  three$response_probabilities <- NULL
  expect_identical(.lta_transition_bound(three),
                   list(natural = "transition_probability[1,2,3]",
                        unconstrained = "transition_logit[1,2,1]"))
})

test_that("bounds `fix` cannot hold are refused under either value", {
  data <- .lta_inference_fixture()
  fit <- .lta_inference_fit(data)
  variance <- fit
  variance$variances[1L, 1L] <- fit$min_variance
  expect_error(parameter_inference(variance, boundary = "fix"),
               class = "latents_boundary_fit")
  # The general engine (here, covariate transitions) holds nothing.
  data$z <- sin(seq_len(nrow(data)))
  general <- quietly(lta(data, c("a", "b"), "person", 2L, time = "wave",
                         transition_covariates = "z", n_starts = 3, seed = 1))
  expect_s3_class(general, "multilpa_lta")
  expect_identical(parameter_inference(general, boundary = "fix"),
                   parameter_inference(general))
  general$boundary <- TRUE
  expect_error(parameter_inference(general), class = "latents_boundary_fit")
  expect_error(parameter_inference(general, boundary = "fix"),
               class = "latents_boundary_fit")
})
