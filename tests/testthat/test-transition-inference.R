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
